#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/interp.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/timer.hxx"

#include <cctk.h>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <optional>
#include <unordered_map>
#include <utility>

namespace CapyrX::MultiPatch {

struct Location {
  int patch{0};
  int level{0};
  int index{0};
  int component{0};
};

using PointList = std::array<std::vector<CCTK_REAL>, dim>;

} // namespace CapyrX::MultiPatch

// See https://stackoverflow.com/a/2595226
static constexpr inline auto hash_combine(std::size_t h1, std::size_t h2)
    -> std::size_t {
  return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}

namespace std {
using namespace CapyrX;

template <> struct equal_to<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::equal_to<std::array<int, 4> >()(
        std::array<int, 4>{x.patch, x.level, x.index, x.component},
        std::array<int, 4>{y.patch, y.level, y.index, y.component});
  }
};

template <> struct hash<MultiPatch::Location> {
  std::size_t operator()(const MultiPatch::Location &x) const {
    return hash_combine(hash_combine(hash_combine(std::hash<int>()(x.patch),
                                                  std::hash<int>()(x.level)),
                                     std::hash<int>()(x.index)),
                        std::hash<int>()(x.component));
  }
};

} // namespace std

namespace CapyrX::MultiPatch {

struct ComponentSlice {
  std::size_t offset;
  std::size_t length;
};

struct InterpolationCache {
  // Epoch at which this cache was built. -1 = never built.
  CCTK_INT epoch{-1};

  // Ghost-zone source points in deterministic insertion order.
  std::vector<std::pair<Location, PointList> > ordered_components;
  // O(1) slot lookup used during the parallel fill phase.
  std::unordered_map<Location, std::size_t> component_index;

  // Flat coordinate arrays fed to InterpolationSetup.
  PointList coords;

  // The CarpetX interpolation setup (particle container + distribution).
  std::optional<CarpetX::InterpolationSetup> setup;

  // Per-patch outer-boundary policy for InterpolateFromSetup.
  std::vector<Arith::vect<Arith::vect<bool, 3>, 2> > policy;

  // Per-Location (offset, length) into the flat coords/results arrays.
  // Built once per epoch alongside ordered_components; avoids scatter copy.
  std::unordered_map<Location, ComponentSlice> slices;

  // Result buffer reused across calls; resized only on epoch change.
  std::vector<std::vector<CCTK_REAL> > results;
};

static InterpolationCache g_interp_cache;

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
#ifdef __CUDACC__
  const nvtxRangeId_t range =
      nvtxRangeStartA("CapyrX::MultiPatch1_Interpolate");
#endif

  static CarpetX::Timer timer("CapyrX::MultiPatch1_Interpolate");
  CarpetX::Interval interval(timer);

  // Step 0: Check input

  // Function Input checking
  if (cctkGH_ == nullptr) {
    CCTK_VERROR("The cctkGH_ pointer is null. Unable to continue");
  }

  if (nvars_ < 0) {
    CCTK_VERROR("The nvars_ variable is negatie. Unable to continue");
  }

  if (varinds_ == nullptr) {
    CCTK_VERROR("The varinds_ pointer is null. Unable to continue");
  }

  // Cast GH and wrap varinds
  const auto cctkGH{static_cast<const cGH *>(cctkGH_)};
  const std::vector<CCTK_INT> varinds(varinds_, varinds_ + nvars_);

  // Check input varinds validity
  for (const auto &varind : varinds) {
    if (varind < 0) {
      CCTK_VERROR("The varind %i is negative", varind);
    }

    const auto gi{CCTK_GroupIndexFromVarI(varind)};

    if (gi < 0) {
      CCTK_VERROR("The goup index %i for varind %i is negative", gi, varind);
    }

    const auto v0{CCTK_FirstVarIndexI(gi)};

    if (gi < 0) {
      CCTK_VERROR("The first var index %i for varind %i is negative", v0,
                  varind);
    }

    const auto vi{varind - v0};

    if (vi < 0) {
      CCTK_VERROR("The index %i for varind %i is negative", vi, varind);
    }

    cGroup group_data{};
    const auto result{CCTK_GroupData(gi, &group_data)};

    switch (result) {
    case -1:
      CCTK_VERROR("Error while retrieving group data for group index %i in "
                  "varind %i: The "
                  "group index is invalid",
                  gi, varind);
      break;

    case -2:
      CCTK_VERROR("Error while retrieving group data for group index %i in "
                  "varind %i: The group data buffer is null",
                  gi, varind);
      break;

    default:
      break;
    }

    if (group_data.grouptype != CCTK_GF) {
      CCTK_VERROR("The group data member \"grouptype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_GF\"",
                  gi, varind);
    }

    if (group_data.vartype != CCTK_VARIABLE_REAL) {
      CCTK_VERROR("The group data member \"vartype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_VARIABLE_REAL\"",
                  gi, varind);
    }

    if (group_data.disttype != CCTK_DISTRIB_DEFAULT) {
      CCTK_VERROR("The group data member \"disttype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_DISTRIB_DEFAULT\"",
                  gi, varind);
    }

    if (group_data.dim != dim) {
      CCTK_VERROR("The group data member \"dim\" in group index %i in varinds "
                  "%i is not %lu",
                  gi, varind, dim);
    }

    // TODO: Check centerings table
  }

  // Step 1: Rebuild the interpolation cache if the AMR epoch has changed.

  const CCTK_INT current_epoch = GetEpoch();

  if (current_epoch != g_interp_cache.epoch) {
    CCTK_VINFO("Interpolation cache out of date (cache epoch = %d. Current "
               "epoch = %d). Rebuilding",
               g_interp_cache.epoch, current_epoch);

    // Collect ghost-zone coordinates.
    // Serial pass: assign each component a fixed slot so the parallel fill
    // phase can write without any synchronisation.
    g_interp_cache.ordered_components.clear();
    g_interp_cache.component_index.clear();
    CarpetX::active_levels_t().loop_serially(
        [&](int patch, int level, int index, int component, const cGH *) {
          const Location location{patch, level, index, component};
          g_interp_cache.component_index[location] =
              g_interp_cache.ordered_components.size();
          g_interp_cache.ordered_components.emplace_back(location, PointList{});
        });

    // Parallel pass: each component writes to its pre-assigned slot.
    {
      CarpetX::active_levels_t().loop_parallel([&](int patch, int level,
                                                   int index, int component,
                                                   const cGH *cctkGH) {
        const Loop::GridDescBase grid(cctkGH);
        const std::array<int, dim> centering{0, 0, 0};
        const Loop::GF3D2layout layout(cctkGH, centering);

        const std::array<Loop::GF3D2<const CCTK_REAL>, dim> vcoords{
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(
                            CCTK_VarDataPtr(cctkGH, 0, "CoordinatesX::vcoordx"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(
                            CCTK_VarDataPtr(cctkGH, 0, "CoordinatesX::vcoordy"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordz")))};

        const auto &current_patch{g_patch_system->patches.at(patch)};
        const auto &patch_faces{current_patch.faces};

        const Location location{patch, level, index, component};
        const std::size_t slot = g_interp_cache.component_index.at(location);

        PointList source_points;

        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          // Skip outer boundaries
          for (int d = 0; d < dim; ++d) {
            if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {
              return;
            }

            if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {
              return;
            }
          }

          for (int d = 0; d < dim; ++d) {
            source_points[d].push_back(vcoords[d](p.I));
          }
        });
        g_interp_cache.ordered_components[slot].second =
            std::move(source_points);
      });
    }

    // Flatten into coordinate arrays and record per-Location slices.
    g_interp_cache.coords = {};
    g_interp_cache.slices.clear();
    std::size_t flat_offset = 0;
    for (const auto &[location, point_list] :
         g_interp_cache.ordered_components) {
      const std::size_t length = point_list[0].size();
      g_interp_cache.slices[location] = {flat_offset, length};
      flat_offset += length;
      for (int d = 0; d < dim; ++d) {
        g_interp_cache.coords[d].insert(g_interp_cache.coords[d].end(),
                                        point_list[d].begin(),
                                        point_list[d].end());
      }
    }

    assert(g_interp_cache.coords[0].size() == g_interp_cache.coords[1].size() &&
           g_interp_cache.coords[0].size() == g_interp_cache.coords[2].size() &&
           g_interp_cache.coords[1].size() == g_interp_cache.coords[2].size());

    // Build InterpolationSetup (expensive — skipped when epoch is unchanged)
    const std::size_t npoints_cache = g_interp_cache.coords[0].size();
    for (auto &r : g_interp_cache.results)
      r.resize(npoints_cache);
    g_interp_cache.setup.emplace(cctkGH, static_cast<CCTK_INT>(npoints_cache),
                                 g_interp_cache.coords[0].data(),
                                 g_interp_cache.coords[1].data(),
                                 g_interp_cache.coords[2].data());

    // Build per-patch outer-boundary policy
    const int npatches = cctkGH->cctk_npatches;
    g_interp_cache.policy.resize(npatches);

    static const bool have_boundary_spec =
        CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification2");

    if (have_boundary_spec) {
      std::array<CCTK_INT, 2 * dim> spec;
      for (int p = 0; p < npatches; ++p) {
        MultiPatch_GetBoundarySpecification2(p, 2 * dim, spec.data());
        for (int f = 0; f < 2; ++f)
          for (int d = 0; d < dim; ++d)
            g_interp_cache.policy[p][f][d] = !spec[2 * d + f];
      }
    } else {
      for (int p = 0; p < npatches; ++p)
        for (int f = 0; f < 2; ++f)
          for (int d = 0; d < dim; ++d)
            g_interp_cache.policy[p][f][d] = true;
    }

    g_interp_cache.epoch = current_epoch;
  }

  // Step 2: Interpolate using the cached setup.

  const std::size_t nvars = varinds.size();
  const std::size_t npoints = g_interp_cache.coords[0].size();

  const std::vector<CCTK_INT> operations(nvars, 0);

  auto &results = g_interp_cache.results;
  results.resize(nvars);
  std::vector<CCTK_REAL *> resultptrs(nvars);

  for (size_t n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }

  g_interp_cache.setup.value().Interpolate(
      cctkGH, nvars, varinds.data(), operations.data(), g_interp_cache.policy,
      resultptrs.data());

  // Diagnostic: count NaN values in the interpolated results. Non-zero means
  // source data in neighboring patches contains NaN (e.g. their interior cells
  // were not yet filled before MultiPatch_Interpolate ran, or
  // poison_undefined_values=yes is leaving interior ghosts uninitialized).
  {
    static bool g_nan_result_reported = false;
    if (!g_nan_result_reported && nvars > 0 && npoints > 0) {
      std::size_t n_nan = 0;
      for (std::size_t n = 0; n < nvars; ++n)
        for (std::size_t p = 0; p < npoints; ++p)
          if (std::isnan(results[n][p])) ++n_nan;
      if (n_nan > 0) {
        g_nan_result_reported = true;
        CCTK_VINFO("MultiPatch1_Interpolate: %zu NaN values in interpolated "
                   "results (out of %zu points × %zu vars = %zu total). "
                   "Source patches may have poisoned interior/ghost cells. "
                   "Angular ghost cells filled from NaN will leave corner "
                   "Neumann sources invalid even after MultiPatch.",
                   n_nan, npoints, nvars, npoints * nvars);
      } else {
        CCTK_VINFO("MultiPatch1_Interpolate: all %zu interpolated values "
                   "are non-NaN (%zu pts × %zu vars)",
                   npoints * nvars, npoints, nvars);
        g_nan_result_reported = true; // suppress future reports once clean
      }
    }
  }

  // Step 3: Write back results
  {
    CarpetX::active_levels_t().loop_parallel([&](int patch, int level,
                                                 int index, int component,
                                                 const cGH *cctkGH) {
      const Loop::GridDescBase grid(cctkGH);
      const std::array<int, dim> centering{0, 0, 0};
      const Loop::GF3D2layout layout(cctkGH, centering);

      std::vector<Loop::GF3D2<CCTK_REAL> > vars;
      vars.reserve(nvars);
      for (const auto &varind : varinds) {
        vars.emplace_back(layout, static_cast<CCTK_REAL *>(
                                      CCTK_VarDataPtrI(cctkGH, 0, varind)));
      }

      const auto &current_patch{g_patch_system->patches.at(patch)};
      const auto &patch_faces{current_patch.faces};

      const Location location{patch, level, index, component};
      const ComponentSlice &slice = g_interp_cache.slices.at(location);

      // Count ghost cells skipped due to outer boundary (includes pure outer
      // cells and corner cells at the outer+interpatch intersection).
      // Corner cells are NOT filled here and require a 2nd BC pass in
      // SyncGroupsByDirI after this function returns.
      int n_outer_skipped = 0;
      for (std::size_t n = 0; n < nvars; n++) {
        std::size_t pos = 0;

        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          // Skip outer boundaries (pure outer ghost cells and corner cells at
          // the outer+interpatch face intersection). Corner cells must be
          // filled by the 2nd BC pass in SyncGroupsByDirI.
          for (int d = 0; d < dim; ++d) {
            if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {
              if (n == 0) ++n_outer_skipped;
              return;
            }

            if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {
              if (n == 0) ++n_outer_skipped;
              return;
            }
          }

          vars[n](p.I) = results[n][slice.offset + pos];

          pos++;
        });
        assert(pos == slice.length);
      }
      // Report once globally: a non-zero count confirms corner cells exist
      // and were left unfilled (outer+interpatch intersection).
      {
        static bool g_skip_reported = false;
        if (n_outer_skipped > 0 && component == 0 && !g_skip_reported) {
          g_skip_reported = true;
          CCTK_VINFO("MultiPatch1_Interpolate [patch %d level %d]: "
                     "skipped %d outer-boundary ghost cells including "
                     "corner cells at outer+interpatch face intersections; "
                     "these are filled by the 2nd BC pass in SyncGroupsByDirI",
                     patch, level, n_outer_skipped);
        }
      }
    });
  }

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace CapyrX::MultiPatch
