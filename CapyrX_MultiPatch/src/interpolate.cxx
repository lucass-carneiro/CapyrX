#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/interp.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/timer.hxx"

#include <cctk.h>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <map>
#include <optional>

namespace CapyrX::MultiPatch {

struct Location {
  int patch{0};
  int level{0};
  int index{0};
  int component{0};
};

struct point_location {
  int patch{0};
  int level{0};
  int component{0};
  int i{0};
  int j{0};
  int k{0};
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

template <> struct less<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::less<std::array<int, 4> >()(
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

struct InterpolationCache {
  // Epoch at which this cache was built. 0 = never built.
  CCTK_INT epoch{0};

  // Ghost-zone coordinate mapping: Location -> (x,y,z) point lists.
  std::map<Location, PointList> source_mapping;

  // Flat coordinate arrays fed to InterpolationSetup.
  PointList coords;

  // The CarpetX interpolation setup (particle container + distribution).
  std::optional<CarpetX::InterpolationSetup> setup;

  // Per-patch outer-boundary policy for InterpolateFromSetup.
  std::vector<Arith::vect<Arith::vect<bool, 3>, 2> > policy;
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

    // Collect ghost-zone coordinates
    g_interp_cache.source_mapping.clear();

    CarpetX::active_levels_t().loop_parallel([&](int patch, int level,
                                                 int index, int component,
                                                 const cGH *cctkGH) {
      const Loop::GridDescBase grid(cctkGH);
      const std::array<int, dim> centering{0, 0, 0};
      const Loop::GF3D2layout layout(cctkGH, centering);

      const auto &current_patch{g_patch_system->patches.at(patch)};
      const auto &patch_faces{current_patch.faces};

      const Location location{patch, level, index, component};

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
#pragma omp critical
      {
        g_interp_cache.source_mapping[location] = std::move(source_points);
      }
    });

    // Flatten into coordinate arrays
    g_interp_cache.coords = {};
    for (const auto &[location, point_list] : g_interp_cache.source_mapping) {
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

  std::vector<std::vector<CCTK_REAL> > results(nvars);
  std::vector<CCTK_REAL *> resultptrs(nvars);

  for (size_t n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }

  CarpetX::InterpolateFromSetup(*g_interp_cache.setup, nvars, varinds.data(),
                                operations.data(), g_interp_cache.policy,
                                resultptrs.data());

  // Scatter interpolated values
  std::map<Location, std::vector<std::vector<CCTK_REAL> > > result_mapping;

  std::size_t pos = 0;
  for (const auto &[location, source_points] : g_interp_cache.source_mapping) {
    const std::size_t length = source_points[0].size();
    std::vector<std::vector<CCTK_REAL> > result_values(nvars);
    for (std::size_t n = 0; n < nvars; ++n)
      result_values.at(n).insert(result_values.at(n).begin(),
                                 results.at(n).data() + pos,
                                 results.at(n).data() + pos + length);
    pos += length;

#pragma omp critical
    {
      result_mapping[location] = std::move(result_values);
    }
  }
  assert(pos == results.at(0).size());

  // Step 3: Write back results
  CarpetX::active_levels_t().loop_parallel([&](int patch, int level, int index,
                                               int component,
                                               const cGH *cctkGH) {
    const Loop::GridDescBase grid(cctkGH);
    const std::array<int, dim> centering{0, 0, 0};
    const Loop::GF3D2layout layout(cctkGH, centering);

    const auto &current_patch{g_patch_system->patches.at(patch)};
    const auto &patch_faces{current_patch.faces};

    const Location location{patch, level, index, component};
    const std::vector<std::vector<CCTK_REAL> > &result_values =
        result_mapping.at(location);

    for (std::size_t n = 0; n < nvars; ++n) {
      const std::vector<CCTK_REAL> &result_values_n = result_values.at(n);

      const Loop::GF3D2<CCTK_REAL> var(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, varinds.at(n))));

      std::size_t pos = 0;

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

        var(p.I) = result_values_n[pos];

        pos++;
      });
      assert(pos == result_values.at(0).size());
    }
  });

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace CapyrX::MultiPatch
