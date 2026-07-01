#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/interp.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/timer.hxx"

#include <cctk.h>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <cmath>
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

  // CarpetX's interpolation order. Needed to decide how many radial cells
  // near a patch's own outer-radial face cannot be served by any donor's
  // interpolation stencil once outer faces are excluded from `policy` (see
  // mp_corners_2.md, fix A companion change).
  CCTK_INT interpolation_order{1};

  // Per-Location (offset, length) into the flat coords/results arrays.
  // Built once per epoch alongside ordered_components; avoids scatter copy.
  std::unordered_map<Location, ComponentSlice> slices;

  // Result buffer reused across calls; resized only on epoch change.
  std::vector<std::vector<CCTK_REAL> > results;
};

static InterpolationCache g_interp_cache;

// Fix A companion check (mp_corners_2.md): once outer faces are excluded
// from the interpolation policy below, a target ghost point that is a
// genuine interpatch (seam) ghost in some direction, but sits within
// `nghostzones + interpolation_order` cells of a patch's own outer-radial
// face, can no longer be served by any donor: the donor's radial stencil
// would have to anchor inside its own outer-BC-filled ghost zone, which is
// now disallowed, and the interpolator's anchor placement only ever nudges
// by a single cell. Such points are excluded from MultiPatch1_Interpolate
// below and instead filled by the same-patch fallback extrapolation.
static bool is_within_donor_outer_reach(const Patch &patch,
                                        const Loop::GridDescBase &grid,
                                        const Loop::PointDesc &p,
                                        const CCTK_INT interpolation_order) {
  if (!patch.c_is_radial)
    return false;

  constexpr int d = 2; // c, the patch-local radial direction

  if (p.NI[d] != 0)
    return false; // already excluded as a literal outer-boundary ghost cell

  for (int f = 0; f < 2; ++f) {
    if (!patch.faces[f][d].is_outer_boundary)
      continue;
    if (!grid.bbox[f][d])
      continue; // this component does not touch the outer face

    const int nghost = grid.nghostzones[d];
    const int dist = f == 0 ? p.I[d] - nghost
                            : (grid.lsh[d] - nghost - 1) - p.I[d];
    if (dist < nghost + interpolation_order)
      return true;
  }

  return false;
}

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

  // mp_corners_4.md, fix option 2: a freshly (re)built cache is the signal
  // that no patch's interpatch ghosts are guaranteed populated yet -- most
  // notably true on the very first call of a run, before any patch has ever
  // been written to by this function, but this signal also naturally
  // re-arms after every future regrid. Fix A's policy (below) assumes some
  // prior call already left every patch's interpatch ghosts valid; that
  // assumption is false exactly when the cache was just rebuilt. Recording
  // this here (rather than a dedicated one-shot "first call" flag) keeps
  // the bootstrap self-correcting instead of depending on bookkeeping that
  // could go stale under a future schedule reordering.
  const bool cache_was_rebuilt = current_epoch != g_interp_cache.epoch;

  if (cache_was_rebuilt) {
    CCTK_VINFO("Interpolation cache out of date (cache epoch = %d. Current "
               "epoch = %d). Rebuilding",
               g_interp_cache.epoch, current_epoch);

    // Fix A companion (mp_corners_2.md): read CarpetX's interpolation order
    // once per cache rebuild. Needed below (and in the write-back step) to
    // determine which near-outer seam ghosts cannot be served by any donor
    // once outer faces are excluded from the interpolation policy.
    {
      const auto *const interpolation_order_ptr =
          static_cast<const CCTK_INT *>(
              CCTK_ParameterGet("interpolation_order", "CarpetX", nullptr));
      if (interpolation_order_ptr == nullptr) {
        CCTK_ERROR(
            "Unable to read parameter interpolation_order from CarpetX");
      }
      g_interp_cache.interpolation_order = *interpolation_order_ptr;
    }

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
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordx"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordy"))),
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

          // Fix A companion (mp_corners_2.md): also skip near-outer-radial
          // seam ghosts that no donor can serve now that outer faces are
          // excluded from the interpolation policy. The fallback
          // extrapolation pass in the write-back step fills these instead.
          if (is_within_donor_outer_reach(
                  current_patch, grid, p,
                  g_interp_cache.interpolation_order)) {
            return;
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
            // Allow an interpolation stencil to anchor into the ghost
            // region only on genuine interpatch faces (filled with valid
            // donor data by MultiPatch_Interpolate itself). On outer faces
            // the ghost cells hold low-order BC-extrapolated values that
            // are not valid interior data; push the stencil inward there
            // instead, so interpatch results never depend on a
            // neighbour's outer BC (fix A, mp_corners_2.md).
            g_interp_cache.policy[p][f][d] = spec[2 * d + f];
      }
    } else {
      for (int p = 0; p < npatches; ++p)
        for (int f = 0; f < 2; ++f)
          for (int d = 0; d < dim; ++d)
            g_interp_cache.policy[p][f][d] = false;
    }

    g_interp_cache.epoch = current_epoch;
  }

  // Step 2/3: Interpolate using the cached setup, then write the results
  // back. Factored into a lambda, parameterized on the outer-boundary
  // policy, so it can be run twice in immediate succession on a freshly
  // rebuilt cache (mp_corners_4.md, fix option 2): once with the
  // conservative pre-Fix-A policy to seed every patch's interpatch ghosts
  // with finite (if less accurate) values, then again with Fix A's real
  // policy now that "some prior call already filled the ghosts" is true.

  const std::size_t nvars = varinds.size();
  const std::size_t npoints = g_interp_cache.coords[0].size();

  const std::vector<CCTK_INT> operations(nvars, 0);

  const auto run_interpolation_pass =
      [&](const std::vector<Arith::vect<Arith::vect<bool, 3>, 2> > &policy) {
        auto &results = g_interp_cache.results;
        results.resize(nvars);
        std::vector<CCTK_REAL *> resultptrs(nvars);

        for (size_t n = 0; n < nvars; ++n) {
          results.at(n).resize(npoints);
          resultptrs.at(n) = results.at(n).data();
        }

        g_interp_cache.setup.value().Interpolate(cctkGH, nvars,
                                                   varinds.data(),
                                                   operations.data(), policy,
                                                   resultptrs.data());

// Diagnostic: count NaN values in the interpolated results. Non-zero means
// source data in neighboring patches contains NaN (e.g. their interior cells
// were not yet filled before MultiPatch_Interpolate ran, or
// poison_undefined_values=yes is leaving interior ghosts uninitialized).
#ifdef CCTK_DEBUG
        {
          static bool g_nan_result_reported = false;
          if (!g_nan_result_reported && nvars > 0 && npoints > 0) {
            std::size_t n_nan = 0;
            for (std::size_t n = 0; n < nvars; ++n)
              for (std::size_t p = 0; p < npoints; ++p)
                if (std::isnan(results[n][p])) {
                  ++n_nan;
                }
            if (n_nan > 0) {
              g_nan_result_reported = true;
              CCTK_VINFO(
                  "MultiPatch1_Interpolate: %zu NaN values in interpolated "
                  "results (out of %zu points x %zu vars = %zu total). Source "
                  "patches may have poisoned interior/ghost cells. ",
                  n_nan, npoints, nvars, npoints * nvars);
            }
          }
        }
#endif // CCTK_DEBUG

        // Step 3: Write back results
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
#ifdef CCTK_DEBUG
        int n_outer_skipped = 0;
        int n_donor_reach_skipped = 0;
#endif // CCTK_DEBUG

        for (std::size_t n = 0; n < nvars; n++) {
          std::size_t pos = 0;

          // Note: This includes symmetry points
          grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
            // Skip outer boundaries (pure outer ghost cells and corner cells at
            // the outer+interpatch face intersection). Corner cells must be
            // filled by the 2nd BC pass in SyncGroupsByDirI.
            for (int d = 0; d < dim; ++d) {
              if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {

#ifdef CCTK_DEBUG
                if (n == 0) {
                  ++n_outer_skipped;
                }
#endif // CCTK_DEBUG

                return;
              }

              if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {

#ifdef CCTK_DEBUG
                if (n == 0) {
                  ++n_outer_skipped;
                }
#endif // CCTK_DEBUG

                return;
              }
            }

            // Fix A companion (mp_corners_2.md): skip near-outer-radial seam
            // ghosts that no donor can serve now that outer faces are
            // excluded from the interpolation policy; the fallback
            // extrapolation pass below fills them instead.
            if (is_within_donor_outer_reach(
                    current_patch, grid, p,
                    g_interp_cache.interpolation_order)) {

#ifdef CCTK_DEBUG
              if (n == 0) {
                ++n_donor_reach_skipped;
              }
#endif // CCTK_DEBUG

              return;
            }

            vars[n](p.I) = results[n][slice.offset + pos];

            pos++;
          });
          assert(pos == slice.length);
        }

// Report once globally: a non-zero count confirms corner cells exist
// and were left unfilled (outer+interpatch intersection).
#ifdef CCTK_DEBUG
        {
          if (n_outer_skipped > 0 && component == 0) {
            CCTK_VINFO("MultiPatch1_Interpolate [patch %d level %d]: "
                       "skipped %d outer-boundary ghost cells including "
                       "corner cells at outer+interpatch face intersections; "
                       "these are filled by the 2nd BC pass in SyncGroupsByDirI",
                       patch, level, n_outer_skipped);
          }
          if (n_donor_reach_skipped > 0 && component == 0) {
            CCTK_VINFO(
                "MultiPatch1_Interpolate [patch %d level %d]: skipped %d "
                "near-outer-radial seam ghost cells that no donor could serve "
                "under the flipped interpolation policy (fix A, "
                "mp_corners_2.md); these are filled by the same-patch "
                "fallback extrapolation below",
                patch, level, n_donor_reach_skipped);
          }
        }
#endif // CCTK_DEBUG

        // Step 3b (fix A companion, mp_corners_2.md): fill the near-outer-
        // radial seam ghosts skipped above via same-patch quadratic
        // extrapolation along the radial direction. Layers are filled one at
        // a time, marching outward from already-valid data (dist descending),
        // so each layer only ever reads strictly-more-interior data that is
        // either genuinely donor-interpolated above or was extrapolated by an
        // earlier iteration of this same loop. This keeps a neighbour's
        // BC-contaminated outer-radial ghost data out of the interpatch
        // channel entirely, while still providing a finite, continuous value
        // for these seam cells.
        if (current_patch.c_is_radial) {
          constexpr int d_rad = 2;
          const int nghost_r = grid.nghostzones[d_rad];
          const CCTK_INT order = g_interp_cache.interpolation_order;

          for (int f_rad = 0; f_rad < 2; ++f_rad) {
            if (!patch_faces[f_rad][d_rad].is_outer_boundary)
              continue;
            if (!grid.bbox[f_rad][d_rad])
              continue;

            const int threshold = nghost_r + order;
            const int last_interior = grid.lsh[d_rad] - nghost_r - 1;
            const int sign = f_rad == 0 ? +1 : -1;

            if (grid.lsh[d_rad] < 2 * nghost_r + order + 3) {
              CCTK_VERROR(
                  "MultiPatch1_Interpolate [patch %d level %d]: radial "
                  "extent (%d cells) is too small to fall back to "
                  "same-patch extrapolation for near-outer seam ghosts "
                  "(need at least %d cells); increase the number of radial "
                  "cells or decrease interpolation_order.",
                  patch, level, grid.lsh[d_rad], 2 * nghost_r + order + 3);
            }

            for (int k = threshold - 1; k >= 0; --k) {
              const int idx_rad = f_rad == 0 ? nghost_r + k : last_interior - k;

              grid.loop_bnd<0, 0, 0>(
                  grid.nghostzones, [&](const Loop::PointDesc &p) {
                    if (p.I[d_rad] != idx_rad || p.NI[d_rad] != 0)
                      return;

                    // Stay clear of cells that are themselves a literal
                    // outer-boundary ghost in some other direction; those
                    // remain the responsibility of the true outer BC / 2nd
                    // BC pass.
                    for (int dd = 0; dd < dim; ++dd) {
                      if (p.NI[dd] < 0 && patch_faces[0][dd].is_outer_boundary)
                        return;
                      if (p.NI[dd] > 0 && patch_faces[1][dd].is_outer_boundary)
                        return;
                    }

                    auto nbr = p.I;
                    for (std::size_t n = 0; n < nvars; ++n) {
                      nbr[d_rad] = idx_rad + 1 * sign;
                      const CCTK_REAL y1 = vars[n](nbr);
                      nbr[d_rad] = idx_rad + 2 * sign;
                      const CCTK_REAL y2 = vars[n](nbr);
                      nbr[d_rad] = idx_rad + 3 * sign;
                      const CCTK_REAL y3 = vars[n](nbr);

                      // Diagnostic (mp_corners_3.md investigation): this
                      // fallback assumes y1/y2/y3 are already finite -- either
                      // genuinely donor-interpolated (Step 1/3 above, for
                      // radial distances outside the near-outer band) or
                      // extrapolated by an earlier, more-interior iteration of
                      // this same k-loop. If that assumption is violated, this
                      // formula silently manufactures and then propagates NaN
                      // outward through every remaining layer and into the
                      // true outer ghost zone. Report the first few
                      // occurrences with full location context so the actual
                      // upstream cause (bad donor data vs. a coverage gap in
                      // the skip/fallback bookkeeping) can be identified
                      // without re-deriving it from the much-larger poison
                      // footprint that SyncGroupsByDirI reports downstream.
                      if (!(std::isfinite(y1) && std::isfinite(y2) &&
                            std::isfinite(y3))) {
                        static int n_reports = 0;
#pragma omp critical
                        if (n_reports < 20) {
                          ++n_reports;
                          CCTK_VWARN(
                              CCTK_WARN_ALERT,
                              "MultiPatch1_Interpolate same-patch fallback "
                              "(fix A companion, mp_corners_2.md/mp_corners_3.md): "
                              "non-finite seed for var '%s' at patch %d level "
                              "%d component %d, target index [%d,%d,%d] "
                              "(f_rad=%d, idx_rad=%d, threshold-relative "
                              "layer=%d), seeds at radial offsets 1,2,3 from "
                              "target = %.17g, %.17g, %.17g. The seed must be "
                              "finite before this fallback runs; a non-finite "
                              "seed here means the bug is upstream (the "
                              "normal cross-patch donor interpolation for "
                              "this lateral ghost cell, or a gap between the "
                              "is_within_donor_outer_reach skip and this "
                              "fallback's coverage), not in the extrapolation "
                              "formula itself.",
                              CCTK_VarName(varinds[n]), patch, level, component,
                              p.I[0], p.I[1], p.I[2], f_rad, idx_rad,
                              threshold - 1 - k, y1, y2, y3);
                        }
                      }

                      // Quadratic extrapolation (the same formula fix B
                      // recommends for the true outer BC), sourced entirely
                      // from this patch's own already-valid interior or
                      // previously-extrapolated data -- never from a
                      // neighbour's outer-BC ghosts.
                      vars[n](p.I) = 3 * y1 - 3 * y2 + y3;
                    }
                  });
            }
          }
        }
        });
  };

  if (cache_was_rebuilt) {
    // mp_corners_4.md, fix option 2: bootstrap pass with the conservative,
    // pre-Fix-A policy (never anchor at any patch face), so every
    // interpatch ghost the accurate pass below could possibly read is
    // already finite, however inaccurate. This is what makes Fix A's "some
    // prior call already filled the ghosts" assumption true on this call.
    CCTK_VINFO("MultiPatch1_Interpolate: cache was just (re)built; running a "
               "one-time conservative bootstrap interpolation pass before "
               "the accurate pass (mp_corners_4.md, fix option 2)");
    const std::vector<Arith::vect<Arith::vect<bool, 3>, 2> >
        conservative_policy(g_interp_cache.policy.size());
    run_interpolation_pass(conservative_policy);
  }

  run_interpolation_pass(g_interp_cache.policy);

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace CapyrX::MultiPatch
