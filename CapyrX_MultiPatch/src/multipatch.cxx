#include "multipatch.hxx"

#include "patch_systems.hxx"
#include "type_aliases.hxx"

#include "cartesian/cartesian.hxx"
#include "cubed_sphere/cubed_sphere.hxx"

#include "../../../CarpetX/CarpetX/src/timer.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <cmath>

namespace CapyrX::MultiPatch {

std::unique_ptr<PatchSystem> g_patch_system{nullptr};

extern "C" CCTK_INT
MultiPatch1_GetSystemSpecification(CCTK_INT *restrict const npatches) {
  *npatches = g_patch_system->patches.size();
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetPatchSpecification(
    const CCTK_INT ipatch, CCTK_INT *restrict const is_cartesian,
    const CCTK_INT size, CCTK_INT *restrict const ncells,
    CCTK_REAL *restrict const xmin, CCTK_REAL *restrict const xmax) {

  assert(ipatch >= 0 && ipatch < g_patch_system->patches.size());
  assert(size == dim);

  const auto &patch{g_patch_system->patches.at(ipatch)};

  if (is_cartesian != nullptr) {
    *is_cartesian = static_cast<CCTK_INT>(patch.is_cartesian);
  }

  if (ncells != nullptr && xmin != nullptr && xmax != nullptr) {
    for (int d = 0; d < dim; ++d) {
      ncells[d] = patch.ncells[d];
      xmin[d] = patch.xmin[d];
      xmax[d] = patch.xmax[d];
    }
  }

  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetBoundarySpecification2(
    const CCTK_INT ipatch, const CCTK_INT size,
    CCTK_INT *restrict const is_interpatch_boundary) {

  assert(ipatch >= 0 && ipatch < g_patch_system->patches.size());
  assert(size == 2 * dim);

  const auto &patch{g_patch_system->patches.at(ipatch)};

  for (int f = 0; f < 2; ++f) {
    for (int d = 0; d < dim; ++d) {
      is_interpatch_boundary[2 * d + f] = !patch.faces[f][d].is_outer_boundary;
    }
  }

  return 0;
}

template <CapyrX::MultiPatch::PatchSystems sys>
static inline void global2local_kernel(const CCTK_INT npoints,
                                       const CCTK_REAL *restrict const globalsx,
                                       const CCTK_REAL *restrict const globalsy,
                                       const CCTK_REAL *restrict const globalsz,
                                       CCTK_INT *restrict const patches,
                                       CCTK_REAL *restrict const localsx,
                                       CCTK_REAL *restrict const localsy,
                                       CCTK_REAL *restrict const localsz) {

  using namespace CapyrX::MultiPatch;

  DECLARE_CCTK_PARAMETERS;

  for (int n = 0; n < npoints; ++n) {
    const svec_t global_coords{globalsx[n], globalsy[n], globalsz[n]};

    std_tuple<int, svec_t> g2l{};

    if constexpr (sys == PatchSystems::cartesian) {
      const Cartesian::PatchParams p{
          .ncells_i = cartesian_ncells_i,
          .ncells_j = cartesian_ncells_j,
          .ncells_k = cartesian_ncells_k,

          .xmin = cartesian_xmin,
          .ymin = cartesian_xmin,
          .zmin = cartesian_xmin,

          .xmax = cartesian_xmax,
          .ymax = cartesian_xmax,
          .zmax = cartesian_xmax,
      };
      g2l = Cartesian::global2local(p, global_coords);
    } else if constexpr (sys == PatchSystems::cubed_spehre) {
      const CubedSphere::PatchParams p{
          .angular_cells = cubed_sphere_angular_cells + 2 * patch_overlap,
          .radial_cells = cubed_sphere_radial_cells + 2 * patch_overlap,
          .inner_boundary = cubed_sphere_inner_boundary_radius,
          .outer_boundary = cubed_sphere_outer_boundary_radius,
          .patch_overlap = patch_overlap};
      g2l = CubedSphere::global2local(p, global_coords);
    }

    const auto patch{std::get<0>(g2l)};
    const auto local_coords{std::get<1>(g2l)};
    patches[n] = patch;
    localsx[n] = local_coords(0);
    localsy[n] = local_coords(1);
    localsz[n] = local_coords(2);
  }
}

extern "C" void MultiPatch1_GlobalToLocal2(
    const CCTK_INT npoints, const CCTK_REAL *restrict const globalsx,
    const CCTK_REAL *restrict const globalsy,
    const CCTK_REAL *restrict const globalsz, CCTK_INT *restrict const patches,
    CCTK_REAL *restrict const localsx, CCTK_REAL *restrict const localsy,
    CCTK_REAL *restrict const localsz) {

#ifdef __CUDACC__
  const nvtxRangeId_t range =
      nvtxRangeStartA("CapyrX::MultiPatch1_GlobalToLocal2");
#endif

  static CarpetX::Timer timer("CapyrX::MultiPatch1_GlobalToLocal2");
  CarpetX::Interval interval(timer);

  DECLARE_CCTK_PARAMETERS;

  switch (g_patch_system->id_tag) {

  case PatchSystems::none:
    if (verbose) {
      CCTK_VWARN(
          CCTK_WARN_ALERT,
          "MultiPatch1_GlobalToLocal2 called with no active patch system");
    }
    break;

  case PatchSystems::cartesian: {
    global2local_kernel<PatchSystems::cartesian>(npoints, globalsx, globalsy,
                                                 globalsz, patches, localsx,
                                                 localsy, localsz);
    break;
  }

  case PatchSystems::cubed_spehre: {
    global2local_kernel<PatchSystems::cubed_spehre>(npoints, globalsx, globalsy,
                                                    globalsz, patches, localsx,
                                                    localsy, localsz);
    break;
  }

  default:
    CCTK_VERROR(
        "MultiPatch1_GlobalToLocal2 called for an unknown patch system");
    break;
  }

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

template <CapyrX::MultiPatch::PatchSystems sys>
static inline void local2global_kernel(const CCTK_INT npoints,
                                       const CCTK_INT *restrict const patches,
                                       const CCTK_REAL *restrict const localsx,
                                       const CCTK_REAL *restrict const localsy,
                                       const CCTK_REAL *restrict const localsz,
                                       CCTK_REAL *restrict const globalsx,
                                       CCTK_REAL *restrict const globalsy,
                                       CCTK_REAL *restrict const globalsz) {
  using namespace CapyrX::MultiPatch;

  DECLARE_CCTK_PARAMETERS;

  for (CCTK_INT i = 0; i < npoints; i++) {
    const auto patch{patches[i]};
    const svec_t local_coords{localsx[i], localsy[i], localsz[i]};

    svec_t global_vars{};

    if constexpr (sys == PatchSystems::cartesian) {
      const Cartesian::PatchParams p{
          .ncells_i = cartesian_ncells_i,
          .ncells_j = cartesian_ncells_j,
          .ncells_k = cartesian_ncells_k,

          .xmin = cartesian_xmin,
          .ymin = cartesian_xmin,
          .zmin = cartesian_xmin,

          .xmax = cartesian_xmax,
          .ymax = cartesian_xmax,
          .zmax = cartesian_xmax,
      };
      global_vars = Cartesian::local2global(p, patch, local_coords);
    } else if constexpr (sys == PatchSystems::cubed_spehre) {
      const CubedSphere::PatchParams p{
          .angular_cells = cubed_sphere_angular_cells + 2 * patch_overlap,
          .radial_cells = cubed_sphere_radial_cells + 2 * patch_overlap,
          .inner_boundary = cubed_sphere_inner_boundary_radius,
          .outer_boundary = cubed_sphere_outer_boundary_radius,
          .patch_overlap = patch_overlap};
      global_vars = CubedSphere::local2global(p, patch, local_coords);
    }

    globalsx[i] = global_vars(0);
    globalsy[i] = global_vars(1);
    globalsz[i] = global_vars(2);
  }
}

extern "C" void MultiPatch1_LocalToGlobal2(
    const CCTK_INT npoints, const CCTK_INT *restrict const patches,
    const CCTK_REAL *restrict const localsx,
    const CCTK_REAL *restrict const localsy,
    const CCTK_REAL *restrict const localsz, CCTK_REAL *restrict const globalsx,
    CCTK_REAL *restrict const globalsy, CCTK_REAL *restrict const globalsz) {

#ifdef __CUDACC__
  const nvtxRangeId_t range =
      nvtxRangeStartA("CapyrX::MultiPatch1_LocalToGlobal2");
#endif

  static CarpetX::Timer timer("CapyrX::MultiPatch1_LocalToGlobal2");
  CarpetX::Interval interval(timer);

  DECLARE_CCTK_PARAMETERS;

  switch (g_patch_system->id_tag) {

  case PatchSystems::none:
    if (verbose) {
      CCTK_VWARN(
          CCTK_WARN_ALERT,
          "MultiPatch1_LocalToGlobal2 called with no active patch system");
    }
    break;

  case PatchSystems::cartesian: {
    local2global_kernel<PatchSystems::cartesian>(npoints, patches, localsx,
                                                 localsy, localsz, globalsx,
                                                 globalsy, globalsz);
    break;
  }

  case PatchSystems::cubed_spehre: {
    local2global_kernel<PatchSystems::cubed_spehre>(npoints, patches, localsx,
                                                    localsy, localsz, globalsx,
                                                    globalsy, globalsz);
    break;
  }

  default:
    CCTK_VERROR(
        "MultiPatch1_LocalToGlobal2 called for an unknown patch system");
    break;
  }

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

extern "C" int CapyrX_MultiPatch_Setup() {
#ifdef __CUDACC__
  const nvtxRangeId_t range =
      nvtxRangeStartA("CapyrX::CapyrX_MultiPatch_Setup");
#endif

  static CarpetX::Timer timer("CapyrX::CapyrX_MultiPatch_Setup");
  CarpetX::Interval interval(timer);

  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(patch_system, "none")) {
    g_patch_system = nullptr;

  } else if (CCTK_EQUALS(patch_system, "Cartesian")) {
    Cartesian::PatchParams p{
        .ncells_i = cartesian_ncells_i,
        .ncells_j = cartesian_ncells_j,
        .ncells_k = cartesian_ncells_k,

        .xmin = cartesian_xmin,
        .ymin = cartesian_xmin,
        .zmin = cartesian_xmin,

        .xmax = cartesian_xmax,
        .ymax = cartesian_xmax,
        .zmax = cartesian_xmax,
    };

    g_patch_system = std::make_unique<PatchSystem>(Cartesian::make_system(p));

  } else if (CCTK_EQUALS(patch_system, "Cubed sphere")) {
    CubedSphere::PatchParams p{
        .angular_cells = cubed_sphere_angular_cells + 2 * patch_overlap,
        .radial_cells = cubed_sphere_radial_cells + 2 * patch_overlap,
        .inner_boundary = cubed_sphere_inner_boundary_radius,
        .outer_boundary = cubed_sphere_outer_boundary_radius,
        .patch_overlap = patch_overlap};

    g_patch_system = std::make_unique<PatchSystem>(CubedSphere::make_system(p));

  } else {
    CCTK_VERROR("Unbable to setup unknown patch system \"%s\"", patch_system);
  }

  if (verbose) {
    CCTK_VINFO("Using patch system \"%s\"", patch_system);
  }

  return 0;
}

#define COORDINATE_SETUP_KERNEL(d2local_dglobal2_func)                         \
  grid.loop_all_device<0, 0, 0>(                                               \
      grid.nghostzones, [=] CCTK_HOST CCTK_DEVICE(                             \
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE { \
        const svec_t a{p.x, p.y, p.z};                                         \
                                                                               \
        const auto d2J_tuple{d2local_dglobal2_func(par, p.patch, a)};          \
                                                                               \
        const auto &x{std::get<0>(d2J_tuple)};                                 \
        const auto &J{std::get<1>(d2J_tuple)};                                 \
        const auto &dJ{std::get<2>(d2J_tuple)};                                \
                                                                               \
        vcoordx(p.I) = x(0);                                                   \
        vcoordy(p.I) = x(1);                                                   \
        vcoordz(p.I) = x(2);                                                   \
                                                                               \
        vJ_da_dx(p.I) = J(0)(0);                                               \
        vJ_da_dy(p.I) = J(0)(1);                                               \
        vJ_da_dz(p.I) = J(0)(2);                                               \
        vJ_db_dx(p.I) = J(1)(0);                                               \
        vJ_db_dy(p.I) = J(1)(1);                                               \
        vJ_db_dz(p.I) = J(1)(2);                                               \
        vJ_dc_dx(p.I) = J(2)(0);                                               \
        vJ_dc_dy(p.I) = J(2)(1);                                               \
        vJ_dc_dz(p.I) = J(2)(2);                                               \
                                                                               \
        vdJ_d2a_dxdx(p.I) = dJ(0)(0, 0);                                       \
        vdJ_d2a_dxdy(p.I) = dJ(0)(0, 1);                                       \
        vdJ_d2a_dxdz(p.I) = dJ(0)(0, 2);                                       \
        vdJ_d2a_dydy(p.I) = dJ(0)(1, 1);                                       \
        vdJ_d2a_dydz(p.I) = dJ(0)(1, 2);                                       \
        vdJ_d2a_dzdz(p.I) = dJ(0)(2, 2);                                       \
        vdJ_d2b_dxdx(p.I) = dJ(1)(0, 0);                                       \
        vdJ_d2b_dxdy(p.I) = dJ(1)(0, 1);                                       \
        vdJ_d2b_dxdz(p.I) = dJ(1)(0, 2);                                       \
        vdJ_d2b_dydy(p.I) = dJ(1)(1, 1);                                       \
        vdJ_d2b_dydz(p.I) = dJ(1)(1, 2);                                       \
        vdJ_d2b_dzdz(p.I) = dJ(1)(2, 2);                                       \
        vdJ_d2c_dxdx(p.I) = dJ(2)(0, 0);                                       \
        vdJ_d2c_dxdy(p.I) = dJ(2)(0, 1);                                       \
        vdJ_d2c_dxdz(p.I) = dJ(2)(0, 2);                                       \
        vdJ_d2c_dydy(p.I) = dJ(2)(1, 1);                                       \
        vdJ_d2c_dydz(p.I) = dJ(2)(1, 2);                                       \
        vdJ_d2c_dzdz(p.I) = dJ(2)(2, 2);                                       \
      });                                                                      \
                                                                               \
  grid.loop_all_device<1, 1, 1>(                                               \
      grid.nghostzones, [=] CCTK_HOST CCTK_DEVICE(                             \
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE { \
        const svec_t a{p.x, p.y, p.z};                                         \
                                                                               \
        const auto d2J_tuple{d2local_dglobal2_func(par, p.patch, a)};          \
                                                                               \
        const auto &x{std::get<0>(d2J_tuple)};                                 \
        const auto &J{std::get<1>(d2J_tuple)};                                 \
        const auto &dJ{std::get<2>(d2J_tuple)};                                \
                                                                               \
        const auto detJ{J(0)(0) * (J(1)(1) * J(2)(2) - J(1)(2) * J(2)(1)) +    \
                        J(0)(1) * (J(1)(2) * J(2)(0) - J(1)(0) * J(2)(2)) +    \
                        J(0)(2) * (J(1)(0) * J(2)(1) - J(1)(1) * J(2)(0))};    \
                                                                               \
        ccoordx(p.I) = x(0);                                                   \
        ccoordy(p.I) = x(1);                                                   \
        ccoordz(p.I) = x(2);                                                   \
                                                                               \
        cvol(p.I) = (p.dx * p.dy * p.dz) * detJ;                               \
      })

extern "C" void CapyrX_MultiPatch_Coordinates_Setup(CCTK_ARGUMENTS) {
#ifdef __CUDACC__
  const nvtxRangeId_t range =
      nvtxRangeStartA("CapyrX::CapyrX_MultiPatch_Coordinates_Setup");
#endif

  DECLARE_CCTK_ARGUMENTSX_CapyrX_MultiPatch_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;

  switch (g_patch_system->id_tag) {

  case PatchSystems::none:
    break;

  case PatchSystems::cartesian: {
    Cartesian::PatchParams par{
        .ncells_i = cartesian_ncells_i,
        .ncells_j = cartesian_ncells_j,
        .ncells_k = cartesian_ncells_k,

        .xmin = cartesian_xmin,
        .ymin = cartesian_xmin,
        .zmin = cartesian_xmin,

        .xmax = cartesian_xmax,
        .ymax = cartesian_xmax,
        .zmax = cartesian_xmax,
    };
    COORDINATE_SETUP_KERNEL(Cartesian::d2local_dglobal2);
    break;
  }

  case PatchSystems::cubed_spehre: {
    CubedSphere::PatchParams par{
        .angular_cells = cubed_sphere_angular_cells + 2 * patch_overlap,
        .radial_cells = cubed_sphere_radial_cells + 2 * patch_overlap,
        .inner_boundary = cubed_sphere_inner_boundary_radius,
        .outer_boundary = cubed_sphere_outer_boundary_radius,
        .patch_overlap = patch_overlap};

    COORDINATE_SETUP_KERNEL(CubedSphere::d2local_dglobal2);
    break;
  }

  default:
    CCTK_VERROR("Unable to setup coordinates for unknown patch system");
    break;
  }

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

extern "C" void CapyrX_MultiPatch_Check_Parameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_MultiPatch_Check_Parameters;
  DECLARE_CCTK_PARAMETERS;

  const auto cubed_sphere_radius_condition{
      cubed_sphere_outer_boundary_radius >
          4.0 * cubed_sphere_inner_boundary_radius ||
      ((cubed_sphere_outer_boundary_radius -
        4.0 * cubed_sphere_inner_boundary_radius) < 1.0e-8)};

  if (!cubed_sphere_radius_condition) {
    CCTK_VERROR("When using the \"Cubed sphere\" patch system, please make "
                "sure that that cubed_sphere_outer_boundary_radius >= 4 * "
                "cubed_sphere_inner_boundary_radius.");
  }
}

extern "C" void CapyrX_MultiPatch_Run_Unit_Tests(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose && !CCTK_EQUALS(patch_system, "none")) {
    CCTK_VINFO("Running unit tests for patch \"%s\"", patch_system);
  }

  if (CCTK_EQUALS(patch_system, "none")) {
    return;

  } else if (CCTK_EQUALS(patch_system, "Cartesian")) {
    Cartesian::PatchParams par{
        .ncells_i = cartesian_ncells_i,
        .ncells_j = cartesian_ncells_j,
        .ncells_k = cartesian_ncells_k,

        .xmin = cartesian_xmin,
        .ymin = cartesian_xmin,
        .zmin = cartesian_xmin,

        .xmax = cartesian_xmax,
        .ymax = cartesian_xmax,
        .zmax = cartesian_xmax,
    };

    const auto pass{Cartesian::unit_test(test_repetitions, test_seed, par)};

    if (pass) {
      CCTK_VINFO("Unit tests for patch \"%s\" have \033[1;32mPASSED\033[0m",
                 patch_system);
    }

  } else if (CCTK_EQUALS(patch_system, "Cubed sphere")) {
    CubedSphere::PatchParams par{
        .angular_cells = cubed_sphere_angular_cells + 2 * patch_overlap,
        .radial_cells = cubed_sphere_radial_cells + 2 * patch_overlap,
        .inner_boundary = cubed_sphere_inner_boundary_radius,
        .outer_boundary = cubed_sphere_outer_boundary_radius,
        .patch_overlap = patch_overlap};

    const auto pass{CubedSphere::unit_test(test_repetitions, test_seed, par)};

    if (pass) {
      CCTK_VINFO("Unit tests for patch \"%s\" have \033[1;32mPASSED\033[0m",
                 patch_system);
    }

  } else {
    CCTK_VERROR("Cannot run unit tests for unknown patch system \"%s\"",
                patch_system);
  }
}

} // namespace CapyrX::MultiPatch
