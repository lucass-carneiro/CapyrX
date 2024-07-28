#include "multipatch.hxx"
#include "CParameters.h"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>

namespace MultiPatch {

PatchTransformations::PatchTransformations()
    : // Cartesian
      cartesian_xmax([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_xmax;
      }()),
      cartesian_xmin([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_xmin;
      }()),
      cartesian_ymax([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ymax;
      }()),
      cartesian_ymin([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ymin;
      }()),
      cartesian_zmax([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_zmax;
      }()),
      cartesian_zmin([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_zmin;
      }()),
      cartesian_ncells_i([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ncells_i;
      }()),
      cartesian_ncells_j([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ncells_j;
      }()),
      cartesian_ncells_k([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ncells_k;
      }()),

      // Cubed sphere
      cubed_sphere_rmin([] {
        DECLARE_CCTK_PARAMETERS;
        return cubed_sphere_rmin;
      }()),
      cubed_sphere_rmax([] {
        DECLARE_CCTK_PARAMETERS;
        return cubed_sphere_rmax;
      }()),

      // Swirl
      swirl_ncells_i([] {
        DECLARE_CCTK_PARAMETERS;
        return swirl_ncells_i;
      }()),
      swirl_ncells_j([] {
        DECLARE_CCTK_PARAMETERS;
        return swirl_ncells_j;
      }()),
      swirl_ncells_k([] {
        DECLARE_CCTK_PARAMETERS;
        return swirl_ncells_k;
      }()),

      // Cake
      cake_outer_boundary_radius([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_outer_boundary_radius;
      }()),
      cake_inner_boundary_radius([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_inner_boundary_radius;
      }()),
      cake_cartesian_ncells_i([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_cartesian_ncells_i;
      }()),
      cake_cartesian_ncells_j([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_cartesian_ncells_j;
      }()),
      cake_cartesian_ncells_k([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_cartesian_ncells_k;
      }()),
      cake_angular_cells([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_angular_cells;
      }()),
      cake_radial_cells([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_radial_cells;
      }()),

      // Two Cubes
      two_cubes_xmin([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_xmin;
      }()),
      two_cubes_xmax([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_xmax;
      }()),
      two_cubes_ymin([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_ymin;
      }()),
      two_cubes_ymax([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_ymax;
      }()),
      two_cubes_zmin([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_zmin;
      }()),
      two_cubes_zmax([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_zmax;
      }()),
      two_cubes_delta_y([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_delta_y;
      }()),
      two_cubes_ncells_left([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_ncells_left;
      }()),
      two_cubes_ncells_right([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_ncells_right;
      }()),
      two_cubes_ncells_y([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_ncells_y;
      }()),
      two_cubes_ncells_z([] {
        DECLARE_CCTK_PARAMETERS;
        return two_cubes_ncells_z;
      }()),

      // Thornburg06 Patch
      thornburg06_outer_boundary_radius([] {
        DECLARE_CCTK_PARAMETERS;
        return thornburg06_outer_boundary_radius;
      }()),
      thornburg06_inner_boundary_radius([] {
        DECLARE_CCTK_PARAMETERS;
        return thornburg06_inner_boundary_radius;
      }()),
      thornburg06_angular_cells([] {
        DECLARE_CCTK_PARAMETERS;
        return thornburg06_angular_cells;
      }()),
      thornburg06_radial_cells([] {
        DECLARE_CCTK_PARAMETERS;
        return thornburg06_radial_cells;
      }()) {}

std::unique_ptr<PatchSystem> the_patch_system;

namespace {
template <typename T, int D>
CCTK_DEVICE CCTK_HOST inline T det(const vec<vec<T, D>, D> &A) {
  return A(0)(0) * (A(1)(1) * A(2)(2) - A(1)(2) * A(2)(1)) +
         A(0)(1) * (A(1)(2) * A(2)(0) - A(1)(0) * A(2)(2)) +
         A(0)(2) * (A(1)(0) * A(2)(1) - A(1)(1) * A(2)(0));
}
} // namespace

// Aliased functions

extern "C" CCTK_INT
MultiPatch1_GetSystemSpecification(CCTK_INT *restrict const npatches) {
  *npatches = the_patch_system->num_patches();
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetPatchSpecification(
    const CCTK_INT ipatch, const CCTK_INT size, CCTK_INT *restrict const ncells,
    CCTK_REAL *restrict const xmin, CCTK_REAL *restrict const xmax) {
  assert(ipatch >= 0 && ipatch < the_patch_system->num_patches());
  assert(size == dim);
  const Patch &patch = the_patch_system->patches.at(ipatch);
  for (int d = 0; d < dim; ++d) {
    ncells[d] = patch.ncells[d];
    xmin[d] = patch.xmin[d];
    xmax[d] = patch.xmax[d];
  }
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetBoundarySpecification2(
    const CCTK_INT ipatch, const CCTK_INT size,
    CCTK_INT *restrict const is_interpatch_boundary) {
  assert(ipatch >= 0 && ipatch < the_patch_system->num_patches());
  assert(size == 2 * dim);
  const Patch &patch = the_patch_system->patches.at(ipatch);
  for (int f = 0; f < 2; ++f)
    for (int d = 0; d < dim; ++d)
      is_interpatch_boundary[2 * d + f] = !patch.faces[f][d].is_outer_boundary;
  return 0;
}

extern "C" void MultiPatch1_GlobalToLocal2(
    const CCTK_INT npoints, const CCTK_REAL *restrict const globalsx,
    const CCTK_REAL *restrict const globalsy,
    const CCTK_REAL *restrict const globalsz, CCTK_INT *restrict const patches,
    CCTK_REAL *restrict const localsx, CCTK_REAL *restrict const localsy,
    CCTK_REAL *restrict const localsz) {
  const PatchTransformations &transformations =
      the_patch_system->transformations;
  const auto &global2local = *transformations.global2local;

  for (int n = 0; n < npoints; ++n) {
    const vec<CCTK_REAL, dim> x{globalsx[n], globalsy[n], globalsz[n]};
    const auto patch_a = global2local(transformations, x);
    const auto patch = std::get<0>(patch_a);
    const auto a = std::get<1>(patch_a);
    patches[n] = patch;
    localsx[n] = a(0);
    localsy[n] = a(1);
    localsz[n] = a(2);
  }
}

extern "C" void MultiPatch1_LocalToGlobal2(
    const CCTK_INT npoints, const CCTK_INT *restrict const patches,
    const CCTK_REAL *restrict const localsx,
    const CCTK_REAL *restrict const localsy,
    const CCTK_REAL *restrict const localsz, CCTK_REAL *restrict const globalsx,
    CCTK_REAL *restrict const globalsy, CCTK_REAL *restrict const globalsz) {

  const auto &pt{the_patch_system->transformations};
  const auto &local2global{*pt.local2global};

  for (CCTK_INT i = 0; i < npoints; i++) {
    const auto patch{patches[i]};
    const vec<CCTK_REAL, dim> local_vars{localsx[i], localsy[i], localsz[i]};

    const auto global_vars{local2global(pt, patch, local_vars)};
    globalsx[i] = global_vars(0);
    globalsy[i] = global_vars(1);
    globalsz[i] = global_vars(2);
  }
}

// Scheduled functions

extern "C" int MultiPatch_Setup() {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(patch_system, "none"))
    the_patch_system = nullptr;
  else if (CCTK_EQUALS(patch_system, "Cartesian"))
    the_patch_system = std::make_unique<PatchSystem>(SetupCartesian());
  else if (CCTK_EQUALS(patch_system, "Cubed sphere"))
    the_patch_system = std::make_unique<PatchSystem>(SetupCubedSphere());
  else if (CCTK_EQUALS(patch_system, "Swirl"))
    the_patch_system = std::make_unique<PatchSystem>(SetupSwirl());
  else if (CCTK_EQUALS(patch_system, "Cake"))
    the_patch_system = std::make_unique<PatchSystem>(SetupCake());
  else if (CCTK_EQUALS(patch_system, "Two Cubes"))
    the_patch_system = std::make_unique<PatchSystem>(SetupTwoCubes());
  else if (CCTK_EQUALS(patch_system, "Thornburg06"))
    the_patch_system = std::make_unique<PatchSystem>(SetupTwoCubes());
  else
    CCTK_VERROR("Unknown patch system \"%s\"", patch_system);

  return 0;
}

extern "C" void MultiPatch_Coordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatch_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  const PatchTransformations pt = the_patch_system->transformations;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] ARITH_DEVICE(const Loop::PointDesc &p) ARITH_INLINE {
        const vec<CCTK_REAL, dim> a = {p.x, p.y, p.z};

        const auto d2J_tuple = pt.d2local_dglobal2(pt, p.patch, a);
        const auto &x = std::get<0>(d2J_tuple);
        const auto &J = std::get<1>(d2J_tuple);
        const auto &dJ = std::get<2>(d2J_tuple);

        vcoordx(p.I) = x(0);
        vcoordy(p.I) = x(1);
        vcoordz(p.I) = x(2);

        vJ_da_dx(p.I) = J(0)(0);
        vJ_da_dy(p.I) = J(0)(1);
        vJ_da_dz(p.I) = J(0)(2);
        vJ_db_dx(p.I) = J(1)(0);
        vJ_db_dy(p.I) = J(1)(1);
        vJ_db_dz(p.I) = J(1)(2);
        vJ_dc_dx(p.I) = J(2)(0);
        vJ_dc_dy(p.I) = J(2)(1);
        vJ_dc_dz(p.I) = J(2)(2);

        vdJ_d2a_dxdx(p.I) = dJ(0)(0, 0);
        vdJ_d2a_dxdy(p.I) = dJ(0)(0, 1);
        vdJ_d2a_dxdz(p.I) = dJ(0)(0, 2);
        vdJ_d2a_dydy(p.I) = dJ(0)(1, 1);
        vdJ_d2a_dydz(p.I) = dJ(0)(1, 2);
        vdJ_d2a_dzdz(p.I) = dJ(0)(2, 2);
        vdJ_d2b_dxdx(p.I) = dJ(1)(0, 0);
        vdJ_d2b_dxdy(p.I) = dJ(1)(0, 1);
        vdJ_d2b_dxdz(p.I) = dJ(1)(0, 2);
        vdJ_d2b_dydy(p.I) = dJ(1)(1, 1);
        vdJ_d2b_dydz(p.I) = dJ(1)(1, 2);
        vdJ_d2b_dzdz(p.I) = dJ(1)(2, 2);
        vdJ_d2c_dxdx(p.I) = dJ(2)(0, 0);
        vdJ_d2c_dxdy(p.I) = dJ(2)(0, 1);
        vdJ_d2c_dxdz(p.I) = dJ(2)(0, 2);
        vdJ_d2c_dydy(p.I) = dJ(2)(1, 1);
        vdJ_d2c_dydz(p.I) = dJ(2)(1, 2);
        vdJ_d2c_dzdz(p.I) = dJ(2)(2, 2);
      });

  // TODO: Currently, on Debug runs, this code fails to set the ccoord[xyz] and
  // cvol variables correctly by calling the multipatch coordinate
  // transformations
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] ARITH_DEVICE(const Loop::PointDesc &p) ARITH_INLINE {
        const vec<CCTK_REAL, dim> a = {p.x, p.y, p.z};

        const auto d2J_tuple = pt.d2local_dglobal2(pt, p.patch, a);
        const auto &x = std::get<0>(d2J_tuple);
        const auto &J = std::get<1>(d2J_tuple);
        // const auto &dJ = std::get<2>(d2J_tuple);
        const auto detJ = det(J);

        // TODO: Is it detJ or sqrt(detJ)? As far as we
        // know no sqrt is necessary. Note that detJ is
        // negative in the cake patch, so that causes NaNs
        const auto vol = (p.dx * p.dy * p.dz) * detJ;

        ccoordx(p.I) = 0.0; // x(0);
        ccoordy(p.I) = 0.0; // x(1);
        ccoordz(p.I) = 0.0; // x(2);
        cvol(p.I) = 0.0;
      });
}

/**
 * TODO: Fill with more parameter checks, if appropriate
 */
extern "C" void MultiPatch_Check_Parameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatch_Check_Parameters;
  DECLARE_CCTK_PARAMETERS;

  if (!(cake_outer_boundary_radius > 4 * cake_inner_boundary_radius))
    CCTK_PARAMWARN(
        "Make sure that outer boundary radius is larger than 4 times "
        "the inner boudary radius.");
}

} // namespace MultiPatch
