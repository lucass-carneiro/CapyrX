#ifndef CAPYRX_PATCH_LLAMA_HXX
#define CAPYRX_PATCH_LLAMA_HXX

#include "multipatch.hxx"

#include <loop_device.hxx>
#include <tuple.hxx>

namespace CapyrX::MultiPatch::Llama {

/**
 * Parameters that define a Llama patch: a central Cartesian cube surrounded by
 * six Thornburg06 spherical wedges. The cube resolution (cube_ncells_*) is
 * independent of the angular grid; the cube is a true [-R,R]^3 Cartesian patch.
 */
struct PatchParams {
  CCTK_INT angular_cells{10};
  CCTK_INT radial_cells{10};

  CCTK_INT cube_ncells_i{10};
  CCTK_INT cube_ncells_j{10};
  CCTK_INT cube_ncells_k{10};

  CCTK_REAL inner_boundary{1};
  CCTK_REAL outer_boundary{4};

  CCTK_INT patch_overlap{0};
};

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
local2global(const PatchParams &par, int patch, const svec_t &local_coords)
    -> svec_t;

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
dlocal_dglobal(const PatchParams &par, int patch, const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t>;

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
d2local_dglobal2(const PatchParams &par, int patch, const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t>;

} // namespace CapyrX::MultiPatch::Llama

#endif //  CAPYRX_PATCH_LLAMA_HXX
