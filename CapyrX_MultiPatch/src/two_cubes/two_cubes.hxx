#ifndef CAPYRX_PATCH_TWO_CUBES_HXX
#define CAPYRX_PATCH_TWO_CUBES_HXX

#include "multipatch.hxx"

#include <loop_device.hxx>
#include <tuple.hxx>

namespace CapyrX::MultiPatch::TwoCubes {

/**
 * Parameters that define a TwoCubes patch
 */
struct PatchParams {
  CCTK_INT ncells_x{10};
  CCTK_INT ncells_y{10};
  CCTK_INT ncells_z{10};

  CCTK_REAL xmin{-1};
  CCTK_REAL ymin{-1};
  CCTK_REAL zmin{-1};

  CCTK_REAL xmax{1};
  CCTK_REAL ymax{1};
  CCTK_REAL zmax{1};

  CCTK_REAL bend{0};

  CCTK_INT patch_overlap{0};
};

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
global2local(const PatchParams &par, const svec_t &global_coords)
    -> std_tuple<int, svec_t>;

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
local2global(const PatchParams &par, int patch, const svec_t &local_coords)
    -> svec_t;

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
dlocal_dglobal(const PatchParams &par, int patch, const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t>;

CCTK_HOST CCTK_DEVICE CAPYRX_EXTERNAL auto
d2local_dglobal2(const PatchParams &par, int patch, const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t>;

auto make_system(const PatchParams &par) -> PatchSystem;

auto unit_test(std::size_t repetitions, std::size_t seed,
               const PatchParams &par) -> bool;

} // namespace CapyrX::MultiPatch::TwoCubes

#endif // CAPYRX_PATCH_TWO_CUBES_HXX
