#ifndef CAPYRX_PATCH_CUBED_SPHERE_HXX
#define CAPYRX_PATCH_CUBED_SPHERE_HXX

#include "multipatch.hxx"

#include <loop_device.hxx>
#include <tuple.hxx>

namespace CapyrX::MultiPatch::CubedSphere {

/**
 * Parameters that define a CubedSphere patch
 */
struct PatchParams {
  CCTK_INT angular_cells{10};
  CCTK_INT radial_cells{10};

  CCTK_REAL inner_boundary{1};
  CCTK_REAL outer_boundary{4};

  CCTK_INT patch_overlap{0};
};

CCTK_HOST CCTK_DEVICE auto global2local(const PatchParams &par,
                                        const svec_t &global_coords)
    -> std_tuple<int, svec_t>;

CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &par, int patch,
                                        const svec_t &local_coords) -> svec_t;

CCTK_HOST CCTK_DEVICE auto dlocal_dglobal(const PatchParams &par, int patch,
                                          const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t>;

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &par, int patch,
                                            const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t>;

auto make_system(const PatchParams &par) -> PatchSystem;

auto unit_test(std::size_t repetitions, std::size_t seed,
               const PatchParams &par) -> bool;

} // namespace CapyrX::MultiPatch::CubedSphere

#endif // CAPYRX_PATCH_CUBED_SPHERE_HXX
