#ifndef CAPYRX_PATCH_CARTESIAN_HXX
#define CAPYRX_PATCH_CARTESIAN_HXX

#include "multipatch.hxx"

#include <loop_device.hxx>
#include <tuple.hxx>

namespace CapyrX::MultiPatch::Cartesian {

struct PatchParams {
  CCTK_INT ncells_i{10};
  CCTK_INT ncells_j{10};
  CCTK_INT ncells_k{10};

  CCTK_REAL xmin{1};
  CCTK_REAL ymin{1};
  CCTK_REAL zmin{1};

  CCTK_REAL xmax{11};
  CCTK_REAL ymax{11};
  CCTK_REAL zmax{11};
};

CCTK_HOST CCTK_DEVICE auto
global2local(const PatchParams &p,
             const svec_t &global_coords) -> std_tuple<int, svec_t>;

CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &p, int patch,
                                        const svec_t &local_coords) -> svec_t;

CCTK_HOST CCTK_DEVICE auto
dlocal_dglobal(const PatchParams &p, int patch,
               const svec_t &local_coords) -> std_tuple<svec_t, jac_t>;

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2_fun(const PatchParams &p, int patch,
                                                const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t>;

auto make_system(const PatchParams &p) -> PatchSystem;

} // namespace CapyrX::MultiPatch::Cartesian

#endif // CAPYRX_PATCH_CARTESIAN_HXX