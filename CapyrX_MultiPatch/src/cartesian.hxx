#ifndef CAPYRX_PATCH_CARTESIAN_HXX
#define CAPYRX_PATCH_CARTESIAN_HXX

#include "type_aliases.hxx"

#include <loop.hxx>

namespace CapyrX::MultiPatch {

struct CartesianPatchParams {
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

class CartesianPatchSystem {
private:
  const CartesianPatchParams params{};

public:
  CartesianPatchSystem(CartesianPatchParams &&p) : params(p) {}
  CartesianPatchSystem(const CartesianPatchParams &p) : params(p) {}

  CCTK_HOST CCTK_DEVICE auto
  global2local(const svec_t &global_vars) -> std_tuple<int, svec_t>;

  CCTK_HOST CCTK_DEVICE auto local2global(int patch,
                                          const svec_t &local_vars) -> svec_t;

  CCTK_HOST CCTK_DEVICE auto dlocal_dglobal(int patch, const svec_t &local_vars)
      -> std_tuple<svec_t, jac_t>(*);

  CCTK_HOST CCTK_DEVICE auto d2local_dglobal2_fun(
      int patch, const svec_t &local_vars) -> std_tuple<svec_t, jac_t, djac_t>;
};

} // namespace CapyrX::MultiPatch

#endif // CAPYRX_PATCH_CARTESIAN_HXX