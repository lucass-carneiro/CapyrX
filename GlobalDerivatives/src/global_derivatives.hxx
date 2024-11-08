#ifndef CAPYRX_GLOBAL_DERIVATIVES_HXX
#define CAPYRX_GLOBAL_DERIVATIVES_HXX

#include <cctk.h>
#include <loop_device.hxx>

namespace MultiPatch::GlobalDerivatives {

struct FirstDerivativeResults {
  CCTK_REAL dx{0.0};
  CCTK_REAL dy{0.0};
  CCTK_REAL dz{0.0};
};

struct SecondDerivativeResults {
  CCTK_REAL dxdx{0.0};
  CCTK_REAL dxdy{0.0};
  CCTK_REAL dxdz{0.0};
  CCTK_REAL dydy{0.0};
  CCTK_REAL dydz{0.0};
  CCTK_REAL dzdz{0.0};
};

auto CCTK_HOST CCTK_DEVICE c4o_1(
    const CCTK_POINTER_TO_CONST cctkGH_, const Loop::PointDesc &p,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> FirstDerivativeResults;

auto CCTK_HOST CCTK_DEVICE c4o_2(
    const CCTK_POINTER_TO_CONST cctkGH_, const Loop::PointDesc &p,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> SecondDerivativeResults;

} // namespace MultiPatch::GlobalDerivatives

#endif // CAPYRX_GLOBAL_DERIVATIVES_HXX