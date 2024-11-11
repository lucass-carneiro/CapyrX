#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>

namespace MultiPatchWaveToy {

extern "C" void MultiPatchWaveToy_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto dPi{MultiPatch::GlobalDerivatives::c4o_1(cctkGH, p, Pi)};
        const auto dDx{MultiPatch::GlobalDerivatives::c4o_1(cctkGH, p, Dx)};
        const auto dDy{MultiPatch::GlobalDerivatives::c4o_1(cctkGH, p, Dy)};
        const auto dDz{MultiPatch::GlobalDerivatives::c4o_1(cctkGH, p, Dz)};

        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = dDx.dx + dDy.dy + dDz.dz;
        Dx_rhs(p.I) = dPi.dx;
        Dy_rhs(p.I) = dPi.dy;
        Dz_rhs(p.I) = dPi.dz;
      });
}

extern "C" void MultiPatchWaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace MultiPatchWaveToy
