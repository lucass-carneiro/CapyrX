#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include "local_derivatives.hxx"
#include <global_derivatives.hxx>

namespace CapyrX::WaveToy_SecondOrder {

extern "C" void CapyrX_WaveToy_SecondOrder_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_SecondOrder_Energy;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;
  using namespace MultiPatch::GlobalDerivatives;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const LocalFirstDerivatives lu{
            .da = d_x<0>(p, u), .db = d_x<1>(p, u), .dc = d_x<2>(p, u)};

        const LocalSecondDerivatives l2u{
            .dada = d_xx<0>(p, u),
            .dadb = d_xy<0, 1>(p, u),
            .dadc = d_xy<0, 2>(p, u),
            .dbdb = d_xx<1>(p, u),
            .dbdc = d_xy<1, 2>(p, u),
            .dcdc = d_xx<2>(p, u),
        };

        const Jacobians jac{VERTEX_JACOBIANS(p)};
        const JacobianDerivatives djac{VERTEX_DJACOBIANS(p)};

        const auto d2x{project_second(lu, l2u, jac, djac)};

        const auto du2{d2x.dxdx * d2x.dxdx + d2x.dydy * d2x.dydy +
                       d2x.dzdz * d2x.dzdz};
        const auto rho2 = rho(p.I) * rho(p.I);

        energy(p.I) = (rho2 + du2) / 2.0;
      });
}

} // namespace CapyrX::WaveToy_SecondOrder
