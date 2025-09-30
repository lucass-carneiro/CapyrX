#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>
#include <newradx.hxx>

#include "local_derivatives.hxx"

namespace CapyrX::WaveToy_SecondOrder {

extern "C" void CapyrX_WaveToy_SecondOrder_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_SecondOrder_RHS;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;
  using namespace CapyrX::MultiPatch::GlobalDerivatives;

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

        const auto d2ux{project_second(lu, l2u, jac, djac)};

        u_rhs(p.I) = rho(p.I);
        rho_rhs(p.I) = d2ux.dxdx + d2ux.dydy + d2ux.dzdz;
      });

  if (use_newradx) {
    NewRadX::NewRadX_Apply(cctkGH, u, u_rhs, vcoordx, vcoordy, vcoordz, 0.0,
                           1.0, rad_power);
    NewRadX::NewRadX_Apply(cctkGH, rho, rho_rhs, vcoordx, vcoordy, vcoordz, 0.0,
                           1.0, rad_power + 1.0);
  }

  if (add_dissipation) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto diss_u{diss_5<0>(p, u) + diss_5<1>(p, u) +
                            diss_5<2>(p, u)};

          const auto diss_rho{diss_5<0>(p, rho) + diss_5<1>(p, rho) +
                              diss_5<2>(p, rho)};

          u_rhs(p.I) += dissipation_epsilon * diss_u;
          rho_rhs(p.I) += dissipation_epsilon * diss_rho;
        });
  }
}

extern "C" void CapyrX_WaveToy_SecondOrder_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace CapyrX::WaveToy_SecondOrder
