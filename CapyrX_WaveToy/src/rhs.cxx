#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>

namespace CapyrX::WaveToy {

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
    c4o_1_0_0(const Loop::PointDesc &p,
              const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[0] + p.I) - 8 * gf(-p.DI[0] + p.I) +
                 8 * gf(p.DI[0] + p.I) - gf(2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * p.DX[0])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_1_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[1] + p.I) - 8 * gf(-p.DI[1] + p.I) +
                 8 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * p.DX[1])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_0_1(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[2] + p.I) - 8 * gf(-p.DI[2] + p.I) +
                 8 * gf(p.DI[2] + p.I) - gf(2 * p.DI[2] + p.I)};
  const auto den{1.0 / (12 * p.DX[2])};
  return num * den;
}

extern "C" void CapyrX_WaveToy_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;
  using namespace CapyrX::MultiPatch::GlobalDerivatives;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_HOST CCTK_DEVICE(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const LocalFirstDerivatives lPi{c4o_1_0_0(p, Pi), c4o_0_1_0(p, Pi),
                                        c4o_0_0_1(p, Pi)};

        const LocalFirstDerivatives lDx{c4o_1_0_0(p, Dx), c4o_0_1_0(p, Dx),
                                        c4o_0_0_1(p, Dx)};

        const LocalFirstDerivatives lDy{c4o_1_0_0(p, Dy), c4o_0_1_0(p, Dy),
                                        c4o_0_0_1(p, Dy)};

        const LocalFirstDerivatives lDz{c4o_1_0_0(p, Dz), c4o_0_1_0(p, Dz),
                                        c4o_0_0_1(p, Dz)};

        const Jacobians jac{VERTEX_JACOBIANS(p)};

        const auto gPi{project_first(lPi, jac)};
        const auto gDx{project_first(lDx, jac)};
        const auto gDy{project_first(lDy, jac)};
        const auto gDz{project_first(lDz, jac)};

        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = gDx.dx + gDy.dy + gDz.dz;
        Dx_rhs(p.I) = gPi.dx;
        Dy_rhs(p.I) = gPi.dy;
        Dz_rhs(p.I) = gPi.dz;
      });
}

extern "C" void CapyrX_WaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace CapyrX::WaveToy
