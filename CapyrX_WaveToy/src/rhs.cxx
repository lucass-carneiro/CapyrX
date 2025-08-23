#include "cctk_Config.h"
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

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE c4o_0_1_0(
    const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[1] + p.I) - 8 * gf(-p.DI[1] + p.I) +
                 8 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * p.DX[1])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE c4o_0_0_1(
    const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[2] + p.I) - 8 * gf(-p.DI[2] + p.I) +
                 8 * gf(p.DI[2] + p.I) - gf(2 * p.DI[2] + p.I)};
  const auto den{1.0 / (12 * p.DX[2])};
  return num * den;
}

template <std::size_t dir>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE diss_5(
    const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const auto fac{(1.0 / 64.0) * (1.0 / p.DX[dir])};
  const auto stencil{gf(p.I - 3 * p.DI[dir]) - 6.0 * gf(p.I - 2 * p.DI[dir]) +
                     15.0 * gf(p.I - 1 * p.DI[dir]) - 20.0 * gf(p.I) +
                     15.0 * gf(p.I + 1 * p.DI[dir]) -
                     6.0 * gf(p.I + 2 * p.DI[dir]) + gf(p.I + 3 * p.DI[dir])};
  return fac * stencil;
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

        const auto diss_phi{
            dissipation_epsilon *
            (diss_5<0>(p, phi) + diss_5<1>(p, phi) + diss_5<2>(p, phi))};

        const auto diss_Pi{
            dissipation_epsilon *
            (diss_5<0>(p, Pi) + diss_5<1>(p, Pi) + diss_5<2>(p, Pi))};

        const auto diss_Dx{
            dissipation_epsilon *
            (diss_5<0>(p, Dx) + diss_5<1>(p, Dx) + diss_5<2>(p, Dx))};

        const auto diss_Dy{
            dissipation_epsilon *
            (diss_5<0>(p, Dy) + diss_5<1>(p, Dy) + diss_5<2>(p, Dy))};

        const auto diss_Dz{
            dissipation_epsilon *
            (diss_5<0>(p, Dz) + diss_5<1>(p, Dz) + diss_5<2>(p, Dz))};

        phi_rhs(p.I) = Pi(p.I) + diss_phi;
        Pi_rhs(p.I) = gDx.dx + gDy.dy + gDz.dz + diss_Pi;
        Dx_rhs(p.I) = gPi.dx + diss_Dx;
        Dy_rhs(p.I) = gPi.dy + diss_Dy;
        Dz_rhs(p.I) = gPi.dz + diss_Dz;
      });
}

extern "C" void CapyrX_WaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace CapyrX::WaveToy
