#include "cctk_Config.h"
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>
#include <newradx.hxx>

namespace CapyrX::WaveToy {

template <typename PointDescCallable>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
c4o_1_0_0(const Loop::PointDesc &p, PointDescCallable gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[0] + p.I) - 8 * gf(-p.DI[0] + p.I) +
                 8 * gf(p.DI[0] + p.I) - gf(2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * p.DX[0])};
  return num * den;
}

template <typename PointDescCallable>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
c4o_0_1_0(const Loop::PointDesc &p, PointDescCallable gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[1] + p.I) - 8 * gf(-p.DI[1] + p.I) +
                 8 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * p.DX[1])};
  return num * den;
}

template <typename PointDescCallable>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
c4o_0_0_1(const Loop::PointDesc &p, PointDescCallable gf) noexcept
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

  grid.loop_int_device<
      0, 0, 0>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // We follow Eqs.(56)-(58) of https://arxiv.org/abs/gr-qc/0507004

    // flux_a = \partial_i(\alpha \Pi)
    const auto flux_a = [&](const vect<int, dim> &pI) {
      return alp(pI) * Pi(pI);
    };

    // flux_b_i = \partial_i (g^{1/2} \beta^i \Pi + \alpha g^{1/2} H^{ij} d_j)
    // where H^{ij} = g^{ij} - \alpha^{-2} \beta^i \beta^j
    const auto flux_b_x = [&](const vect<int, dim>
                                  &pI) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      return (sqrt(-(pow(gxz(pI), 2) * gyy(pI)) +
                   2 * gxy(pI) * gxz(pI) * gyz(pI) - gxx(pI) * pow(gyz(pI), 2) -
                   pow(gxy(pI), 2) * gzz(pI) + gxx(pI) * gyy(pI) * gzz(pI)) *
              (pow(alp(pI), 2) *
                   (Dz(pI) * gxz(pI) * gyy(pI) - Dz(pI) * gxy(pI) * gyz(pI) -
                    Dy(pI) * gxz(pI) * gyz(pI) + Dx(pI) * pow(gyz(pI), 2) +
                    Dy(pI) * gxy(pI) * gzz(pI) - Dx(pI) * gyy(pI) * gzz(pI)) -
               betax(pI) *
                   (betax(pI) * Dx(pI) + betay(pI) * Dy(pI) +
                    betaz(pI) * Dz(pI)) *
                   (pow(gxz(pI), 2) * gyy(pI) -
                    2 * gxy(pI) * gxz(pI) * gyz(pI) +
                    gxx(pI) * pow(gyz(pI), 2) + pow(gxy(pI), 2) * gzz(pI) -
                    gxx(pI) * gyy(pI) * gzz(pI)) +
               alp(pI) * betax(pI) *
                   (pow(gxz(pI), 2) * gyy(pI) -
                    2 * gxy(pI) * gxz(pI) * gyz(pI) +
                    gxx(pI) * pow(gyz(pI), 2) + pow(gxy(pI), 2) * gzz(pI) -
                    gxx(pI) * gyy(pI) * gzz(pI)) *
                   Pi(pI))) /
             (alp(pI) *
              (pow(gxz(pI), 2) * gyy(pI) - 2 * gxy(pI) * gxz(pI) * gyz(pI) +
               pow(gxy(pI), 2) * gzz(pI) +
               gxx(pI) * (pow(gyz(pI), 2) - gyy(pI) * gzz(pI))));
    };

    const auto flux_b_y = [&](const vect<int, dim>
                                  &pI) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      return (sqrt(-(pow(gxz(pI), 2) * gyy(pI)) +
                   2 * gxy(pI) * gxz(pI) * gyz(pI) - gxx(pI) * pow(gyz(pI), 2) -
                   pow(gxy(pI), 2) * gzz(pI) + gxx(pI) * gyy(pI) * gzz(pI)) *
              (pow(alp(pI), 2) *
                   (-(Dz(pI) * gxy(pI) * gxz(pI)) + Dy(pI) * pow(gxz(pI), 2) +
                    Dz(pI) * gxx(pI) * gyz(pI) - Dx(pI) * gxz(pI) * gyz(pI) -
                    Dy(pI) * gxx(pI) * gzz(pI) + Dx(pI) * gxy(pI) * gzz(pI)) -
               betay(pI) *
                   (betax(pI) * Dx(pI) + betay(pI) * Dy(pI) +
                    betaz(pI) * Dz(pI)) *
                   (pow(gxz(pI), 2) * gyy(pI) -
                    2 * gxy(pI) * gxz(pI) * gyz(pI) +
                    gxx(pI) * pow(gyz(pI), 2) + pow(gxy(pI), 2) * gzz(pI) -
                    gxx(pI) * gyy(pI) * gzz(pI)) +
               alp(pI) * betay(pI) *
                   (pow(gxz(pI), 2) * gyy(pI) -
                    2 * gxy(pI) * gxz(pI) * gyz(pI) +
                    gxx(pI) * pow(gyz(pI), 2) + pow(gxy(pI), 2) * gzz(pI) -
                    gxx(pI) * gyy(pI) * gzz(pI)) *
                   Pi(pI))) /
             (alp(pI) *
              (pow(gxz(pI), 2) * gyy(pI) - 2 * gxy(pI) * gxz(pI) * gyz(pI) +
               pow(gxy(pI), 2) * gzz(pI) +
               gxx(pI) * (pow(gyz(pI), 2) - gyy(pI) * gzz(pI))));
    };

    const auto flux_b_z = [&](const vect<int, dim>
                                  &pI) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      return (sqrt(-(pow(gxz(pI), 2) * gyy(pI)) +
                   2 * gxy(pI) * gxz(pI) * gyz(pI) - gxx(pI) * pow(gyz(pI), 2) -
                   pow(gxy(pI), 2) * gzz(pI) + gxx(pI) * gyy(pI) * gzz(pI)) *
              (pow(alp(pI), 2) *
                   (-(Dy(pI) * gxy(pI) * gxz(pI)) + Dx(pI) * gxz(pI) * gyy(pI) +
                    Dz(pI) * (pow(gxy(pI), 2) - gxx(pI) * gyy(pI)) +
                    Dy(pI) * gxx(pI) * gyz(pI) - Dx(pI) * gxy(pI) * gyz(pI)) -
               betaz(pI) *
                   (betax(pI) * Dx(pI) + betay(pI) * Dy(pI) +
                    betaz(pI) * Dz(pI)) *
                   (pow(gxz(pI), 2) * gyy(pI) -
                    2 * gxy(pI) * gxz(pI) * gyz(pI) +
                    gxx(pI) * pow(gyz(pI), 2) + pow(gxy(pI), 2) * gzz(pI) -
                    gxx(pI) * gyy(pI) * gzz(pI)) +
               alp(pI) * betaz(pI) *
                   (pow(gxz(pI), 2) * gyy(pI) -
                    2 * gxy(pI) * gxz(pI) * gyz(pI) +
                    gxx(pI) * pow(gyz(pI), 2) + pow(gxy(pI), 2) * gzz(pI) -
                    gxx(pI) * gyy(pI) * gzz(pI)) *
                   Pi(pI))) /
             (alp(pI) *
              (pow(gxz(pI), 2) * gyy(pI) - 2 * gxy(pI) * gxz(pI) * gyz(pI) +
               pow(gxy(pI), 2) * gzz(pI) +
               gxx(pI) * (pow(gyz(pI), 2) - gyy(pI) * gzz(pI))));
    };

    // Local derivatives of flux_a
    const LocalFirstDerivatives l_dflux_a{
        c4o_1_0_0(p, flux_a),
        c4o_0_1_0(p, flux_a),
        c4o_0_0_1(p, flux_a),
    };

    // Local derivatives of flux_b
    const LocalFirstDerivatives l_dflux_b_x{
        c4o_1_0_0(p, flux_b_x),
        c4o_0_1_0(p, flux_b_x),
        c4o_0_0_1(p, flux_b_x),
    };

    const LocalFirstDerivatives l_dflux_b_y{
        c4o_1_0_0(p, flux_b_y),
        c4o_0_1_0(p, flux_b_y),
        c4o_0_0_1(p, flux_b_y),
    };

    const LocalFirstDerivatives l_dflux_b_z{
        c4o_1_0_0(p, flux_b_z),
        c4o_0_1_0(p, flux_b_z),
        c4o_0_0_1(p, flux_b_z),
    };

    // Jacobians
    const Jacobians jac{VERTEX_JACOBIANS(p)};

    // Global derivatives of flux_a
    const auto g_dflux_a{project_first(l_dflux_a, jac)};

    // Global derivatives of flux_b
    const auto g_dflux_b_x{project_first(l_dflux_b_x, jac)};
    const auto g_dflux_b_y{project_first(l_dflux_b_y, jac)};
    const auto g_dflux_b_z{project_first(l_dflux_b_z, jac)};

    // RHS
    phi_rhs(p.I) = flux_a(p.I);

    Pi_rhs(p.I) =
        (betax(p.I) * g_dflux_a.dx + betay(p.I) * g_dflux_a.dy +
         betaz(p.I) * g_dflux_a.dz +
         (alp(p.I) * (g_dflux_b_x.dx + g_dflux_b_y.dy + g_dflux_b_z.dz)) /
             sqrt(-(pow(gxz(p.I), 2) * gyy(p.I)) +
                  2 * gxy(p.I) * gxz(p.I) * gyz(p.I) -
                  gxx(p.I) * pow(gyz(p.I), 2) - pow(gxy(p.I), 2) * gzz(p.I) +
                  gxx(p.I) * gyy(p.I) * gzz(p.I))) /
        alp(p.I);

    Dx_rhs(p.I) = g_dflux_a.dx;
    Dy_rhs(p.I) = g_dflux_a.dy;
    Dz_rhs(p.I) = g_dflux_a.dz;
  });
}

extern "C" void CapyrX_WaveToy_Dissipation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_Dissipation;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
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

        phi_rhs(p.I) += diss_phi;
        Pi_rhs(p.I) += diss_Pi;
        Dx_rhs(p.I) += diss_Dx;
        Dy_rhs(p.I) += diss_Dy;
        Dz_rhs(p.I) += diss_Dz;
      });
}

extern "C" void CapyrX_WaveToy_ApplyNewRadX(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_ApplyNewRadX;
  DECLARE_CCTK_PARAMETERS;

  using namespace NewRadX;

  NewRadX_Apply(cctkGH, phi, phi_rhs, NEWRADX_MULTIPATCH_QUANTITIES, 0.0, 1.0,
                rad_power);
  NewRadX_Apply(cctkGH, Pi, Pi_rhs, NEWRADX_MULTIPATCH_QUANTITIES, 0.0, 1.0,
                rad_power);

  NewRadX_Apply(cctkGH, Dx, Dx_rhs, NEWRADX_MULTIPATCH_QUANTITIES, 0.0, 1.0,
                rad_power);

  NewRadX_Apply(cctkGH, Dy, Dy_rhs, NEWRADX_MULTIPATCH_QUANTITIES, 0.0, 1.0,
                rad_power);

  NewRadX_Apply(cctkGH, Dz, Dz_rhs, NEWRADX_MULTIPATCH_QUANTITIES, 0.0, 1.0,
                rad_power);
}

extern "C" void CapyrX_WaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace CapyrX::WaveToy
