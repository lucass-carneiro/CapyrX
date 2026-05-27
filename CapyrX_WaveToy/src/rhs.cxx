#include "cctk_Config.h"
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>
#include <newradx.hxx>
#include <mat.hxx>

namespace CapyrX::WaveToy {

template <typename PointDescCallable>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
c4o_1_0_0(int dir, const Loop::PointDesc &p, PointDescCallable gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(dir, -2 * p.DI[0] + p.I) - 8 * gf(dir, -p.DI[0] + p.I) +
                 8 * gf(dir, p.DI[0] + p.I) - gf(dir, 2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * p.DX[0])};
  return num * den;
}

template <typename PointDescCallable>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
c4o_0_1_0(int dir, const Loop::PointDesc &p, PointDescCallable gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(dir, -2 * p.DI[1] + p.I) - 8 * gf(dir, -p.DI[1] + p.I) +
                 8 * gf(dir, p.DI[1] + p.I) - gf(dir, 2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * p.DX[1])};
  return num * den;
}

template <typename PointDescCallable>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
c4o_0_0_1(int dir, const Loop::PointDesc &p, PointDescCallable gf) noexcept
    -> CCTK_REAL {
  const auto num{gf(dir, -2 * p.DI[2] + p.I) - 8 * gf(dir, -p.DI[2] + p.I) +
                 8 * gf(dir, p.DI[2] + p.I) - gf(dir, 2 * p.DI[2] + p.I)};
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
  using namespace Arith;
  using namespace CapyrX::MultiPatch::GlobalDerivatives;
  using std::sqrt, std::pow;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /*
         * We assemble the metric matrix object using Arith::smat and use
         * it to compute determinants and inverses
         */
        const auto get_metric =
            [&](const vect<int, dim> &pI) CCTK_ATTRIBUTE_ALWAYS_INLINE {
              const smat<CCTK_REAL, 3> g_dd{gxx(pI), gxy(pI), gxz(pI),
                                            gyy(pI), gyz(pI), gzz(pI)};
              return g_dd;
            };

        /*
         * We use the Arith::vec object to represent rank-1
         * pysical tensors (AKA vectors)
         */
        const auto get_shift =
            [&](const vect<int, dim> &pI) CCTK_ATTRIBUTE_ALWAYS_INLINE {
              const vec<CCTK_REAL, 3> beta{betax(pI), betay(pI), betaz(pI)};
              return beta;
            };

        // Flux vector
        const auto flux =
            [&](int dir,
                const vect<int, dim> &pI) CCTK_ATTRIBUTE_ALWAYS_INLINE {
              const auto g_dd = get_metric(pI);
              const auto detg = calc_det(g_dd);
              const auto sqrtdetg = sqrt(detg);
              const auto g_uu = calc_inv(g_dd, detg);

              const vec<CCTK_REAL, 3> beta = get_shift(pI);
              const vec<CCTK_REAL, 3> D{Dx(pI), Dy(pI), Dz(pI)};

              // Now we compute the flux
              return sqrtdetg *
                     (beta(dir) * Pi(pI) +
                      alp(pI) * (g_uu(dir, 0) * D(0) + g_uu(dir, 1) * D(1) +
                                 g_uu(dir, 2) * D(2)) -
                      (1.0 / alp(pI)) * beta(dir) *
                          (beta(0) * D(0) + beta(1) * D(1) + beta(2) * D(2)));
            };

        /*
         * alpha * \Pi product
         * Because of how we defined our derivative functions, we will have to
         * pass a dummy direction argument to this scalar
         */
        const auto alpha_Pi =
            [&](int, const vect<int, dim> &pI)
                CCTK_ATTRIBUTE_ALWAYS_INLINE { return alp(pI) * Pi(pI); };

        // Local derivatives of the flux vector
        const LocalFirstDerivatives l_dFx = {
            c4o_1_0_0(0, p, flux),
            c4o_0_1_0(0, p, flux),
            c4o_0_0_1(0, p, flux),
        };

        const LocalFirstDerivatives l_dFy = {
            c4o_1_0_0(1, p, flux),
            c4o_0_1_0(1, p, flux),
            c4o_0_0_1(1, p, flux),
        };

        const LocalFirstDerivatives l_dFz = {
            c4o_1_0_0(2, p, flux),
            c4o_0_1_0(2, p, flux),
            c4o_0_0_1(2, p, flux),
        };

        // local derivatives of alpha * \Pi
        const LocalFirstDerivatives l_dalphaPi = {
            c4o_1_0_0(0, p, alpha_Pi),
            c4o_0_1_0(0, p, alpha_Pi),
            c4o_0_0_1(0, p, alpha_Pi),
        };

        // Jacobians
        const Jacobians jac{VERTEX_JACOBIANS(p)};

        // Global derivatives of the flux vector
        const auto g_dFx = project_first(l_dFx, jac);
        const auto g_dFy = project_first(l_dFy, jac);
        const auto g_dFz = project_first(l_dFz, jac);

        // Global derivatives of alpha * Pi
        const auto g_dalphaPi = project_first(l_dalphaPi, jac);

        /*
         * We follow Eqs.(56)-(58) of https://arxiv.org/abs/gr-qc/0507004.
         * IMPORTANT: These equations are not valid if the metric is not static.
         */
        const auto g_dd = get_metric(p.I);
        const auto detg = calc_det(g_dd);
        const auto sqrtdetg = sqrt(detg);

        const auto beta = get_shift(p.I);

        phi_rhs(p.I) = alpha_Pi(0, p.I);

        Pi_rhs(p.I) = (1 / alp(p.I)) *
                          (beta(0) * g_dalphaPi.dx + beta(1) * g_dalphaPi.dy +
                           beta(2) * g_dalphaPi.dz) +
                      (1 / sqrtdetg) * (g_dFx.dx + g_dFy.dy + g_dFz.dz);

        Dx_rhs(p.I) = g_dalphaPi.dx;
        Dy_rhs(p.I) = g_dalphaPi.dy;
        Dz_rhs(p.I) = g_dalphaPi.dz;
      });
}

extern "C" void CapyrX_WaveToy_Dissipation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_Dissipation;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;
  using std::sqrt;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL jac_norms[3] = {
            sqrt(vJ_da_dx(p.I) * vJ_da_dx(p.I) + vJ_da_dy(p.I) * vJ_da_dy(p.I) +
                 vJ_da_dz(p.I) * vJ_da_dz(p.I)),
            sqrt(vJ_db_dx(p.I) * vJ_db_dx(p.I) + vJ_db_dy(p.I) * vJ_db_dy(p.I) +
                 vJ_db_dz(p.I) * vJ_db_dz(p.I)),
            sqrt(vJ_dc_dx(p.I) * vJ_dc_dx(p.I) + vJ_dc_dy(p.I) * vJ_dc_dy(p.I) +
                 vJ_dc_dz(p.I) * vJ_dc_dz(p.I))};

        const auto diss_phi{dissipation_epsilon *
                            (jac_norms[0] * diss_5<0>(p, phi) +
                             jac_norms[1] * diss_5<1>(p, phi) +
                             jac_norms[2] * diss_5<2>(p, phi))};

        const auto diss_Pi{dissipation_epsilon *
                           (jac_norms[0] * diss_5<0>(p, Pi) +
                            jac_norms[1] * diss_5<1>(p, Pi) +
                            jac_norms[2] * diss_5<2>(p, Pi))};

        const auto diss_Dx{dissipation_epsilon *
                           (jac_norms[0] * diss_5<0>(p, Dx) +
                            jac_norms[1] * diss_5<1>(p, Dx) +
                            jac_norms[2] * diss_5<2>(p, Dx))};

        const auto diss_Dy{dissipation_epsilon *
                           (jac_norms[0] * diss_5<0>(p, Dy) +
                            jac_norms[1] * diss_5<1>(p, Dy) +
                            jac_norms[2] * diss_5<2>(p, Dy))};

        const auto diss_Dz{dissipation_epsilon *
                           (jac_norms[0] * diss_5<0>(p, Dz) +
                            jac_norms[1] * diss_5<1>(p, Dz) +
                            jac_norms[2] * diss_5<2>(p, Dz))};

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
