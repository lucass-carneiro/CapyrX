#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <mat.hxx>

namespace CapyrX::WaveToy {

extern "C" void CapyrX_WaveToy_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_Energy;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;
  using namespace Arith;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const smat<CCTK_REAL, 3> g_dd{gxx(p.I), gxy(p.I), gxz(p.I),
                                      gyy(p.I), gyz(p.I), gzz(p.I)};
        const auto detg = calc_det(g_dd);
        const auto sqrtdetg = sqrt(detg);
        const auto g_uu = calc_inv(g_dd, detg);

        energy(p.I) = 0.5 * sqrtdetg *
                      (Pi(p.I) * Pi(p.I) + g_uu(0, 0) * Dx(p.I) * Dx(p.I) +
                       2.0 * g_uu(0, 1) * Dx(p.I) * Dy(p.I) +
                       2.0 * g_uu(0, 2) * Dx(p.I) * Dz(p.I) +
                       g_uu(1, 1) * Dy(p.I) * Dy(p.I) +
                       2.0 * g_uu(1, 2) * Dy(p.I) * Dz(p.I) +
                       g_uu(2, 2) * Dz(p.I) * Dz(p.I));
      });
}

} // namespace CapyrX::WaveToy
