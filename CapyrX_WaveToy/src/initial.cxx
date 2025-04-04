#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <limits>
#include <loop_device.hxx>

#include "standing_wave.hxx"
#include "gaussian.hxx"
#include "quad_gaussian.hxx"

namespace MultiPatchWaveToy {

using namespace Arith;

extern "C" void MultiPatchWaveToy_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          phi(p.I) = sw::phi(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
          Pi(p.I) = sw::Pi(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
          Dx(p.I) = sw::Dx(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
          Dy(p.I) = sw::Dy(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
          Dz(p.I) = sw::Dz(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          phi(p.I) = gauss::phi(amplitude, gaussian_width, t, x, y, z);
          Pi(p.I) = gauss::Pi(amplitude, gaussian_width, t, x, y, z);
          Dx(p.I) = gauss::Dx(amplitude, gaussian_width, t, x, y, z);
          Dy(p.I) = gauss::Dy(amplitude, gaussian_width, t, x, y, z);
          Dz(p.I) = gauss::Dz(amplitude, gaussian_width, t, x, y, z);
        });
  } else if (CCTK_EQUALS(initial_condition, "Quadrupolar Gaussian")) {
    constexpr auto eps{std::numeric_limits<double>::epsilon()};
    const auto sigma{gaussian_width};
    const auto R0{quad_gaussian_R0};
    const auto x0{quad_gaussian_x0};
    const auto y0{quad_gaussian_y0};
    const auto z0{quad_gaussian_z0};

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::sqrt;
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto rxy2{(x - x0) * (x - x0) + (y - y0) * (y - y0)};
          const auto r{sqrt(rxy2 + (z - z0) * (z - z0))};

          if (r < eps || rxy2 < eps) {
            phi(p.I) = 0.0;
            Pi(p.I) = 0.0;
            Dx(p.I) = 0.0;
            Dy(p.I) = 0.0;
            Dz(p.I) = 0.0;
          } else {
            phi(p.I) = quad_gauss::phi(sigma, R0, x0, y0, z0, x, y, z);
            Pi(p.I) = quad_gauss::Pi(sigma, R0, x0, y0, z0, x, y, z);
            Dx(p.I) = quad_gauss::Dx(sigma, R0, x0, y0, z0, x, y, z);
            Dy(p.I) = quad_gauss::Dy(sigma, R0, x0, y0, z0, x, y, z);
            Dz(p.I) = quad_gauss::Dz(sigma, R0, x0, y0, z0, x, y, z);
          }
        });
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy
