// clang-format off
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
// clang-format on

#include "standing_wave.hxx"
#include "gaussian.hxx"

namespace MultiPatchWaveToy {

using namespace Arith;

extern "C" void MultiPatchWaveToy_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int<0, 0, 0>(
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
    grid.loop_int<0, 0, 0>(
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
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy
