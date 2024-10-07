// clang-format off
#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
// clang-format on

#include "standing_wave.hxx"

namespace MultiPatchWaveToy {

extern "C" void MultiPatchWaveToy_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
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
    CCTK_ERROR("Unimplemented initial condition");

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy