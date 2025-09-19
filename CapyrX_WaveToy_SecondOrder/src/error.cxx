#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include "standing_wave.hxx"
#include "gaussian.hxx"

namespace CapyrX::WaveToy_SecondOrder {

extern "C" void CapyrX_WaveToy_SecondOrder_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_SecondOrder_Error;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto u_exact{
              sw::u(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};
          const auto rho_exact{
              sw::rho(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};

          u_error(p.I) = u(p.I) - u_exact;
          rho_error(p.I) = rho(p.I) - rho_exact;
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto u_exact{gauss::u(amplitude, gaussian_width, t, x, y, z)};
          const auto rho_exact{
              gauss::rho(amplitude, gaussian_width, t, x, y, z)};

          u_error(p.I) = u(p.I) - u_exact;
          rho_error(p.I) = rho(p.I) - rho_exact;
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace CapyrX::WaveToy_SecondOrder
