#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include "standing_wave.hxx"
#include "gaussian.hxx"
#include "quad_gaussian.hxx"

namespace CapyrX::WaveToy_SecondOrder {

extern "C" void CapyrX_WaveToy_SecondOrder_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_SecondOrder_Initial;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          u(p.I) = sw::u(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
          rho(p.I) = sw::rho(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z);
        });

  } else if (CCTK_EQUALS(initial_condition, "gaussian")) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          u(p.I) = gauss::u(amplitude, gaussian_width, t, x, y, z);
          rho(p.I) = gauss::rho(amplitude, gaussian_width, t, x, y, z);
        });

  } else if (CCTK_EQUALS(initial_condition, "Quadrupolar Gaussian")) {
    constexpr auto eps{std::numeric_limits<double>::epsilon()};
    const auto sigma{gaussian_width};
    const auto R0{quad_gaussian_R0};
    const auto x0{quad_gaussian_x0};
    const auto y0{quad_gaussian_y0};
    const auto z0{quad_gaussian_z0};

    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::sqrt;
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto rxy2{(x - x0) * (x - x0) + (y - y0) * (y - y0)};
          const auto r{sqrt(rxy2 + (z - z0) * (z - z0))};

          if (r < eps || rxy2 < eps) {
            u(p.I) = 0.0;
            rho(p.I) = 0.0;
          } else {
            u(p.I) = quad_gauss::u(amplitude, sigma, R0, x0, y0, z0, x, y, z);
            rho(p.I) =
                quad_gauss::rho(amplitude, sigma, R0, x0, y0, z0, x, y, z);
          }
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace CapyrX::WaveToy_SecondOrder
