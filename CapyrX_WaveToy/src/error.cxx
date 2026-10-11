#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include "standing_wave.hxx"
#include "plane_wave.hxx"
#include "gaussian.hxx"

namespace CapyrX::WaveToy {

extern "C" void CapyrX_WaveToy_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_WaveToy_Error;
  DECLARE_CCTK_PARAMETERS;

  using namespace Loop;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto expected_phi{
              sw::phi(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};
          const auto expected_Pi{
              sw::Pi(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};
          const auto expected_Dx{
              sw::Dx(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};
          const auto expected_Dy{
              sw::Dy(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};
          const auto expected_Dz{
              sw::Dz(amplitude, wave_kx, wave_ky, wave_kz, t, x, y, z)};

          const auto actual_phi{phi(p.I)};
          const auto actual_Pi{Pi(p.I)};
          const auto actual_Dx{Dx(p.I)};
          const auto actual_Dy{Dy(p.I)};
          const auto actual_Dz{Dz(p.I)};

          phi_error(p.I) = fabs(expected_phi - actual_phi);
          Pi_error(p.I) = fabs(expected_Pi - actual_Pi);
          Dx_error(p.I) = fabs(expected_Dx - actual_Dx);
          Dy_error(p.I) = fabs(expected_Dy - actual_Dy);
          Dz_error(p.I) = fabs(expected_Dz - actual_Dz);
        });
  } else if (CCTK_EQUALS(initial_condition, "plane wave")) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto expected_phi{pw::phi(amplitude, wave_kx, wave_ky, wave_kz,
                                          wave_x0, wave_y0, wave_z0, wave_t0, t,
                                          x, y, z)};
          const auto expected_Pi{pw::Pi(amplitude, wave_kx, wave_ky, wave_kz,
                                        wave_x0, wave_y0, wave_z0, wave_t0, t, x,
                                        y, z)};
          const auto expected_Dx{pw::Dx(amplitude, wave_kx, wave_ky, wave_kz,
                                        wave_x0, wave_y0, wave_z0, wave_t0, t, x,
                                        y, z)};
          const auto expected_Dy{pw::Dy(amplitude, wave_kx, wave_ky, wave_kz,
                                        wave_x0, wave_y0, wave_z0, wave_t0, t, x,
                                        y, z)};
          const auto expected_Dz{pw::Dz(amplitude, wave_kx, wave_ky, wave_kz,
                                        wave_x0, wave_y0, wave_z0, wave_t0, t, x,
                                        y, z)};

          phi_error(p.I) = fabs(expected_phi - phi(p.I));
          Pi_error(p.I) = fabs(expected_Pi - Pi(p.I));
          Dx_error(p.I) = fabs(expected_Dx - Dx(p.I));
          Dy_error(p.I) = fabs(expected_Dy - Dy(p.I));
          Dz_error(p.I) = fabs(expected_Dz - Dz(p.I));
        });
  }
}

} // namespace CapyrX::WaveToy
