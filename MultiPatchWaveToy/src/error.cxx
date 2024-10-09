// clang-format off
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
// clang-format on

#include "standing_wave.hxx"

namespace MultiPatchWaveToy {

using namespace Arith;

extern "C" void MultiPatchWaveToy_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_all<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
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

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    CCTK_ERROR("Unimplemented initial condition");

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy