#include "mpwavetoy.hxx"

namespace MultiPatchWaveToy {

extern "C" void MultiPatchWaveToy_Error(CCTK_ARGUMENTS) {
  using std::pow;

  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0{0}, rho0{0};
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, vcoordx(p.I), vcoordy(p.I),
                        vcoordz(p.I), u0, rho0);
          u_err(p.I) = u(p.I) - u0;
          rho_err(p.I) = rho(p.I) - rho0;
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0{0}, rho0{0};
          gaussian(amplitude, gaussian_width, cctk_time, vcoordx(p.I),
                   vcoordy(p.I), vcoordz(p.I), u0, rho0);
          u_err(p.I) = u(p.I) - u0;
          rho_err(p.I) = rho(p.I) - rho0;
        });

  } else if (CCTK_EQUALS(initial_condition, "plane wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0{0}, rho0{0};
          plane_wave(amplitude, plane_wave_frequency, plane_wave_nx,
                     plane_wave_ny, plane_wave_nz, cctk_time, vcoordx(p.I),
                     vcoordy(p.I), vcoordz(p.I), u0, rho0);
          u_err(p.I) = u(p.I) - u0;
          rho_err(p.I) = rho(p.I) - rho0;
        });

  } else {
    CCTK_ERROR("Unable to compute error: unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy