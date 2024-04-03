#include "mpwavetoy.hxx"

#include <cmath>
#include <limits>

namespace MultiPatchWaveToy {

void standing_wave(const CCTK_REAL A, const CCTK_REAL kx, const CCTK_REAL ky,
                   const CCTK_REAL kz, const CCTK_REAL t, const CCTK_REAL x,
                   const CCTK_REAL y, const CCTK_REAL z, CCTK_REAL &u,
                   CCTK_REAL &rho) noexcept {
  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const CCTK_REAL pi = acos(-CCTK_REAL(1));
  const CCTK_REAL omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
      cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  rho = A * (-2 * pi * omega) * sin(2 * pi * omega * t) * cos(2 * pi * kx * x) *
        cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

void gaussian(const CCTK_REAL A, const CCTK_REAL W, const CCTK_REAL t,
              const CCTK_REAL x, const CCTK_REAL y, const CCTK_REAL z,
              CCTK_REAL &u, CCTK_REAL &rho) noexcept {
  using std::exp, std::pow, std::sqrt;

  const CCTK_REAL r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  const auto f = [&](const CCTK_REAL v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<CCTK_REAL>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
    rho = -2 / pow(W, 4) * f(t) * (pow(t, 2) - pow(W, 2));
  } else {
    u = (f(t - r) - f(t + r)) / r;
    rho = -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow(W, 2) * r);
  }
}

extern "C" void MultiPatchWaveToy_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, vcoordx(p.I), vcoordy(p.I),
                        vcoordz(p.I), u(p.I), rho(p.I));
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gaussian(amplitude, gaussian_width, cctk_time, vcoordx(p.I),
                   vcoordy(p.I), vcoordz(p.I), u(p.I), rho(p.I));
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace MultiPatchWaveToy