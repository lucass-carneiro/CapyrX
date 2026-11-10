#ifndef CAPYRX_WAVETOY_MULTIPOLE_GAUSSIAN_HXX
#define CAPYRX_WAVETOY_MULTIPOLE_GAUSSIAN_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy::mp_gauss {

// Pi of LlamaWaveToy's "multipole l=2, m=2" initial data at t = 0: a Gaussian
// shell in the time derivative with a real Y_22 angular profile, while
// phi = Dx = Dy = Dz = 0. The formula and its regularizations replicate
// LlamaWaveToy/src/init.F90 so that both codes start from identical data.
template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T A, T sigma, T r0, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::exp;

  const T eps{1.0e-10};
  const T pi{static_cast<T>(3.14159265358979323846)};

  const auto r1{sqrt(x * x + y * y + z * z + eps * eps)};

  auto rxy{sqrt(x * x + y * y)};
  if (rxy < 1.0e-12) {
    rxy = 1.0e-12;
  }

  const auto sin_theta{rxy / r1};
  const auto cos_phi{x / rxy};
  const auto sin_phi{y / rxy};

  const auto y22{0.25 * sqrt(15.0 / (2.0 * pi)) * sin_theta * sin_theta *
                 (cos_phi * cos_phi - sin_phi * sin_phi)};

  const auto r_eff{(r1 - r0) / sigma};

  return A * exp(-r_eff * r_eff) * y22;
}

} // namespace CapyrX::WaveToy::mp_gauss

#endif // CAPYRX_WAVETOY_MULTIPOLE_GAUSSIAN_HXX
