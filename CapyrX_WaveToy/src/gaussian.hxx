#ifndef CAPYRX_WAVETOY_GAUSSIAN_HXX
#define CAPYRX_WAVETOY_GAUSSIAN_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy::gauss {

// Eq. (74) of https://arxiv.org/abs/gr-qc/0507004
template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T A, T sigma, T r0, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::exp;

  const auto r = sqrt(x * x + y * y + z * z);
  const auto r_eff = r - r0;
  const auto r_eff2 = r_eff * r_eff;
  const auto sigma2 = sigma * sigma;

  return A * exp(-r_eff2 / sigma2);
}

} // namespace CapyrX::WaveToy::gauss

#endif // CAPYRX_WAVETOY_GAUSSIAN_HXX
