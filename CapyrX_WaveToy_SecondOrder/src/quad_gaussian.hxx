#ifndef CAPYRX_WAVETOY_SECOND_ORDER_QUAD_GAUSSIAN_HXX
#define CAPYRX_WAVETOY_SECOND_ORDER_QUAD_GAUSSIAN_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy_SecondOrder::quad_gauss {

template <typename T>
static inline auto CCTK_DEVICE u(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                                 T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
         (4. *
          exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)),
                  2) /
              (2. * pow(sigma, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2)) *
          (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_DEVICE rho(T A, T sigma, T R0, T x0, T y0, T z0, T x,
                                   T y, T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

} // namespace CapyrX::WaveToy_SecondOrder::quad_gauss

#endif // CAPYRX_WAVETOY_SECOND_ORDER_QUAD_GAUSSIAN_HXX