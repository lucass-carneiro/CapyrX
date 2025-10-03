#ifndef CAPYRX_WAVETOY_SECOND_ORDER_GAUSSIAN_HXX
#define CAPYRX_WAVETOY_SECOND_ORDER_GAUSSIAN_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy_SecondOrder::gauss {

template <typename T>
static inline auto CCTK_DEVICE u(T A, T W, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return (2 * A * t) / (exp(pow(t, 2) / (2. * pow(W, 2))) * pow(W, 2));
  } else {
    return (A *
            (exp(-0.5 * pow(t - sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)), 2) /
                 pow(W, 2)) -
             exp(-0.5 * pow(t + sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)), 2) /
                 pow(W, 2)))) /
           sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  }
}

template <typename T>
static inline auto CCTK_DEVICE rho(T A, T W, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return (2 * A * (-t + W) * (t + W)) /
           (exp(pow(t, 2) / (2. * pow(W, 2))) * pow(W, 4));
  } else {
    return (2 * A *
            (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 cosh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)) -
             t * sinh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) /
                (2. * pow(W, 2))) *
            pow(W, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
  }
}

} // namespace CapyrX::WaveToy_SecondOrder::gauss

#endif // CAPYRX_WAVETOY_SECOND_ORDER_GAUSSIAN_HXX