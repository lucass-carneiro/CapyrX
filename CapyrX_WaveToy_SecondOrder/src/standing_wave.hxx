#ifndef CAPYRX_WAVETOY_SECOND_ORDER_STANDING_WAVE_HXX
#define CAPYRX_WAVETOY_SECOND_ORDER_STANDING_WAVE_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy_SecondOrder::sw {

template <typename T>
static inline auto CCTK_DEVICE u(T A, T kx, T ky, T kz, T t, T x, T y,
                                 T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return A * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_DEVICE rho(T A, T kx, T ky, T kz, T t, T x, T y,
                                   T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * omega * pi * cos(2 * kx * pi * x) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * omega * pi * t);
}

} // namespace CapyrX::WaveToy_SecondOrder::sw

#endif // CAPYRX_WAVETOY_SECOND_ORDER_STANDING_WAVE_HXX