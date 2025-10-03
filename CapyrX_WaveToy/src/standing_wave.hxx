#ifndef CAPYRX_WAVETOY_STANDING_WAVE_HXX
#define CAPYRX_WAVETOY_STANDING_WAVE_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy::sw {

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phi(T A, T kx, T ky, T kz, T t, T x,
                                             T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return A * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * omega * pi * cos(2 * kx * pi * x) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * omega * pi * t);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dx(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * kx * pi * cos(2 * omega * pi * t) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * kx * pi * x);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dy(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * ky * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * kz * pi * z) * sin(2 * ky * pi * y);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dz(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * kz * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidx(T A, T kx, T ky, T kz, T t, T x,
                                               T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 4 * A * kx * omega * pow(pi, 2) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * omega * pi * t) * sin(2 * kx * pi * x);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidy(T A, T kx, T ky, T kz, T t, T x,
                                               T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 4 * A * ky * omega * pow(pi, 2) * cos(2 * kx * pi * x) *
         cos(2 * kz * pi * z) * sin(2 * omega * pi * t) * sin(2 * ky * pi * y);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidz(T A, T kx, T ky, T kz, T t, T x,
                                               T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 4 * A * kz * omega * pow(pi, 2) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * sin(2 * omega * pi * t) * sin(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDxdx(T A, T kx, T ky, T kz, T t, T x,
                                               T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -4 * A * pow(kx, 2) * pow(pi, 2) * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDydy(T A, T kx, T ky, T kz, T t, T x,
                                               T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -4 * A * pow(ky, 2) * pow(pi, 2) * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDzdz(T A, T kx, T ky, T kz, T t, T x,
                                               T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -4 * A * pow(kz, 2) * pow(pi, 2) * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

} // namespace CapyrX::WaveToy::sw

#endif // CAPYRX_WAVETOY_STANDING_WAVE_HXX