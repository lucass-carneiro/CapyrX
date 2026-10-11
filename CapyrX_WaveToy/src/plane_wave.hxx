#ifndef CAPYRX_WAVETOY_PLANE_WAVE_HXX
#define CAPYRX_WAVETOY_PLANE_WAVE_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

// Travelling plane wave matching LlamaWaveToy's "plane" initial data
// (src/init.F90): phi = A cos(2*pi*(k.(x-x0) + omega*(t-t0))),
// omega = |k|, so that it is an exact Minkowski solution of the first-order
// wave system evolved here (Pi = d phi/dt, Di = d phi/dxi). The +omega*t sign
// is kept identical to LlamaWaveToy so both codes share one analytic solution.

namespace CapyrX::WaveToy::pw {

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phase(T kx, T ky, T kz, T x0, T y0,
                                               T z0, T t0, T t, T x, T y,
                                               T z) noexcept -> T {
  using std::sqrt, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 2 * pi *
         (kx * (x - x0) + ky * (y - y0) + kz * (z - z0) + omega * (t - t0));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phi(T A, T kx, T ky, T kz, T x0, T y0,
                                             T z0, T t0, T t, T x, T y,
                                             T z) noexcept -> T {
  using std::cos;
  return A * cos(phase(kx, ky, kz, x0, y0, z0, t0, t, x, y, z));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T A, T kx, T ky, T kz, T x0, T y0,
                                            T z0, T t0, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sqrt, std::sin, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * pi * omega *
         sin(phase(kx, ky, kz, x0, y0, z0, t0, t, x, y, z));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dx(T A, T kx, T ky, T kz, T x0, T y0,
                                            T z0, T t0, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sin, std::acos;
  const auto pi{acos(T{-1})};
  return -2 * A * pi * kx *
         sin(phase(kx, ky, kz, x0, y0, z0, t0, t, x, y, z));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dy(T A, T kx, T ky, T kz, T x0, T y0,
                                            T z0, T t0, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sin, std::acos;
  const auto pi{acos(T{-1})};
  return -2 * A * pi * ky *
         sin(phase(kx, ky, kz, x0, y0, z0, t0, t, x, y, z));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dz(T A, T kx, T ky, T kz, T x0, T y0,
                                            T z0, T t0, T t, T x, T y,
                                            T z) noexcept -> T {
  using std::sin, std::acos;
  const auto pi{acos(T{-1})};
  return -2 * A * pi * kz *
         sin(phase(kx, ky, kz, x0, y0, z0, t0, t, x, y, z));
}

} // namespace CapyrX::WaveToy::pw

#endif // CAPYRX_WAVETOY_PLANE_WAVE_HXX
