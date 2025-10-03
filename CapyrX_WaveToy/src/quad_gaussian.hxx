#ifndef CAPYRX_WAVETOY_QUAD_GAUSSIAN_HXX
#define CAPYRX_WAVETOY_QUAD_GAUSSIAN_HXX

#include <cctk.h>
#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy::quad_gauss {

template <typename T>
static inline auto phi(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                       T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
         (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)),
                  2) /
              (2. * pow(sigma, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2)) *
          (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto Pi(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto Dx(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (2 * A * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto Dy(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (2 * A * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto Dz(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (-2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) * z0) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) * (z - z0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (z - z0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5));
}

template <typename T>
static inline auto dPidx(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto dPidy(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto dPidz(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto dDxdx(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (-8 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (4 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (4 * A * (x - x0) * (x - x0 + y - y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (4 * A * (x - x0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (8 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) -
         (4 * A * (x - x0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (4 * A * (x - x0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (8 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (2 * A * (x - x0) * (x - x0 + y - y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * A * (x - x0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (4 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * A * (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (4 * A * (x - x0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (4 * A * (x - x0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * A * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 3) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto dDydy(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (-8 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (4 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (4 * A * (y - y0) * (x - x0 + y - y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (4 * A * (y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (8 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) +
         (4 * A * (y - y0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (4 * A * (y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (8 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (2 * A * (y - y0) * (x - x0 + y - y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * A * (y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (4 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (2 * A * (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (4 * A * (y - y0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (4 * A * (y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * A * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 3) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto dDzdz(T A, T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (8 * A * (x - x0 + y - y0) * (x - x0 - y + y0) * (z - z0) * z0) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (4 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (z - z0) * z0) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (2 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (8 * A * (x - x0 + y - y0) * (x - x0 - y + y0) * pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) -
         (A * (x - x0 + y - y0) * (x - x0 - y + y0) * pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2));
}

} // namespace CapyrX::WaveToy::quad_gauss

#endif // CAPYRX_WAVETOY_QUAD_GAUSSIAN_HXX
