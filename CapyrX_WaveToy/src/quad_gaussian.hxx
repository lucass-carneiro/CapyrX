#ifndef CAPYRX_WAVETOY_QUAD_GAUSSIAN_HXX
#define CAPYRX_WAVETOY_QUAD_GAUSSIAN_HXX

#include <cctk.h>
#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy::quad_gauss {

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phi(T A, T sigma, T R0, T x0, T y0,
                                             T z0, T x, T y, T z) noexcept
    -> T {
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
static inline auto CCTK_HOST CCTK_DEVICE Pi(T A, T sigma, T R0, T x0, T y0,
                                            T z0, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dx(T A, T sigma, T R0, T x0, T y0,
                                            T z0, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (x - x0 - y + y0)) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dy(T A, T sigma, T R0, T x0, T y0,
                                            T z0, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (x - x0 - y + y0)) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dz(T A, T sigma, T R0, T x0, T y0,
                                            T z0, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return -0.5 *
             (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
              z0) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (z - z0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (z - z0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidx(T A, T sigma, T R0, T x0, T y0,
                                               T z0, T x, T y, T z) noexcept
    -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidy(T A, T sigma, T R0, T x0, T y0,
                                               T z0, T x, T y, T z) noexcept
    -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidz(T A, T sigma, T R0, T x0, T y0,
                                               T z0, T x, T y, T z) noexcept
    -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDxdx(T A, T sigma, T R0, T x0, T y0,
                                               T z0, T x, T y, T z) noexcept
    -> T {
  using std::sqrt, std::pow, std::exp;
  return (-2 * A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) -
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (2 * A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * A * sqrt(15 / M_PI) * pow(x - x0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 3) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDydy(T A, T sigma, T R0, T x0, T y0,
                                               T z0, T x, T y, T z) noexcept
    -> T {
  using std::sqrt, std::pow, std::exp;
  return (-2 * A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) +
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (2 * A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * sqrt(15 / M_PI) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * A * sqrt(15 / M_PI) * pow(y - y0, 2) * (x - x0 + y - y0) *
          (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 3) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDzdz(T A, T sigma, T R0, T x0, T y0,
                                               T z0, T x, T y, T z) noexcept
    -> T {
  using std::sqrt, std::pow, std::exp;
  return (2 * A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (z - z0) * z0) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (z - z0) * z0) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (2. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) -
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (A * sqrt(15 / M_PI) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          pow(z - z0, 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (4. *
              exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2));
}

} // namespace CapyrX::WaveToy::quad_gauss

#endif // CAPYRX_WAVETOY_QUAD_GAUSSIAN_HXX
