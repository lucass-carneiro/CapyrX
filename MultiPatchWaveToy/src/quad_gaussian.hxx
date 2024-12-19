#ifndef MULTIPATCH_WAVE_TOY_QUADRUPOLAR_GAUSSIAN_HXX
#define MULTIPATCH_WAVE_TOY_QUADRUPOLAR_GAUSSIAN_HXX

#include <cctk.h>

#include <cmath>

namespace MultiPatchWaveToy::quad_gauss {

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phi(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                       T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return ((x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
         (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)),
                  2) /
              (2. * pow(sigma, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2)) *
          (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dx(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (2 * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         ((x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         ((x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         ((x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * (x - x0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dy(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (2 * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         ((y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         ((x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         ((x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * (y - y0) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dz(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                      T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return ((x - x0 + y - y0) * (x - x0 - y + y0) *
          ((-2 * z0) / (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) -
           (2 * (z - z0) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2) +
           ((R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
            (z - z0) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               (pow(sigma, 2) *
                pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)))) /
         (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)),
                  2) /
              (2. * pow(sigma, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidx(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidy(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dPidz(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return 0;
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDxdx(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (-8 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (4 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (4 * (x - x0) * (x - x0 + y - y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (4 * (x - x0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (8 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) -
         (4 * (x - x0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (4 * (x - x0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (8 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (2 * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (2 * (x - x0) * (x - x0 + y - y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * (x - x0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (4 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         ((x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (4 * (x - x0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (4 * (x - x0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * pow(x - x0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 3) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDydy(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return (-8 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (4 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (4 * (y - y0) * (x - x0 + y - y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (4 * (y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (2 * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (8 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0)) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3)) -
         (5 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) +
         (4 * (y - y0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (4 * (y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (2 * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (8 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) +
         (pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)), 2) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 4) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)) -
         (2 * (y - y0) * (x - x0 + y - y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         (2 * (y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
         ((x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (4 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
         (2 * (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (4 * (y - y0) * (x - x0 + y - y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (4 * (y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) -
         (2 * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 2) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         (8 * pow(y - y0, 2) * (x - x0 + y - y0) * (x - x0 - y + y0) *
          (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(pow(x - x0, 2) + pow(y - y0, 2), 3) *
              (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)));
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE dDzdz(T sigma, T R0, T x0, T y0, T z0, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::pow, std::exp;
  return ((x - x0 + y - y0) * (x - x0 - y + y0) *
          (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
          (z - z0) *
          ((-2 * z0) / (pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2)) -
           (2 * (z - z0) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2) +
           ((R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
            (z - z0) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               (pow(sigma, 2) *
                pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              pow(sigma, 2) * (pow(x - x0, 2) + pow(y - y0, 2)) *
              sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) +
         ((x - x0 + y - y0) * (x - x0 - y + y0) *
          ((8 * (z - z0) * z0) /
               pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2) -
           (2 * (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
            (z - z0) * z0) /
               (pow(sigma, 2) *
                pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) -
           (2 * (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2) +
           ((R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               (pow(sigma, 2) *
                pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 1.5)) +
           (8 * pow(z - z0, 2) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 3) -
           (3 * (R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2))) *
            pow(z - z0, 2) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               (pow(sigma, 2) *
                pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2.5)) -
           (pow(z - z0, 2) *
            (pow(x - x0, 2) + pow(y - y0, 2) - 2 * z * z0 + pow(z0, 2))) /
               (pow(sigma, 2) *
                pow(pow(x - x0, 2) + pow(y - y0, 2) + pow(z - z0, 2), 2)))) /
             (exp(pow(R0 - sqrt(pow(x - x0, 2) + pow(y - y0, 2) +
                                pow(z - z0, 2)),
                      2) /
                  (2. * pow(sigma, 2))) *
              (pow(x - x0, 2) + pow(y - y0, 2)));
}

} // namespace MultiPatchWaveToy::quad_gauss

#endif // MULTIPATCH_WAVE_TOY_QUADRUPOLAR_GAUSSIAN_HXX
