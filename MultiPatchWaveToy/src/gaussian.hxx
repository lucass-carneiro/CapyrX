#ifndef MULTIPATCH_WAVE_TOY_GAUSSIAN_HXX
#define MULTIPATCH_WAVE_TOY_GAUSSIAN_HXX

#include <cmath>

namespace MultiPatchWaveToy::gauss {

template <typename T>
static inline auto phi(T W, T A, T t, T x, T y, T z) noexcept -> T {
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
static inline auto Pi(T W, T A, T t, T x, T y, T z) noexcept -> T {
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

template <typename T>
static inline auto Dx(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A * x *
            (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
             t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto Dy(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A * y *
            (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
             t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto Dz(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A * z *
            (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
             t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
             exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                 pow(W, 2)) *
                 (pow(W, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto dPidx(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (-2 * A * x *
            (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 (pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 cosh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)) -
             t * (pow(W, 2) + 2 * (pow(x, 2) + pow(y, 2) + pow(z, 2))) *
                 sinh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) /
                (2. * pow(W, 2))) *
            pow(W, 4) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto dPidy(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (-2 * A * y *
            (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 (pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 cosh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)) -
             t * (pow(W, 2) + 2 * (pow(x, 2) + pow(y, 2) + pow(z, 2))) *
                 sinh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) /
                (2. * pow(W, 2))) *
            pow(W, 4) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto dPidz(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (-2 * A * z *
            (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 (pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 cosh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)) -
             t * (pow(W, 2) + 2 * (pow(x, 2) + pow(y, 2) + pow(z, 2))) *
                 sinh((t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                      pow(W, 2)))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) /
                (2. * pow(W, 2))) *
            pow(W, 4) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  }
}

template <typename T>
static inline auto dDxdx(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return (A *
            ((-1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                       pow(W, 2))) *
                 pow(W, 4) * (2 * pow(x, 2) - pow(y, 2) - pow(z, 2)) -
             2 *
                 (1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                          pow(W, 2))) *
                 t * pow(x, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5) +
             (-1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                       pow(W, 2))) *
                 pow(x, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) *
                 (pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) +
             pow(W, 2) *
                 ((-1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                            pow(W, 2))) *
                      pow(x, 4) -
                  2 *
                      (1 +
                       exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                           pow(W, 2))) *
                      t * pow(x, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  (pow(y, 2) + pow(z, 2)) *
                      ((-1 +
                        exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                            pow(W, 2))) *
                           pow(y, 2) -
                       pow(z, 2) - t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
                       exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                           pow(W, 2)) *
                           (pow(z, 2) -
                            t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))))))) /
           (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
                 2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                (2. * pow(W, 2))) *
            pow(W, 4) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.5));
  }
}

template <typename T>
static inline auto dDydy(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return -(
        (A *
         ((-1 +
           exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / pow(W, 2))) *
              pow(W, 4) * (pow(x, 2) - 2 * pow(y, 2) + pow(z, 2)) -
          pow(y, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) *
              ((-1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                         pow(W, 2))) *
                   pow(t, 2) -
               2 *
                   (1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                            pow(W, 2))) *
                   t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
               (-1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                         pow(W, 2))) *
                   (pow(x, 2) + pow(y, 2) + pow(z, 2))) +
          pow(W, 2) *
              ((-1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                         pow(W, 2))) *
                   pow(x, 4) +
               pow(y, 4) - pow(z, 4) +
               t * (2 * pow(y, 2) - pow(z, 2)) *
                   sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
               pow(x, 2) *
                   (2 *
                        (-1 +
                         exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                             pow(W, 2))) *
                        pow(z, 2) -
                    (1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                             pow(W, 2))) *
                        t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) +
               exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                   pow(W, 2)) *
                   (-pow(y, 4) + pow(z, 4) +
                    2 * t * pow(y, 2) *
                        sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                    t * pow(z, 2) *
                        sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)))))) /
        (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
              2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
             (2. * pow(W, 2))) *
         pow(W, 4) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.5)));
  }
}

template <typename T>
static inline auto dDzdz(T W, T A, T t, T x, T y, T z) noexcept -> T {
  using std::sqrt, std::cosh, std::sinh;

  const auto r{sqrt(x * x + y * y + z * z)};

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    return 0;
  } else {
    return -(
        (A *
         ((-1 +
           exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / pow(W, 2))) *
              pow(W, 4) * (pow(x, 2) + pow(y, 2) - 2 * pow(z, 2)) +
          2 *
              (1 + exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                       pow(W, 2))) *
              t * pow(z, 2) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5) -
          (-1 +
           exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / pow(W, 2))) *
              pow(z, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) *
              (pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2)) +
          pow(W, 2) *
              (-pow(pow(x, 2) + pow(y, 2), 2) + pow(z, 4) -
               t * pow(x, 2) *
                   sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                   pow(W, 2)) *
                   t *
                   pow(x, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               t * pow(y, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                   pow(W, 2)) *
                   t *
                   pow(y, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
               2 * t *
                   pow(z, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
               2 * exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / pow(W, 2)) *
                   t * pow(z, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
               exp((2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
                   pow(W, 2)) *
                   (pow(pow(x, 2) + pow(y, 2), 2) - pow(z, 4))))) /
        (exp((pow(t, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2) +
              2 * t * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) /
             (2. * pow(W, 2))) *
         pow(W, 4) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 2.5)));
  }
}

} // namespace MultiPatchWaveToy::gauss

#endif // MULTIPATCH_WAVE_TOY_GAUSSIAN_HXX
