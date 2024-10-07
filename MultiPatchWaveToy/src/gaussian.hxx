#ifndef MULTIPATCH_WAVE_TOY_GAUSSIAN_HXX
#define MULTIPATCH_WAVE_TOY_GAUSSIAN_HXX

namespace MultiPatchWaveToy::gaussian {

template <typename T>
static inline auto phi(T A, T gaussian_width, T t, T x, T y,
                       T z) noexcept -> T {
  return 0.0; // TODO
}

template <typename T>
static inline auto Pi(T A, T gaussian_width, T t, T x, T y, T z) noexcept -> T {
  return 0.0; // TODO
}

template <typename T>
static inline auto Dx(T A, T gaussian_width, T t, T x, T y, T z) noexcept -> T {
  return 0.0; // TODO
}

template <typename T>
static inline auto Dy(T A, T gaussian_width, T t, T x, T y, T z) noexcept -> T {
  return 0.0; // TODO
}

template <typename T>
static inline auto Dz(T A, T gaussian_width, T t, T x, T y, T z) noexcept -> T {
  return 0.0; // TODO
}

} // namespace MultiPatchWaveToy::gaussian

#endif // MULTIPATCH_WAVE_TOY_GAUSSIAN_HXX