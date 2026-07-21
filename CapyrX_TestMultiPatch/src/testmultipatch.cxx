#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>

#include <array>
#include <cmath>
#include <limits>

namespace CapyrX::TestMultiPatch {

// Must live at namespace scope (not inside a function): nvcc forbids a
// function-local type from being used as the type of a variable captured by
// an extended __device__ lambda.
enum class SmoothKind { z_global, parabola, one_over_r };

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
              CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
              CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
         cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dx(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                 CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                 CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -2 * A * kx * pi * cos(2 * omega * pi * t) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * kx * pi * x);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dy(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                 CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                 CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -2 * A * ky * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * kz * pi * z) * sin(2 * ky * pi * y);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dz(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                 CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                 CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -2 * A * kz * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dx2(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                  CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                  CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -4 * A * kx * kx * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dy2(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                  CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                  CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -4 * A * ky * ky * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dz2(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                  CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                  CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -4 * A * kz * kz * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dxy(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                  CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                  CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return 4 * A * kx * ky * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kz * pi * z) * sin(2 * kx * pi * x) * sin(2 * ky * pi * y);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dxz(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                  CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                  CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return 4 * A * kx * kz * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * ky * pi * y) * sin(2 * kx * pi * x) * sin(2 * kz * pi * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
standing_wave_dyz(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                  CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                  CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return 4 * A * ky * kz * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * sin(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) -> CCTK_REAL {
  return x * x + y * y + z * z;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dx(CCTK_REAL x, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0 * x;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dy(CCTK_REAL, CCTK_REAL y, CCTK_REAL) -> CCTK_REAL {
  return 2.0 * y;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dz(CCTK_REAL, CCTK_REAL, CCTK_REAL z) -> CCTK_REAL {
  return 2.0 * z;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dx2(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dy2(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dz2(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dxy(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 0.0;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dxz(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 0.0;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
parabola_dyz(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 0.0;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
    c4o_1_0_0(const Loop::PointDesc &p,
              const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[0] + p.I) - 8 * gf(-p.DI[0] + p.I) +
                 8 * gf(p.DI[0] + p.I) - gf(2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * p.DX[0])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_1_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[1] + p.I) - 8 * gf(-p.DI[1] + p.I) +
                 8 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * p.DX[1])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_0_1(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[2] + p.I) - 8 * gf(-p.DI[2] + p.I) +
                 8 * gf(p.DI[2] + p.I) - gf(2 * p.DI[2] + p.I)};
  const auto den{1.0 / (12 * p.DX[2])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_2_0_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{-30 * gf(p.I) - gf(-2 * p.DI[0] + p.I) +
                 16 * gf(-p.DI[0] + p.I) + 16 * gf(p.DI[0] + p.I) -
                 gf(2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * (p.DX[0] * p.DX[0]))};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_1_1_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{
      gf(-2 * p.DI[0] - 2 * p.DI[1] + p.I) -
      8 * gf(-p.DI[0] - 2 * p.DI[1] + p.I) +
      8 * gf(p.DI[0] - 2 * p.DI[1] + p.I) -
      gf(2 * p.DI[0] - 2 * p.DI[1] + p.I) -
      8 * gf(-2 * p.DI[0] - p.DI[1] + p.I) + 64 * gf(-p.DI[0] - p.DI[1] + p.I) -
      64 * gf(p.DI[0] - p.DI[1] + p.I) + 8 * gf(2 * p.DI[0] - p.DI[1] + p.I) +
      8 * gf(-2 * p.DI[0] + p.DI[1] + p.I) - 64 * gf(-p.DI[0] + p.DI[1] + p.I) +
      64 * gf(p.DI[0] + p.DI[1] + p.I) - 8 * gf(2 * p.DI[0] + p.DI[1] + p.I) -
      gf(-2 * p.DI[0] + 2 * p.DI[1] + p.I) +
      8 * gf(-p.DI[0] + 2 * p.DI[1] + p.I) -
      8 * gf(p.DI[0] + 2 * p.DI[1] + p.I) +
      gf(2 * p.DI[0] + 2 * p.DI[1] + p.I)};
  const auto den{1.0 / (144 * p.DX[0] * p.DX[1])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_1_0_1(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{
      gf(-2 * p.DI[0] - 2 * p.DI[2] + p.I) -
      8 * gf(-p.DI[0] - 2 * p.DI[2] + p.I) +
      8 * gf(p.DI[0] - 2 * p.DI[2] + p.I) -
      gf(2 * p.DI[0] - 2 * p.DI[2] + p.I) -
      8 * gf(-2 * p.DI[0] - p.DI[2] + p.I) + 64 * gf(-p.DI[0] - p.DI[2] + p.I) -
      64 * gf(p.DI[0] - p.DI[2] + p.I) + 8 * gf(2 * p.DI[0] - p.DI[2] + p.I) +
      8 * gf(-2 * p.DI[0] + p.DI[2] + p.I) - 64 * gf(-p.DI[0] + p.DI[2] + p.I) +
      64 * gf(p.DI[0] + p.DI[2] + p.I) - 8 * gf(2 * p.DI[0] + p.DI[2] + p.I) -
      gf(-2 * p.DI[0] + 2 * p.DI[2] + p.I) +
      8 * gf(-p.DI[0] + 2 * p.DI[2] + p.I) -
      8 * gf(p.DI[0] + 2 * p.DI[2] + p.I) +
      gf(2 * p.DI[0] + 2 * p.DI[2] + p.I)};
  const auto den{1.0 / (144 * p.DX[0] * p.DX[2])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_2_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{-30 * gf(p.I) - gf(-2 * p.DI[1] + p.I) +
                 16 * gf(-p.DI[1] + p.I) + 16 * gf(p.DI[1] + p.I) -
                 gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * (p.DX[1] * p.DX[1]))};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_1_1(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{
      gf(-2 * p.DI[1] - 2 * p.DI[2] + p.I) -
      8 * gf(-p.DI[1] - 2 * p.DI[2] + p.I) +
      8 * gf(p.DI[1] - 2 * p.DI[2] + p.I) -
      gf(2 * p.DI[1] - 2 * p.DI[2] + p.I) -
      8 * gf(-2 * p.DI[1] - p.DI[2] + p.I) + 64 * gf(-p.DI[1] - p.DI[2] + p.I) -
      64 * gf(p.DI[1] - p.DI[2] + p.I) + 8 * gf(2 * p.DI[1] - p.DI[2] + p.I) +
      8 * gf(-2 * p.DI[1] + p.DI[2] + p.I) - 64 * gf(-p.DI[1] + p.DI[2] + p.I) +
      64 * gf(p.DI[1] + p.DI[2] + p.I) - 8 * gf(2 * p.DI[1] + p.DI[2] + p.I) -
      gf(-2 * p.DI[1] + 2 * p.DI[2] + p.I) +
      8 * gf(-p.DI[1] + 2 * p.DI[2] + p.I) -
      8 * gf(p.DI[1] + 2 * p.DI[2] + p.I) +
      gf(2 * p.DI[1] + 2 * p.DI[2] + p.I)};
  const auto den{1.0 / (144 * p.DX[1] * p.DX[2])};
  return num * den;
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_0_0_2(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{-30 * gf(p.I) - gf(-2 * p.DI[2] + p.I) +
                 16 * gf(-p.DI[2] + p.I) + 16 * gf(p.DI[2] + p.I) -
                 gf(2 * p.DI[2] + p.I)};
  const auto den{1.0 / (12 * (p.DX[2] * p.DX[2]))};
  return num * den;
}

extern "C" void CapyrX_TestMultiPatch_write_test_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestMultiPatch_write_test_data;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(test_data, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          u(p.I) = standing_wave(A, kx, ky, kz, t, x, y, z);
        });

  } else if (CCTK_Equals(test_data, "parabola")) {
    grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                  [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        const auto x{vcoordx(p.I)};
                                        const auto y{vcoordy(p.I)};
                                        const auto z{vcoordz(p.I)};

                                        u(p.I) = parabola(x, y, z);
                                      });
  } else {
    CCTK_VERROR("Unknown test data type %s", test_data);
  }
}

extern "C" void CapyrX_TestMultiPatch_compute_interp_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestMultiPatch_compute_interp_error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(test_data, "standing wave")) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto interpd_u{u(p.I)};
          const auto exact_u{standing_wave(A, kx, ky, kz, t, x, y, z)};

          interp(p.I) = fabs(interpd_u - exact_u);
        });

  } else if (CCTK_Equals(test_data, "parabola")) {
    grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                  [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        using std::fabs;

                                        const auto x{vcoordx(p.I)};
                                        const auto y{vcoordy(p.I)};
                                        const auto z{vcoordz(p.I)};

                                        const auto interpd_u{u(p.I)};
                                        const auto exact_u{parabola(x, y, z)};

                                        interp(p.I) = fabs(interpd_u - exact_u);
                                      });
  }
}

extern "C" void CapyrX_TestMultiPatch_compute_deriv_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestMultiPatch_compute_deriv_error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(test_data, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using namespace CapyrX::MultiPatch::GlobalDerivatives;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto true_dfdx{standing_wave_dx(A, kx, ky, kz, t, x, y, z)};
          const auto true_dfdy{standing_wave_dy(A, kx, ky, kz, t, x, y, z)};
          const auto true_dfdz{standing_wave_dz(A, kx, ky, kz, t, x, y, z)};
          const auto true_d2fdx2{standing_wave_dx2(A, kx, ky, kz, t, x, y, z)};
          const auto true_d2fdy2{standing_wave_dy2(A, kx, ky, kz, t, x, y, z)};
          const auto true_d2fdz2{standing_wave_dz2(A, kx, ky, kz, t, x, y, z)};
          const auto true_d2fdxy{standing_wave_dxy(A, kx, ky, kz, t, x, y, z)};
          const auto true_d2fdxz{standing_wave_dxz(A, kx, ky, kz, t, x, y, z)};
          const auto true_d2fdyz{standing_wave_dyz(A, kx, ky, kz, t, x, y, z)};

          const LocalFirstDerivatives ldu{c4o_1_0_0(p, u), c4o_0_1_0(p, u),
                                          c4o_0_0_1(p, u)};

          const LocalSecondDerivatives ld2u{c4o_2_0_0(p, u), c4o_1_1_0(p, u),
                                            c4o_1_0_1(p, u), c4o_0_2_0(p, u),
                                            c4o_0_1_1(p, u), c4o_0_0_2(p, u)};

          const Jacobians jac{VERTEX_JACOBIANS(p)};
          const JacobianDerivatives djac{VERTEX_DJACOBIANS(p)};

          const auto first_derivs{project_first(ldu, jac)};
          const auto second_derivs{project_second(ldu, ld2u, jac, djac)};

          dfdx(p.I) = fabs(true_dfdx - first_derivs.dx);
          dfdy(p.I) = fabs(true_dfdy - first_derivs.dy);
          dfdz(p.I) = fabs(true_dfdz - first_derivs.dz);
          d2fdx2(p.I) = fabs(true_d2fdx2 - second_derivs.dxdx);
          d2fdy2(p.I) = fabs(true_d2fdy2 - second_derivs.dydy);
          d2fdz2(p.I) = fabs(true_d2fdz2 - second_derivs.dzdz);
          d2fdxy(p.I) = fabs(true_d2fdxy - second_derivs.dxdy);
          d2fdxz(p.I) = fabs(true_d2fdxz - second_derivs.dxdz);
          d2fdyz(p.I) = fabs(true_d2fdyz - second_derivs.dydz);
        });

  } else if (CCTK_Equals(test_data, "parabola")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using namespace CapyrX::MultiPatch::GlobalDerivatives;

          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto true_dfdx{parabola_dx(x, y, z)};
          const auto true_dfdy{parabola_dy(x, y, z)};
          const auto true_dfdz{parabola_dz(x, y, z)};
          const auto true_d2fdx2{parabola_dx2(x, y, z)};
          const auto true_d2fdy2{parabola_dy2(x, y, z)};
          const auto true_d2fdz2{parabola_dz2(x, y, z)};
          const auto true_d2fdxy{parabola_dxy(x, y, z)};
          const auto true_d2fdxz{parabola_dxz(x, y, z)};
          const auto true_d2fdyz{parabola_dyz(x, y, z)};

          const LocalFirstDerivatives ldu{c4o_1_0_0(p, u), c4o_0_1_0(p, u),
                                          c4o_0_0_1(p, u)};

          const LocalSecondDerivatives ld2u{c4o_2_0_0(p, u), c4o_1_1_0(p, u),
                                            c4o_1_0_1(p, u), c4o_0_2_0(p, u),
                                            c4o_0_1_1(p, u), c4o_0_0_2(p, u)};

          const Jacobians jac{VERTEX_JACOBIANS(p)};
          const JacobianDerivatives djac{VERTEX_DJACOBIANS(p)};

          const auto first_derivs{project_first(ldu, jac)};
          const auto second_derivs{project_second(ldu, ld2u, jac, djac)};

          dfdx(p.I) = fabs(true_dfdx - first_derivs.dx);
          dfdy(p.I) = fabs(true_dfdy - first_derivs.dy);
          dfdz(p.I) = fabs(true_dfdz - first_derivs.dz);
          d2fdx2(p.I) = fabs(true_d2fdx2 - second_derivs.dxdx);
          d2fdy2(p.I) = fabs(true_d2fdy2 - second_derivs.dydy);
          d2fdz2(p.I) = fabs(true_d2fdz2 - second_derivs.dzdz);
          d2fdxy(p.I) = fabs(true_d2fdxy - second_derivs.dxdy);
          d2fdxz(p.I) = fabs(true_d2fdxz - second_derivs.dxdz);
          d2fdyz(p.I) = fabs(true_d2fdyz - second_derivs.dydz);
        });
  }
}

extern "C" void CapyrX_TestMultiPatch_write_color(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestMultiPatch_write_color;
  DECLARE_CCTK_PARAMETERS;

  // Coloring is a 2 digit number. The first digit indicates the region
  // (interior, boundary, ghost, overlap band). The second digit indicates
  // the patch index. For this to hold, we need to have patch systems with 9
  // (or less) patches
  assert(cctk_npatches < 10);

  constexpr auto interior_marker = 10;
  constexpr auto boundary_marker = 20;
  constexpr auto ghost_marker = 30;
  constexpr auto overlap_marker = 40;

  // CapyrX_MultiPatch::patch_overlap is not this thorn's own parameter; read
  // it cross-thorn the same way CapyrX_MultiPatch itself reads
  // CarpetX::interpolation_order (see CapyrX_MultiPatch_Check_Parameters in
  // multipatch.cxx).
  const auto patch_overlap_param_ptr =
      CCTK_ParameterGet("patch_overlap", "CapyrX_MultiPatch", nullptr);
  if (patch_overlap_param_ptr == nullptr)
    CCTK_ERROR("Unable to read parameter patch_overlap from CapyrX_MultiPatch");
  const auto patch_overlap =
      *static_cast<const CCTK_INT *>(patch_overlap_param_ptr);

  // Per axis, the global-index thresholds beyond which an interior point
  // falls in the overlap band: the patch_overlap-wide slice of this patch's
  // own interior that exists purely so a neighboring patch has enough source
  // data to interpolate its ghost zone. Faces that are the true
  // physical/outer boundary (not another patch) never grow an overlap band,
  // so their threshold is set to never trigger.
  Loop::vect<int, Loop::dim> lo_threshold, hi_threshold;
  {
    std::array<CCTK_INT, 2 * Loop::dim> is_interpatch_face{};
    if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification2"))
      MultiPatch_GetBoundarySpecification2(grid.patch, 2 * Loop::dim,
                                           is_interpatch_face.data());

    for (int d = 0; d < Loop::dim; ++d) {
      lo_threshold[d] = is_interpatch_face[2 * d + 0]
                            ? grid.nghostzones[d] + patch_overlap
                            : -1;
      hi_threshold[d] = is_interpatch_face[2 * d + 1]
                            ? grid.gsh[d] - grid.nghostzones[d] - patch_overlap
                            : std::numeric_limits<int>::max();
    }
  }

  const auto lbnd{grid.lbnd};

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool in_overlap_band = false;
        for (int d = 0; d < Loop::dim; ++d) {
          const auto gI{p.I[d] + lbnd[d]};
          if (gI < lo_threshold[d] || gI >= hi_threshold[d]) {
            in_overlap_band = true;
            break;
          }
        }
        color(p.I) =
            (in_overlap_band ? overlap_marker : interior_marker) + p.patch;
      });

  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      color(p.I) = boundary_marker + p.patch;
                                    });

  grid.loop_ghosts_device<0, 0, 0>(grid.nghostzones,
                                   [=] CCTK_DEVICE(const Loop::PointDesc &p)

                                       CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                         color(p.I) = ghost_marker + p.patch;
                                       });

  // Snapshot into color_pre before returning: the driver runs "SYNC: color"
  // immediately after this routine, which overwrites color's ghost zones via
  // interpatch interpolation. color_pre lives in a separate storage group
  // that nothing ever syncs, so it preserves exactly what color looked like
  // beforehand (in particular, ghost points still holding their own patch's
  // never-filled ghost_marker).
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { color_pre(p.I) = color(p.I); });
}

extern "C" void CapyrX_TestMultiPatch_write_nan_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestMultiPatch_write_nan_test;
  DECLARE_CCTK_PARAMETERS;

  const auto nan_value{std::numeric_limits<CCTK_REAL>::quiet_NaN()};

  const bool inject_ghost{CCTK_Equals(nan_injection_mode, "ghost") ||
                          CCTK_Equals(nan_injection_mode, "ghost_and_interior")};
  const bool inject_interior{
      CCTK_Equals(nan_injection_mode, "interior") ||
      CCTK_Equals(nan_injection_mode, "ghost_and_interior")};

  // Per-face flag: true where this patch's face is the genuine physical
  // outer boundary, i.e. not shared with another patch. Same call
  // write_color already makes to size its overlap band.
  std::array<CCTK_INT, 2 * Loop::dim> is_interpatch_face{};
  if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification2"))
    MultiPatch_GetBoundarySpecification2(grid.patch, 2 * Loop::dim,
                                         is_interpatch_face.data());

  Loop::vect<bool, Loop::dim> is_outer_lo, is_outer_hi;
  for (int d = 0; d < Loop::dim; ++d) {
    is_outer_lo[d] = !is_interpatch_face[2 * d + 0];
    is_outer_hi[d] = !is_interpatch_face[2 * d + 1];
  }

  const auto gsh{grid.gsh};
  const auto lbnd{grid.lbnd};
  const auto nghostzones{grid.nghostzones};

  // Baseline fill and NaN injection in a single pass: each cell's final value
  // is fully determined by its own position (either the 1.0 sentinel or
  // NaN), so there is no need to fill first and overwrite afterwards.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool inject = false;
        for (int d = 0; d < Loop::dim; ++d) {
          const auto gI{p.I[d] + lbnd[d]};
          if (is_outer_lo[d]) {
            if (inject_ghost && gI < nghostzones[d])
              inject = true;
            if (inject_interior && gI == nghostzones[d])
              inject = true;
          }
          if (is_outer_hi[d]) {
            if (inject_ghost && gI >= gsh[d] - nghostzones[d])
              inject = true;
            if (inject_interior && gI == gsh[d] - nghostzones[d] - 1)
              inject = true;
          }
        }
        nan_test(p.I) = inject ? nan_value : 1.0;
      });

  // Snapshot into nan_test_pre before returning, mirroring color_pre: "SYNC:
  // nan_test" runs immediately after this routine and overwrites nan_test's
  // ghost zones via interpatch interpolation, so nan_test_pre (never synced)
  // preserves exactly what was injected.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        nan_test_pre(p.I) = nan_test(p.I);
      });
}

extern "C" void CapyrX_TestMultiPatch_write_smooth_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_TestMultiPatch_write_smooth_test;
  DECLARE_CCTK_PARAMETERS;

  // Which analytic field to use, resolved once on the host so the device
  // lambda only carries an int (CCTK_Equals on the raw keyword string inside a
  // device kernel is neither needed nor cheap).
  SmoothKind kind;
  if (CCTK_Equals(smooth_field, "z_global"))
    kind = SmoothKind::z_global;
  else if (CCTK_Equals(smooth_field, "parabola"))
    kind = SmoothKind::parabola;
  else if (CCTK_Equals(smooth_field, "one_over_r"))
    kind = SmoothKind::one_over_r;
  else
    CCTK_VERROR("Unknown smooth_field type %s", smooth_field);

  // Fill EVERYWHERE (loop_all), including the outer-BC ghost cells: this is
  // what lets the boundary_* = none run keep exact analytic values in the
  // outer ghosts (measurement (c)). The vertex_coords field is valid in the
  // ghost zones too (write_color/write_nan_test already read it there), so the
  // exact value in an outer ghost uses that ghost's genuine global coordinate.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::sqrt;
        const auto x{vcoordx(p.I)};
        const auto y{vcoordy(p.I)};
        const auto z{vcoordz(p.I)};

        CCTK_REAL val;
        switch (kind) {
        case SmoothKind::z_global:
          val = z;
          break;
        case SmoothKind::parabola:
          val = x * x + y * y + z * z;
          break;
        case SmoothKind::one_over_r:
          val = 1.0 / sqrt(x * x + y * y + z * z);
          break;
        default:
          val = 0.0; // unreachable; keeps the compiler quiet about `val`
          break;
        }
        smooth(p.I) = val;
      });

  // Snapshot into smooth_pre before returning, mirroring nan_test_pre /
  // color_pre: "SYNC: smooth" runs immediately after this routine and
  // overwrites smooth's ghost zones (interpatch cells via interpolation,
  // outer-BC cells via the boundary kernel unless boundary_* = none), so
  // smooth_pre (never synced) preserves the exact analytic value everywhere.
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        smooth_pre(p.I) = smooth(p.I);
      });
}

extern "C" void CapyrX_TestMultiPatch_sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace CapyrX::TestMultiPatch