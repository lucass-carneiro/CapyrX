#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <global_derivatives.hxx>

#include <cmath>

namespace TestMultiPatch {

static inline auto standing_wave(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                 CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                 CCTK_REAL y,
                                 CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
         cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

static inline auto standing_wave_dx(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                    CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                    CCTK_REAL y,
                                    CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -2 * A * kx * pi * cos(2 * omega * pi * t) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * kx * pi * x);
}

static inline auto standing_wave_dy(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                    CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                    CCTK_REAL y,
                                    CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -2 * A * ky * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * kz * pi * z) * sin(2 * ky * pi * y);
}

static inline auto standing_wave_dz(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                    CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                    CCTK_REAL y,
                                    CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -2 * A * kz * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

static inline auto standing_wave_dx2(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                     CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                     CCTK_REAL y,
                                     CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -4 * A * kx * kx * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

static inline auto standing_wave_dy2(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                     CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                     CCTK_REAL y,
                                     CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -4 * A * ky * ky * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

static inline auto standing_wave_dz2(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                     CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                     CCTK_REAL y,
                                     CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return -4 * A * kz * kz * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

static inline auto standing_wave_dxy(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                     CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                     CCTK_REAL y,
                                     CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return 4 * A * kx * ky * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kz * pi * z) * sin(2 * kx * pi * x) * sin(2 * ky * pi * y);
}

static inline auto standing_wave_dxz(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                     CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                     CCTK_REAL y,
                                     CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return 4 * A * kx * kz * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * ky * pi * y) * sin(2 * kx * pi * x) * sin(2 * kz * pi * z);
}

static inline auto standing_wave_dyz(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky,
                                     CCTK_REAL kz, CCTK_REAL t, CCTK_REAL x,
                                     CCTK_REAL y,
                                     CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return 4 * A * ky * kz * pi * pi * cos(2 * omega * pi * t) *
         cos(2 * kx * pi * x) * sin(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

static inline auto parabola(CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) -> CCTK_REAL {
  return x * x + y * y + z * z;
}

static inline auto parabola_dx(CCTK_REAL x, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0 * x;
}

static inline auto parabola_dy(CCTK_REAL, CCTK_REAL y, CCTK_REAL) -> CCTK_REAL {
  return 2.0 * y;
}

static inline auto parabola_dz(CCTK_REAL, CCTK_REAL, CCTK_REAL z) -> CCTK_REAL {
  return 2.0 * z;
}

static inline auto parabola_dx2(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0;
}

static inline auto parabola_dy2(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0;
}

static inline auto parabola_dz2(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 2.0;
}

static inline auto parabola_dxy(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 0.0;
}

static inline auto parabola_dxz(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 0.0;
}

static inline auto parabola_dyz(CCTK_REAL, CCTK_REAL, CCTK_REAL) -> CCTK_REAL {
  return 0.0;
}

extern "C" void TestMultiPatch_write_test_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_write_test_data;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(test_data, "standing wave")) {
    grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          u(p.I) = standing_wave(A, kx, ky, kz, t, x, y, z);
        });

  } else if (CCTK_Equals(test_data, "parabola")) {
    grid.loop_int<0, 0, 0>(grid.nghostzones,
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

extern "C" void TestMultiPatch_compute_interp_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_compute_interp_error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(test_data, "standing wave")) {
    grid.loop_all<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};
          const auto x{vcoordx(p.I)};
          const auto y{vcoordy(p.I)};
          const auto z{vcoordz(p.I)};

          const auto evolved_u{u(p.I)};
          const auto real_u{standing_wave(A, kx, ky, kz, t, x, y, z)};

          interp(p.I) = fabs(evolved_u - real_u);
        });

  } else if (CCTK_Equals(test_data, "parabola")) {
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [=] CCTK_DEVICE(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 using std::fabs;

                                 const auto x{vcoordx(p.I)};
                                 const auto y{vcoordy(p.I)};
                                 const auto z{vcoordz(p.I)};

                                 const auto evolved_u{u(p.I)};
                                 const auto real_u{parabola(x, y, z)};

                                 interp(p.I) = fabs(evolved_u - real_u);
                               });
  }
}

extern "C" void TestMultiPatch_compute_deriv_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_compute_deriv_error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_Equals(test_data, "standing wave")) {
    grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using namespace MultiPatch::GlobalDerivatives;

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

          // TODO Get from projections
          const auto first_derivs{c4o_1(cctkGH, p, u)};
          const auto second_derivs{c4o_2(cctkGH, p, u)};

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
    grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using namespace MultiPatch::GlobalDerivatives;

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

          // TODO: Get from projections
          const auto first_derivs{c4o_1(cctkGH, p, u)};
          const auto second_derivs{c4o_2(cctkGH, p, u)};

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

extern "C" void TestMultiPatch_sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace TestMultiPatch
