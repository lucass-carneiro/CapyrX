#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace TestMultiPatch {

template <typename T>
static constexpr void standing_wave(const T A, const T kx, const T ky,
                                    const T kz, const T x, const T y, const T z,
                                    T &u) {

  using std::cos;
  using std::sqrt;

  u = A * cos(2 * M_PI * kx * x) * cos(2 * M_PI * ky * y) *
      cos(2 * M_PI * kz * z);
}

extern "C" void TestMultiPatch_TestSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestSync;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        standing_wave(A, kx, ky, kz, vcoordx(p.I), vcoordy(p.I), vcoordz(p.I),
                      test_gf(p.I));
      });
}

extern "C" void TestMultiPatch_TestGhostInterp(CCTK_ARGUMENTS) {
  using std::abs;

  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestGhostInterp;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Testing interpolated values at ghost zones");

  grid.loop_all<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        const auto x{vcoordx(p.I)};
        const auto y{vcoordy(p.I)};
        const auto z{vcoordz(p.I)};
        const auto actual_gf{test_gf(p.I)};

        CCTK_REAL expected_gf{0};
        standing_wave(A, kx, ky, kz, x, y, z, expected_gf);

        // If the actual value does not match the expected and it is not 1138
        // (the boundary value), the test is a failure.
        if (!(abs(expected_gf - actual_gf) < tolerance)) {
          if (!(abs(actual_gf - 1138.0) < tolerance)) {
            CCTK_VINFO("\033[31;1mFAILED\033[0m:\n"
                       "  Local coords: (%.16f, %.16f, %.16f).\n"
                       "  Patch index: %i.\n"
                       "  Grid index (%i, %i, %i).\n"
                       "  Expected: %.16f.\n"
                       "  Obtained: %.16f",
                       x, y, z, p.patch, p.I[0], p.I[1], p.I[2], expected_gf,
                       actual_gf);
          }
        }
      });
}

// TestMultiPatch_TestCoordsSync
extern "C" void TestMultiPatch_TestCoordsSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestCoordsSync;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           test_x(p.I) = vcoordx(p.I);
                           test_y(p.I) = vcoordy(p.I);
                           test_z(p.I) = vcoordz(p.I);
                         });
}

extern "C" void TestMultiPatch_TestCoordsGhostInterp(CCTK_ARGUMENTS) {
  using std::fabs;

  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestCoordsGhostInterp;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Testing interpolated coordinate values at ghost zones");

  grid.loop_all<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        const auto expected_x{vcoordx(p.I)};
        const auto expected_y{vcoordy(p.I)};
        const auto expected_z{vcoordz(p.I)};

        const auto obtained_x{test_x(p.I)};
        const auto obtained_y{test_y(p.I)};
        const auto obtained_z{test_z(p.I)};

        if (!(fabs(expected_x - obtained_x) < exact_tolerance)) {
          if (!(fabs(obtained_x - 1138.0) < exact_tolerance)) {
            CCTK_VINFO("\033[31;1mFAILED\033[0m:\n"
                       "  Local coords: (%.16f, %.16f, %.16f).\n"
                       "  Patch index: %i.\n"
                       "  Grid index (%i, %i, %i).\n"
                       "  Expected x: %.16f.\n"
                       "  Obtained x: %.16f.\n"
                       "  Error: %.16f\n",
                       expected_x, expected_y, expected_z, p.patch, p.I[0],
                       p.I[1], p.I[2], expected_x, obtained_x,
                       expected_x - obtained_x);
          }
        }

        if (!(abs(expected_y - obtained_y) < exact_tolerance)) {
          if (!(abs(obtained_y - 1138.0) < exact_tolerance)) {
            CCTK_VINFO("\033[31;1mFAILED\033[0m:\n"
                       "  Local coords: (%.16f, %.16f, %.16f).\n"
                       "  Patch index: %i.\n"
                       "  Grid index (%i, %i, %i).\n"
                       "  Expected y: %.16f.\n"
                       "  Obtained y: %.16f"
                       "  Error: %.16f\n",
                       expected_x, expected_y, expected_z, p.patch, p.I[0],
                       p.I[1], p.I[2], expected_y, obtained_y,
                       expected_y - obtained_y);
          }
        }

        if (!(abs(expected_z - obtained_z) < exact_tolerance)) {
          if (!(abs(obtained_z - 1138.0) < exact_tolerance)) {
            CCTK_VINFO("\033[31;1mFAILED\033[0m:\n"
                       "  Local coords: (%.16f, %.16f, %.16f).\n"
                       "  Patch index: %i.\n"
                       "  Grid index (%i, %i, %i).\n"
                       "  Expected z: %.16f.\n"
                       "  Obtained z: %.16f"
                       "  Error: %.16f\n",
                       expected_x, expected_y, expected_z, p.patch, p.I[0],
                       p.I[1], p.I[2], expected_z, obtained_z,
                       expected_z - obtained_z);
          }
        }
      });
}

} // namespace TestMultiPatch
