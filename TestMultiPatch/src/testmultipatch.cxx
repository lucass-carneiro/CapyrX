#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace TestMultiPatch {

/*
 * u(t,x,y,z) = A cos(2 pi omega t) sin(2 pi kx x) sin(2 pi ky y) sin(2 pi kz z)
 */
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

        if (!(abs(expected_gf - actual_gf) < tolerance)) {
          CCTK_VINFO("Test FAILED:\n"
                     "  Local coords: (%.16f, %.16f, %.16f).\n"
                     "  Patch index: %i.\n"
                     "  Grid index (%i, %i, %i).\n"
                     "  Expected: %.16f.\n"
                     "  Obtained %.16f",
                     x, y, z, p.patch, p.I[0], p.I[1], p.I[2], expected_gf,
                     actual_gf);
        }
      });
}

} // namespace TestMultiPatch
