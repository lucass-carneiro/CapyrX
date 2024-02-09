#include <loop.hxx>
#include <loop_device.hxx>

#include <arith.hxx>
#include <vec.hxx>

#include "../../MultiPatch/src/tests.hxx"

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

  // u = A * cos(2 * M_PI * kx * x) * cos(2 * M_PI * ky * y) * cos(2 * M_PI * kz
  // * z);
  u = T{0};
}

extern "C" void TestMultiPatch_TestSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestSync;
  DECLARE_CCTK_PARAMETERS;

  const auto loop_lambda = [=](const Loop::PointDesc &p) ARITH_INLINE {
    standing_wave(A, kx, ky, kz, vcoordx(p.I), vcoordy(p.I), vcoordz(p.I),
                  test_gf(p.I));
  };

  grid.loop_int<0, 0, 0>(grid.nghostzones, loop_lambda);
}

extern "C" void TestMultiPatch_TestGhostInterp(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestGhostInterp;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Testing interpolated values at ghost zones");

  const auto loop_lambda = [=](const Loop::PointDesc &p) ARITH_INLINE {
    const auto x{vcoordx(p.I)};
    const auto y{vcoordy(p.I)};
    const auto z{vcoordz(p.I)};
    const auto actual_gf{test_gf(p.I)};

    CCTK_REAL expected_gf{0};
    standing_wave(A, kx, ky, kz, x, y, z, expected_gf);

    if (!MultiPatchTests::isapprox(actual_gf, expected_gf)) {
      CCTK_VINFO("Test FAILED at (%.16f, %.16f, %.16f). Patch index %i. Grid "
                 "index (%i, %i, %i). "
                 "Expected: %.16f. "
                 "Obtained %.16f",
                 x, y, z, p.patch, p.I[0], p.I[1], p.I[2], expected_gf,
                 actual_gf);
    }
  };

  grid.loop_all<0, 0, 0>(grid.nghostzones, loop_lambda);
}

} // namespace TestMultiPatch
