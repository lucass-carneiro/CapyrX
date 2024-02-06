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

  u = A * cos(2 * M_PI * kx * x) * cos(2 * M_PI * ky * y) *
      cos(2 * M_PI * kz * z);
}

extern "C" void TestMultiPatch_TestSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestSync;
  DECLARE_CCTK_PARAMETERS;

  const auto loop_lambda =
      [=] ARITH_DEVICE(const Loop::PointDesc &p) ARITH_INLINE {
        standing_wave(A, kx, ky, kz, vcoordx(p.I), vcoordy(p.I), vcoordz(p.I),
                      test_gf(p.I));
      };

  grid.loop_int_device<0, 0, 0>(grid.nghostzones, loop_lambda);
}

extern "C" void TestMultiPatch_TestGhostInterp(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_TestGhostInterp;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Testing interpolated coordinate values at ghost zones");

  /*
   * 1. For each local point + patch number, compute the global point using the
   * aliased function MultiPatch_LocalToGlobal2
   * 2. Compare the result with the value stored in vcoord(xyz)(p.I)
   * 3. If these values do not match, report error
   */
  const auto loop_lambda = [=] ARITH_DEVICE(
                               const Loop::PointDesc &p) ARITH_INLINE {
    const CCTK_INT patches[1] = {p.patch};

    const CCTK_REAL localsx[1] = {p.x};
    const CCTK_REAL localsy[1] = {p.y};
    const CCTK_REAL localsz[1] = {p.z};

    CCTK_REAL globalsx[1] = {0};
    CCTK_REAL globalsy[1] = {0};
    CCTK_REAL globalsz[1] = {0};

    MultiPatch_LocalToGlobal2(1, patches, localsx, localsy, localsz, globalsx,
                              globalsy, globalsz);

    if (!MultiPatchTests::isapprox(vcoordx(p.I), globalsx[0])) {
      CCTK_VINFO(
          "Test FAILED at x direction, patch %i. Got: %.16f. Expected %.16f",
          p.patch, vcoordx(p.I), globalsx[0]);
    }

    if (!MultiPatchTests::isapprox(vcoordy(p.I), globalsy[0])) {
      CCTK_VINFO(
          "Test FAILED at x direction, patch %i. Got: %.16f. Expected %.16f",
          p.patch, vcoordy(p.I), globalsy[0]);
    }

    if (!MultiPatchTests::isapprox(vcoordz(p.I), globalsz[0])) {
      CCTK_VINFO(
          "Test FAILED at x direction, patch %i. Got: %.16f. Expected %.16f",
          p.patch, vcoordz(p.I), globalsz[0]);
    }
  };

  grid.loop_all_device<0, 0, 0>(grid.nghostzones, loop_lambda);
}

} // namespace TestMultiPatch
