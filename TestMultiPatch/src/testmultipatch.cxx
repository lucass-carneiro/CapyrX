#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace TestMultiPatch {

static auto standing_wave(CCTK_REAL A, CCTK_REAL kx, CCTK_REAL ky, CCTK_REAL kz,
                          CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                          CCTK_REAL z) noexcept -> CCTK_REAL {
  using std::cos, std::sin, std::sqrt;

  const auto pi{acos(-1.0)};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};

  return A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
         cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

extern "C" void TestMultiPatch_write_state(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_write_state;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(grid.nghostzones,
                         [=](const Loop::PointDesc &p) ARITH_INLINE {
                           const auto t{cctk_time};
                           const auto x{vcoordx(p.I)};
                           const auto y{vcoordy(p.I)};
                           const auto z{vcoordz(p.I)};

                           u(p.I) = standing_wave(A, kx, ky, kz, t, x, y, z);
                         });
}

extern "C" void TestMultiPatch_write_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestMultiPatch_write_error;
  DECLARE_CCTK_PARAMETERS;

  using std::fabs;

  grid.loop_all<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        const auto t{cctk_time};
        const auto x{vcoordx(p.I)};
        const auto y{vcoordy(p.I)};
        const auto z{vcoordz(p.I)};

        const auto evolved_u{u(p.I)};
        const auto real_u{standing_wave(A, kx, ky, kz, t, x, y, z)};

        u_err(p.I) = fabs(evolved_u - real_u);
      });
}

extern "C" void TestMultiPatch_sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace TestMultiPatch
