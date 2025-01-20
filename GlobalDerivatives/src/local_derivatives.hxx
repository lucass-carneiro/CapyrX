#ifndef CAPYRX_LOCAL_DERIVATIVES_HXX
#define CAPYRX_LOCAL_DERIVATIVES_HXX

#include <cctk.h>
#include <loop_device.hxx>

namespace MultiPatch::GlobalDerivatives {

static inline auto
c4o_1_0_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[0] + p.I) - 8 * gf(-p.DI[0] + p.I) +
                 8 * gf(p.DI[0] + p.I) - gf(2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * p.DX[0])};
  return num * den;
}

static inline auto
c4o_0_1_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[1] + p.I) - 8 * gf(-p.DI[1] + p.I) +
                 8 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * p.DX[1])};
  return num * den;
}

static inline auto
c4o_0_0_1(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{gf(-2 * p.DI[2] + p.I) - 8 * gf(-p.DI[2] + p.I) +
                 8 * gf(p.DI[2] + p.I) - gf(2 * p.DI[2] + p.I)};
  const auto den{1.0 / (12 * p.DX[2])};
  return num * den;
}

static inline auto
c4o_2_0_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{-30 * gf(p.I) - gf(-2 * p.DI[0] + p.I) +
                 16 * gf(-p.DI[0] + p.I) + 16 * gf(p.DI[0] + p.I) -
                 gf(2 * p.DI[0] + p.I)};
  const auto den{1.0 / (12 * (p.DX[0] * p.DX[0]))};
  return num * den;
}

static inline auto
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

static inline auto
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

static inline auto
c4o_0_2_0(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{-30 * gf(p.I) - gf(-2 * p.DI[1] + p.I) +
                 16 * gf(-p.DI[1] + p.I) + 16 * gf(p.DI[1] + p.I) -
                 gf(2 * p.DI[1] + p.I)};
  const auto den{1.0 / (12 * (p.DX[1] * p.DX[1]))};
  return num * den;
}

static inline auto
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

static inline auto
c4o_0_0_2(const Loop::PointDesc &p,
          const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto num{-30 * gf(p.I) - gf(-2 * p.DI[2] + p.I) +
                 16 * gf(-p.DI[2] + p.I) + 16 * gf(p.DI[2] + p.I) -
                 gf(2 * p.DI[2] + p.I)};
  const auto den{1.0 / (12 * (p.DX[2] * p.DX[2]))};
  return num * den;
}

} // namespace MultiPatch::GlobalDerivatives

#endif // CAPYRX_LOCAL_DERIVATIVES_HXX