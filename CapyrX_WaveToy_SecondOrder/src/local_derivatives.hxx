#ifndef CAPYRX_WAVETOY_SECOND_ORDER_LOCAL_DERIVATIVES_HXX
#define CAPYRX_WAVETOY_SECOND_ORDER_LOCAL_DERIVATIVES_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <cmath>

namespace CapyrX::WaveToy_SecondOrder {

template <std::size_t dir>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
d_x(const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const CCTK_REAL invdx = 1.0 / p.DX[dir];

  const auto Im1dx = p.I - 1 * p.DI[dir];
  const auto Im2dx = p.I - 2 * p.DI[dir];
  const auto Ip1dx = p.I + 1 * p.DI[dir];
  const auto Ip2dx = p.I + 2 * p.DI[dir];

  return invdx / 12.0 *
         (-8 * gf(Im1dx) + gf(Im2dx) + 8 * gf(Ip1dx) - gf(Ip2dx));
}

template <std::size_t dir>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
d_xx(const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const CCTK_REAL invdxsq = 1.0 / (p.DX[dir] * p.DX[dir]);
  const auto Im2dx = p.I - 2 * p.DI[dir];
  const auto Im1dx = p.I - 1 * p.DI[dir];
  const auto Ip1dx = p.I + 1 * p.DI[dir];
  const auto Ip2dx = p.I + 2 * p.DI[dir];
  const auto Ip3dx = p.I + 3 * p.DI[dir];
  const auto Ip4dx = p.I + 4 * p.DI[dir];

  return invdxsq / 12.0 *
         (16 * gf(Im1dx) - gf(Im2dx) + 16 * gf(Ip1dx) - gf(Ip2dx) -
          30 * gf(p.I));
}

template <std::size_t dir1, std::size_t dir2>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
d_xy(const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const auto invdxdy = 1.0 / (p.DX[dir1] * p.DX[dir2]);

  const auto Ip2dxp2dy = p.I + 2 * p.DI[dir1] + 2 * p.DI[dir2];
  const auto Ip2dxp1dy = p.I + 2 * p.DI[dir1] + 1 * p.DI[dir2];
  const auto Ip2dxm2dy = p.I + 2 * p.DI[dir1] - 2 * p.DI[dir2];
  const auto Ip2dxm1dy = p.I + 2 * p.DI[dir1] - 1 * p.DI[dir2];

  const auto Ip1dxp2dy = p.I + 1 * p.DI[dir1] + 2 * p.DI[dir2];
  const auto Ip1dxp1dy = p.I + 1 * p.DI[dir1] + 1 * p.DI[dir2];
  const auto Ip1dxm2dy = p.I + 1 * p.DI[dir1] - 2 * p.DI[dir2];
  const auto Ip1dxm1dy = p.I + 1 * p.DI[dir1] - 1 * p.DI[dir2];

  const auto Im2dxp2dy = p.I - 2 * p.DI[dir1] + 2 * p.DI[dir2];
  const auto Im2dxp1dy = p.I - 2 * p.DI[dir1] + 1 * p.DI[dir2];
  const auto Im2dxm2dy = p.I - 2 * p.DI[dir1] - 2 * p.DI[dir2];
  const auto Im2dxm1dy = p.I - 2 * p.DI[dir1] - 1 * p.DI[dir2];

  const auto Im1dxp2dy = p.I - 1 * p.DI[dir1] + 2 * p.DI[dir2];
  const auto Im1dxp1dy = p.I - 1 * p.DI[dir1] + 1 * p.DI[dir2];
  const auto Im1dxm2dy = p.I - 1 * p.DI[dir1] - 2 * p.DI[dir2];
  const auto Im1dxm1dy = p.I - 1 * p.DI[dir1] - 1 * p.DI[dir2];

  return invdxdy / 144.0 *
         (64 * gf(Im1dxm1dy) - 8 * gf(Im1dxm2dy) - 64 * gf(Im1dxp1dy) +
          8 * gf(Im1dxp2dy) - 8 * gf(Im2dxm1dy) + gf(Im2dxm2dy) +
          8 * gf(Im2dxp1dy) - gf(Im2dxp2dy) - 64 * gf(Ip1dxm1dy) +
          8 * gf(Ip1dxm2dy) + 64 * gf(Ip1dxp1dy) - 8 * gf(Ip1dxp2dy) +
          8 * gf(Ip2dxm1dy) - gf(Ip2dxm2dy) - 8 * gf(Ip2dxp1dy) +
          gf(Ip2dxp2dy));
}

template <std::size_t dir>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE diss_5(
    const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> CCTK_REAL {
  const auto fac{(1.0 / 64.0) * (1.0 / p.DX[dir])};
  const auto stencil{gf(p.I - 3 * p.DI[dir]) - 6.0 * gf(p.I - 2 * p.DI[dir]) +
                     15.0 * gf(p.I - 1 * p.DI[dir]) - 20.0 * gf(p.I) +
                     15.0 * gf(p.I + 1 * p.DI[dir]) -
                     6.0 * gf(p.I + 2 * p.DI[dir]) + gf(p.I + 3 * p.DI[dir])};
  return fac * stencil;
}

} // namespace CapyrX::WaveToy_SecondOrder

#endif // CAPYRX_WAVETOY_SECOND_ORDER_LOCAL_DERIVATIVES_HXX