#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstddef>
#include <loop_device.hxx>
#include <global_derivatives.hxx>
#include <local_derivatives.hxx>

namespace MultiPatchWaveToy {

enum class fd_dir : std::size_t { a = 0, b = 1, c = 2 };

template <fd_dir dir>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o(const Loop::PointDesc &p,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  const auto d{static_cast<std::size_t>(dir)};
  const auto num{gf(-2 * p.DI[d] + p.I) - 8 * gf(-p.DI[d] + p.I) +
                 8 * gf(p.DI[d] + p.I) - gf(2 * p.DI[d] + p.I)};
  const auto den{1.0 / (12 * p.DX[d])};
  return num * den;
}

extern "C" void MultiPatchWaveToy_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  using namespace MultiPatch::GlobalDerivatives;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const LocalFirstDerivatives lPi{c4o_1_0_0(p, Pi), c4o_0_1_0(p, Pi),
                                        c4o_0_0_1(p, Pi)};

        const LocalFirstDerivatives lDx{c4o_1_0_0(p, Dx), c4o_0_1_0(p, Dx),
                                        c4o_0_0_1(p, Dx)};

        const LocalFirstDerivatives lDy{c4o_1_0_0(p, Dy), c4o_0_1_0(p, Dy),
                                        c4o_0_0_1(p, Dy)};

        const LocalFirstDerivatives lDz{c4o_1_0_0(p, Dz), c4o_0_1_0(p, Dz),
                                        c4o_0_0_1(p, Dz)};

        const Jacobians jac{VERTEX_JACOBIANS_FIRST(p)};

        const auto gPi{project_first(lPi, jac)};
        const auto gDx{project_first(lDx, jac)};
        const auto gDy{project_first(lDy, jac)};
        const auto gDz{project_first(lDz, jac)};

        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = gDx.dx + gDy.dy + gDz.dz;
        Dx_rhs(p.I) = gPi.dx;
        Dy_rhs(p.I) = gPi.dy;
        Dz_rhs(p.I) = gPi.dz;
      });
}

extern "C" void MultiPatchWaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace MultiPatchWaveToy
