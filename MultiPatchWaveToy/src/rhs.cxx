#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstddef>
#include <loop_device.hxx>
#include <global_derivatives.hxx>

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
        const FirstDerivativeInputs local_dPi{c4o<fd_dir::a>(p, Pi),
                                              c4o<fd_dir::b>(p, Pi),
                                              c4o<fd_dir::c>(p, Pi)};
        const auto dPi{project_first(local_dPi, p, VERTEX_JACOBIANS_FIRST)};

        const FirstDerivativeInputs local_dDx{c4o<fd_dir::a>(p, Dx),
                                              c4o<fd_dir::b>(p, Dx),
                                              c4o<fd_dir::c>(p, Dx)};
        const auto dDx{project_first(local_dDx, p, VERTEX_JACOBIANS_FIRST)};

        const FirstDerivativeInputs local_dDy{c4o<fd_dir::a>(p, Dy),
                                              c4o<fd_dir::b>(p, Dy),
                                              c4o<fd_dir::c>(p, Dy)};
        const auto dDy{project_first(local_dDy, p, VERTEX_JACOBIANS_FIRST)};

        const FirstDerivativeInputs local_dDz{c4o<fd_dir::a>(p, Dz),
                                              c4o<fd_dir::b>(p, Dz),
                                              c4o<fd_dir::c>(p, Dz)};
        const auto dDz{project_first(local_dDz, p, VERTEX_JACOBIANS_FIRST)};

        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = dDx.dx + dDy.dy + dDz.dz;
        Dx_rhs(p.I) = dPi.dx;
        Dy_rhs(p.I) = dPi.dy;
        Dz_rhs(p.I) = dPi.dz;
      });
}

extern "C" void MultiPatchWaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace MultiPatchWaveToy
