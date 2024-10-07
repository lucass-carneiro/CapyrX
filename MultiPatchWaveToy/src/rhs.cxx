// clang-format off
#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
// clang-format on

namespace MultiPatchWaveToy {

// Finite difference directions
enum class local_fd_dir : std::size_t { a = 0, b = 1, c = 2 };

// Local central 1st order derivative (4th order accurate)
template <local_fd_dir direction, typename T>
static inline auto fd_l_1_4(const Loop::PointDesc &p,
                            const Loop::GF3D2<T> &gf) noexcept -> T {
  constexpr auto d{static_cast<size_t>(direction)};

  const auto num{gf(p.I - 2 * p.DI[d]) - 8.0 * gf(p.I - 1 * p.DI[d]) +
                 8.0 * gf(p.I + 1 * p.DI[d]) - 1.0 * gf(p.I + 2 * p.DI[d])};
  const auto den{1.0 / (12.0 * p.DX[d])};

  return den * num;
}

extern "C" void MultiPatchWaveToy_RHS(CCTK_ARGUMENTS) {
  using std::pow;

  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<0, 0, 0>(
      grid.nghostzones, [=](const Loop::PointDesc &p) ARITH_INLINE {
        const auto d_Pi_dx{vJ_da_dx(p.I) * fd_l_1_4<local_fd_dir::a>(p, Pi) +
                           vJ_db_dx(p.I) * fd_l_1_4<local_fd_dir::b>(p, Pi) +
                           vJ_dc_dx(p.I) * fd_l_1_4<local_fd_dir::c>(p, Pi)};

        const auto d_Pi_dy{vJ_da_dy(p.I) * fd_l_1_4<local_fd_dir::a>(p, Pi) +
                           vJ_db_dy(p.I) * fd_l_1_4<local_fd_dir::b>(p, Pi) +
                           vJ_dc_dy(p.I) * fd_l_1_4<local_fd_dir::c>(p, Pi)};

        const auto d_Pi_dz{vJ_da_dz(p.I) * fd_l_1_4<local_fd_dir::a>(p, Pi) +
                           vJ_db_dz(p.I) * fd_l_1_4<local_fd_dir::b>(p, Pi) +
                           vJ_dc_dz(p.I) * fd_l_1_4<local_fd_dir::c>(p, Pi)};

        const auto d_Dx_dx{vJ_da_dx(p.I) * fd_l_1_4<local_fd_dir::a>(p, Dx) +
                           vJ_db_dx(p.I) * fd_l_1_4<local_fd_dir::b>(p, Dx) +
                           vJ_dc_dx(p.I) * fd_l_1_4<local_fd_dir::c>(p, Dx)};

        const auto d_Dy_dy{vJ_da_dy(p.I) * fd_l_1_4<local_fd_dir::a>(p, Dy) +
                           vJ_db_dy(p.I) * fd_l_1_4<local_fd_dir::b>(p, Dy) +
                           vJ_dc_dy(p.I) * fd_l_1_4<local_fd_dir::c>(p, Dy)};

        const auto d_Dz_dz{vJ_da_dz(p.I) * fd_l_1_4<local_fd_dir::a>(p, Dz) +
                           vJ_db_dz(p.I) * fd_l_1_4<local_fd_dir::b>(p, Dz) +
                           vJ_dc_dz(p.I) * fd_l_1_4<local_fd_dir::c>(p, Dz)};

        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = d_Dx_dx + d_Dy_dy + d_Dz_dz;
        Dx_rhs(p.I) = d_Pi_dx;
        Dy_rhs(p.I) = d_Pi_dy;
        Dz_rhs(p.I) = d_Pi_dz;
      });
}

extern "C" void MultiPatchWaveToy_Sync(CCTK_ARGUMENTS) {
  // Do nothing
}

} // namespace MultiPatchWaveToy
