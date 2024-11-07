#include "global_derivatives.hxx"

#include <array>

/*
CCTK_REAL vertex_dJacobians TYPE=gf CENTERING={vvv} TAGS='checkpoint="no"'
{
  vdJ_d2a_dxdx, vdJ_d2a_dxdy, vdJ_d2a_dxdz, vdJ_d2a_dydy, vdJ_d2a_dydz,
vdJ_d2a_dzdz, vdJ_d2b_dxdx, vdJ_d2b_dxdy, vdJ_d2b_dxdz, vdJ_d2b_dydy,
vdJ_d2b_dydz, vdJ_d2b_dzdz, vdJ_d2c_dxdx, vdJ_d2c_dxdy, vdJ_d2c_dxdz,
vdJ_d2c_dydy, vdJ_d2c_dydz, vdJ_d2c_dzdz } "The vertex centered Jacobian matrix
derivatives of the patch coordinate transformations"
*/

namespace MultiPatch::GlobalDerivatives {

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
get_var(const CCTK_POINTER_TO_CONST cctkGH_,
        const char *var_name) noexcept -> Loop::GF3D2<const CCTK_REAL> {
  using namespace Loop;

  const auto cctkGH{static_cast<const cGH *>(cctkGH_)};
  const std::array<int, dim> c{0, 0, 0};
  const GF3D2layout l(cctkGH, c);

  const auto var_ptr{CCTK_VarDataPtr(cctkGH, 0, var_name)};
  return GF3D2<const CCTK_REAL>(l, static_cast<const CCTK_REAL *>(var_ptr));
}

// Finite difference directions
enum class local_fd_dir : std::size_t { a = 0, b = 1, c = 2 };

// Local central first derivative, fourth order accurae operator
template <local_fd_dir direction>
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c_l_1_4(const Loop::PointDesc &p,
        const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> CCTK_REAL {
  constexpr auto d{static_cast<size_t>(direction)};

  const auto num{gf(p.I - 2 * p.DI[d]) - 8.0 * gf(p.I - 1 * p.DI[d]) +
                 8.0 * gf(p.I + 1 * p.DI[d]) - 1.0 * gf(p.I + 2 * p.DI[d])};
  const auto den{1.0 / (12.0 * p.DX[d])};

  return den * num;
}

// Global Central first derivative, fourth order accurae operator
auto CCTK_HOST CCTK_DEVICE c_1_4(
    const CCTK_POINTER_TO_CONST cctkGH_, const Loop::PointDesc &p,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> FirstDerivativeResults {
  using namespace Loop;

  // Recover Jacobians
  const auto da_dx{get_var(cctkGH_, "MultiPatch::vJ_da_dx")(p.I)};
  const auto da_dy{get_var(cctkGH_, "MultiPatch::vJ_da_dy")(p.I)};
  const auto da_dz{get_var(cctkGH_, "MultiPatch::vJ_da_dz")(p.I)};
  const auto db_dx{get_var(cctkGH_, "MultiPatch::vJ_db_dx")(p.I)};
  const auto db_dy{get_var(cctkGH_, "MultiPatch::vJ_db_dy")(p.I)};
  const auto db_dz{get_var(cctkGH_, "MultiPatch::vJ_db_dz")(p.I)};
  const auto dc_dx{get_var(cctkGH_, "MultiPatch::vJ_dc_dx")(p.I)};
  const auto dc_dy{get_var(cctkGH_, "MultiPatch::vJ_dc_dy")(p.I)};
  const auto dc_dz{get_var(cctkGH_, "MultiPatch::vJ_dc_dz")(p.I)};

  // Compute local derivatives
  const auto dgf_da{c_l_1_4<local_fd_dir::a>(p, gf)};
  const auto dgf_db{c_l_1_4<local_fd_dir::a>(p, gf)};
  const auto dgf_dc{c_l_1_4<local_fd_dir::a>(p, gf)};

  // Projections
  const auto dgf_dx{dgf_db * db_dx + dgf_dc * dc_dx + da_dx * dgf_da};
  const auto dgf_dy{dgf_dc * dc_dy + db_dy * dgf_db + da_dy * dgf_da};
  const auto dgf_dz{dc_dz * dgf_dc + db_dz * dgf_db + da_dz * dgf_da};

  return FirstDerivativeResults{dgf_dx, dgf_dy, dgf_dz};
}

} // namespace MultiPatch::GlobalDerivatives