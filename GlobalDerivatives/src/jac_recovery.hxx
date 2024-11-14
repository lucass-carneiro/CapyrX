#ifndef CAPYRX_GLOBAL_DERIVATIVES_JAC_RECOVERY_HXX
#define CAPYRX_GLOBAL_DERIVATIVES_JAC_RECOVERY_HXX

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

struct JacobianComponents {
  CCTK_REAL da_dx{};
  CCTK_REAL da_dy{};
  CCTK_REAL da_dz{};
  CCTK_REAL db_dx{};
  CCTK_REAL db_dy{};
  CCTK_REAL db_dz{};
  CCTK_REAL dc_dx{};
  CCTK_REAL dc_dy{};
  CCTK_REAL dc_dz{};

  CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE inline JacobianComponents(
      const CCTK_POINTER_TO_CONST cctkGH_, const Loop::PointDesc &p) noexcept {
    da_dx = get_var(cctkGH_, "MultiPatch::vJ_da_dx")(p.I);
    da_dy = get_var(cctkGH_, "MultiPatch::vJ_da_dy")(p.I);
    da_dz = get_var(cctkGH_, "MultiPatch::vJ_da_dz")(p.I);
    db_dx = get_var(cctkGH_, "MultiPatch::vJ_db_dx")(p.I);
    db_dy = get_var(cctkGH_, "MultiPatch::vJ_db_dy")(p.I);
    db_dz = get_var(cctkGH_, "MultiPatch::vJ_db_dz")(p.I);
    dc_dx = get_var(cctkGH_, "MultiPatch::vJ_dc_dx")(p.I);
    dc_dy = get_var(cctkGH_, "MultiPatch::vJ_dc_dy")(p.I);
    dc_dz = get_var(cctkGH_, "MultiPatch::vJ_dc_dz")(p.I);
  }
};

struct JacobianDerivativeComponents {
  CCTK_REAL d2a_dxdx{};
  CCTK_REAL d2a_dxdy{};
  CCTK_REAL d2a_dxdz{};
  CCTK_REAL d2a_dydy{};
  CCTK_REAL d2a_dydz{};
  CCTK_REAL d2a_dzdz{};
  CCTK_REAL d2b_dxdx{};
  CCTK_REAL d2b_dxdy{};
  CCTK_REAL d2b_dxdz{};
  CCTK_REAL d2b_dydy{};
  CCTK_REAL d2b_dydz{};
  CCTK_REAL d2b_dzdz{};
  CCTK_REAL d2c_dxdx{};
  CCTK_REAL d2c_dxdy{};
  CCTK_REAL d2c_dxdz{};
  CCTK_REAL d2c_dydy{};
  CCTK_REAL d2c_dydz{};
  CCTK_REAL d2c_dzdz{};

  CCTK_ATTRIBUTE_ALWAYS_INLINE
      CCTK_HOST CCTK_DEVICE inline JacobianDerivativeComponents(
          const CCTK_POINTER_TO_CONST cctkGH_,
          const Loop::PointDesc &p) noexcept {
    d2a_dxdx = get_var(cctkGH_, "MultiPatch::vdJ_d2a_dxdx")(p.I);
    d2a_dxdy = get_var(cctkGH_, "MultiPatch::vdJ_d2a_dxdy")(p.I);
    d2a_dxdz = get_var(cctkGH_, "MultiPatch::vdJ_d2a_dxdz")(p.I);
    d2a_dydy = get_var(cctkGH_, "MultiPatch::vdJ_d2a_dydy")(p.I);
    d2a_dydz = get_var(cctkGH_, "MultiPatch::vdJ_d2a_dydz")(p.I);
    d2a_dzdz = get_var(cctkGH_, "MultiPatch::vdJ_d2a_dzdz")(p.I);
    d2b_dxdx = get_var(cctkGH_, "MultiPatch::vdJ_d2b_dxdx")(p.I);
    d2b_dxdy = get_var(cctkGH_, "MultiPatch::vdJ_d2b_dxdy")(p.I);
    d2b_dxdz = get_var(cctkGH_, "MultiPatch::vdJ_d2b_dxdz")(p.I);
    d2b_dydy = get_var(cctkGH_, "MultiPatch::vdJ_d2b_dydy")(p.I);
    d2b_dydz = get_var(cctkGH_, "MultiPatch::vdJ_d2b_dydz")(p.I);
    d2b_dzdz = get_var(cctkGH_, "MultiPatch::vdJ_d2b_dzdz")(p.I);
    d2c_dxdx = get_var(cctkGH_, "MultiPatch::vdJ_d2c_dxdx")(p.I);
    d2c_dxdy = get_var(cctkGH_, "MultiPatch::vdJ_d2c_dxdy")(p.I);
    d2c_dxdz = get_var(cctkGH_, "MultiPatch::vdJ_d2c_dxdz")(p.I);
    d2c_dydy = get_var(cctkGH_, "MultiPatch::vdJ_d2c_dydy")(p.I);
    d2c_dydz = get_var(cctkGH_, "MultiPatch::vdJ_d2c_dydz")(p.I);
    d2c_dzdz = get_var(cctkGH_, "MultiPatch::vdJ_d2c_dzdz")(p.I);
  }
};

} // namespace MultiPatch::GlobalDerivatives

#endif // CAPYRX_GLOBAL_DERIVATIVES_JAC_RECOVERY_HXX