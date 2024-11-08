#include "global_derivatives.hxx"
#include "local_derivatives.hxx"

#include <array>

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

  JacobianComponents(const CCTK_POINTER_TO_CONST cctkGH_,
                     const Loop::PointDesc &p) noexcept {
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

  JacobianDerivativeComponents(const CCTK_POINTER_TO_CONST cctkGH_,
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

// Global Central first derivative, fourth order accurae operator
auto CCTK_HOST CCTK_DEVICE c4o_1(
    const CCTK_POINTER_TO_CONST cctkGH_, const Loop::PointDesc &p,
    const Loop::GF3D2<const CCTK_REAL> &gf) noexcept -> FirstDerivativeResults {
  using namespace Loop;

  // Compute local derivatives
  const auto dgf_da{c4o_1_0_0(p, gf)};
  const auto dgf_db{c4o_0_1_0(p, gf)};
  const auto dgf_dc{c4o_0_0_1(p, gf)};

  // Recover Jacobians
  const auto [da_dx, da_dy, da_dz, db_dx, db_dy, db_dz, dc_dx, dc_dy, dc_dz] =
      JacobianComponents(cctkGH_, p);

  // Projections
  const auto dgf_dx{dgf_db * db_dx + dgf_dc * dc_dx + da_dx * dgf_da};
  const auto dgf_dy{dgf_dc * dc_dy + db_dy * dgf_db + da_dy * dgf_da};
  const auto dgf_dz{dc_dz * dgf_dc + db_dz * dgf_db + da_dz * dgf_da};

  return FirstDerivativeResults{dgf_dx, dgf_dy, dgf_dz};
}

auto CCTK_HOST CCTK_DEVICE c4o_2(const CCTK_POINTER_TO_CONST cctkGH_,
                                 const Loop::PointDesc &p,
                                 const Loop::GF3D2<const CCTK_REAL> &gf)
    noexcept -> SecondDerivativeResults {
  using namespace Loop;

  // Compute local derivatives
  const auto dgf_da{c4o_1_0_0(p, gf)};
  const auto dgf_db{c4o_0_1_0(p, gf)};
  const auto dgf_dc{c4o_0_0_1(p, gf)};

  const auto d2gf_dada{c4o_2_0_0(p, gf)};
  const auto d2gf_dadb{c4o_1_1_0(p, gf)};
  const auto d2gf_dadc{c4o_1_0_1(p, gf)};
  const auto d2gf_dbdb{c4o_0_2_0(p, gf)};
  const auto d2gf_dbdc{c4o_0_1_1(p, gf)};
  const auto d2gf_dcdc{c4o_0_0_2(p, gf)};

  // Recover Jacobians
  const auto [da_dx, da_dy, da_dz, db_dx, db_dy, db_dz, dc_dx, dc_dy, dc_dz] =
      JacobianComponents(cctkGH_, p);

  const auto [d2a_dxdx, d2a_dxdy, d2a_dxdz, d2a_dydy, d2a_dydz, d2a_dzdz,
              d2b_dxdx, d2b_dxdy, d2b_dxdz, d2b_dydy, d2b_dydz, d2b_dzdz,
              d2c_dxdx, d2c_dxdy, d2c_dxdz, d2c_dydy, d2c_dydz, d2c_dzdz] =
      JacobianDerivativeComponents(cctkGH_, p);

  // Projections
  const auto dgf_dxdx{
      dc_dx * (d2gf_dbdc * db_dx + d2gf_dcdc * dc_dx + da_dx * d2gf_dadc) +
      db_dx * (d2gf_dbdb * db_dx + d2gf_dbdc * dc_dx + da_dx * d2gf_dadb) +
      dgf_da * d2a_dxdx + dgf_db * d2b_dxdx + dgf_dc * d2c_dxdx +
      da_dx * (dc_dx * d2gf_dadc + db_dx * d2gf_dadb + da_dx * d2gf_dada)};

  const auto dgf_dxdy{
      dc_dx * (d2gf_dcdc * dc_dy + db_dy * d2gf_dbdc + da_dy * d2gf_dadc) +
      dgf_da * d2a_dxdy + dgf_db * d2b_dxdy + dgf_dc * d2c_dxdy +
      db_dx * (dc_dy * d2gf_dbdc + db_dy * d2gf_dbdb + da_dy * d2gf_dadb) +
      da_dx * (dc_dy * d2gf_dadc + db_dy * d2gf_dadb + da_dy * d2gf_dada)};

  const auto dgf_dxdz{
      dgf_da * d2a_dxdz + dgf_db * d2b_dxdz + dgf_dc * d2c_dxdz +
      dc_dx * (dc_dz * d2gf_dcdc + db_dz * d2gf_dbdc + da_dz * d2gf_dadc) +
      db_dx * (dc_dz * d2gf_dbdc + db_dz * d2gf_dbdb + da_dz * d2gf_dadb) +
      da_dx * (dc_dz * d2gf_dadc + db_dz * d2gf_dadb + da_dz * d2gf_dada)};

  const auto dgf_dydy{
      dgf_db * d2b_dydy + dgf_dc * d2c_dydy + d2a_dydy * dgf_da +
      dc_dy * (d2gf_dcdc * dc_dy + db_dy * d2gf_dbdc + da_dy * d2gf_dadc) +
      db_dy * (dc_dy * d2gf_dbdc + db_dy * d2gf_dbdb + da_dy * d2gf_dadb) +
      da_dy * (dc_dy * d2gf_dadc + db_dy * d2gf_dadb + da_dy * d2gf_dada)};

  const auto dgf_dydz{
      dgf_db * d2b_dydz + dgf_dc * d2c_dydz + d2a_dydz * dgf_da +
      dc_dy * (dc_dz * d2gf_dcdc + db_dz * d2gf_dbdc + da_dz * d2gf_dadc) +
      db_dy * (dc_dz * d2gf_dbdc + db_dz * d2gf_dbdb + da_dz * d2gf_dadb) +
      da_dy * (dc_dz * d2gf_dadc + db_dz * d2gf_dadb + da_dz * d2gf_dada)};

  const auto dgf_dzdz{
      dgf_dc * d2c_dzdz + d2b_dzdz * dgf_db + d2a_dzdz * dgf_da +
      dc_dz * (dc_dz * d2gf_dcdc + db_dz * d2gf_dbdc + da_dz * d2gf_dadc) +
      db_dz * (dc_dz * d2gf_dbdc + db_dz * d2gf_dbdb + da_dz * d2gf_dadb) +
      da_dz * (dc_dz * d2gf_dadc + db_dz * d2gf_dadb + da_dz * d2gf_dada)};

  return SecondDerivativeResults{dgf_dxdx, dgf_dxdy, dgf_dxdz,
                                 dgf_dydy, dgf_dydz, dgf_dzdz};
}

} // namespace MultiPatch::GlobalDerivatives