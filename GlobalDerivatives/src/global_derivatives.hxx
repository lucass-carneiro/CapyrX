#ifndef CAPYRX_GLOBAL_DERIVATIVES_HXX
#define CAPYRX_GLOBAL_DERIVATIVES_HXX

#include <cctk.h>
#include <loop_device.hxx>

#include "local_derivatives.hxx"
#include "jac_recovery.hxx"

namespace MultiPatch::GlobalDerivatives {

struct FirstDerivativeResults {
  CCTK_REAL dx{0.0};
  CCTK_REAL dy{0.0};
  CCTK_REAL dz{0.0};
};

struct FirstDerivativeInputs {
  CCTK_REAL da{0.0};
  CCTK_REAL db{0.0};
  CCTK_REAL dc{0.0};
};

struct SecondDerivativeResults {
  CCTK_REAL dxdx{0.0};
  CCTK_REAL dxdy{0.0};
  CCTK_REAL dxdz{0.0};
  CCTK_REAL dydy{0.0};
  CCTK_REAL dydz{0.0};
  CCTK_REAL dzdz{0.0};
};

// Global Central first derivative, fourth order accurae operator
static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE c4o_1(
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

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
c4o_2(const CCTK_POINTER_TO_CONST cctkGH_, const Loop::PointDesc &p,
      const Loop::GF3D2<const CCTK_REAL> &gf) noexcept
    -> SecondDerivativeResults {
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

#define VERTEX_JACOBIANS_FIRST                                                 \
  vJ_da_dx, vJ_da_dy, vJ_da_dz, vJ_db_dx, vJ_db_dy, vJ_db_dz, vJ_dc_dx,        \
      vJ_dc_dy, vJ_dc_dz

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
project_first(
    const CCTK_REAL &dgf_da, const CCTK_REAL &dgf_db, const CCTK_REAL &dgf_dc,
    const Loop::PointDesc &p, const Loop::GF3D2<const CCTK_REAL> &vJ_da_dx,
    const Loop::GF3D2<const CCTK_REAL> &vJ_da_dy,
    const Loop::GF3D2<const CCTK_REAL> &vJ_da_dz,
    const Loop::GF3D2<const CCTK_REAL> &vJ_db_dx,
    const Loop::GF3D2<const CCTK_REAL> &vJ_db_dy,
    const Loop::GF3D2<const CCTK_REAL> &vJ_db_dz,
    const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dx,
    const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dy,
    const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dz) -> FirstDerivativeResults {

  // Recover Jacobians
  const auto da_dx{vJ_da_dx(p.I)};
  const auto da_dy{vJ_da_dy(p.I)};
  const auto da_dz{vJ_da_dz(p.I)};
  const auto db_dx{vJ_db_dx(p.I)};
  const auto db_dy{vJ_db_dy(p.I)};
  const auto db_dz{vJ_db_dz(p.I)};
  const auto dc_dx{vJ_dc_dx(p.I)};
  const auto dc_dy{vJ_dc_dy(p.I)};
  const auto dc_dz{vJ_dc_dz(p.I)};

  // Projections
  const auto dgf_dx{dgf_db * db_dx + dgf_dc * dc_dx + da_dx * dgf_da};
  const auto dgf_dy{dgf_dc * dc_dy + db_dy * dgf_db + da_dy * dgf_da};
  const auto dgf_dz{dc_dz * dgf_dc + db_dz * dgf_db + da_dz * dgf_da};

  return FirstDerivativeResults{dgf_dx, dgf_dy, dgf_dz};
}

static inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE
project_first(const FirstDerivativeInputs &dgf, const Loop::PointDesc &p,
              const Loop::GF3D2<const CCTK_REAL> &vJ_da_dx,
              const Loop::GF3D2<const CCTK_REAL> &vJ_da_dy,
              const Loop::GF3D2<const CCTK_REAL> &vJ_da_dz,
              const Loop::GF3D2<const CCTK_REAL> &vJ_db_dx,
              const Loop::GF3D2<const CCTK_REAL> &vJ_db_dy,
              const Loop::GF3D2<const CCTK_REAL> &vJ_db_dz,
              const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dx,
              const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dy,
              const Loop::GF3D2<const CCTK_REAL> &vJ_dc_dz)
    -> FirstDerivativeResults {
  return project_first(dgf.da, dgf.db, dgf.dc, p, vJ_da_dx, vJ_da_dy, vJ_da_dz,
                       vJ_db_dx, vJ_db_dy, vJ_db_dz, vJ_dc_dx, vJ_dc_dy,
                       vJ_dc_dz);
}

} // namespace MultiPatch::GlobalDerivatives

#endif // CAPYRX_GLOBAL_DERIVATIVES_HXX