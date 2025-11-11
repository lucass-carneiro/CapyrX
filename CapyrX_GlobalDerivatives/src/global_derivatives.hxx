#ifndef CAPYRX_GLOBAL_DERIVATIVES_HXX
#define CAPYRX_GLOBAL_DERIVATIVES_HXX

#include <cctk.h>

namespace CapyrX::MultiPatch::GlobalDerivatives {

#define VERTEX_JACOBIANS(point_descriptor)                                     \
  vJ_da_dx(point_descriptor.I), vJ_da_dy(point_descriptor.I),                  \
      vJ_da_dz(point_descriptor.I), vJ_db_dx(point_descriptor.I),              \
      vJ_db_dy(point_descriptor.I), vJ_db_dz(point_descriptor.I),              \
      vJ_dc_dx(point_descriptor.I), vJ_dc_dy(point_descriptor.I),              \
      vJ_dc_dz(point_descriptor.I)

#define VERTEX_DJACOBIANS(point_descriptor)                                    \
  vdJ_d2a_dxdx(point_descriptor.I), vdJ_d2a_dxdy(point_descriptor.I),          \
      vdJ_d2a_dxdz(point_descriptor.I), vdJ_d2a_dydy(point_descriptor.I),      \
      vdJ_d2a_dydz(point_descriptor.I), vdJ_d2a_dzdz(point_descriptor.I),      \
      vdJ_d2b_dxdx(point_descriptor.I), vdJ_d2b_dxdy(point_descriptor.I),      \
      vdJ_d2b_dxdz(point_descriptor.I), vdJ_d2b_dydy(point_descriptor.I),      \
      vdJ_d2b_dydz(point_descriptor.I), vdJ_d2b_dzdz(point_descriptor.I),      \
      vdJ_d2c_dxdx(point_descriptor.I), vdJ_d2c_dxdy(point_descriptor.I),      \
      vdJ_d2c_dxdz(point_descriptor.I), vdJ_d2c_dydy(point_descriptor.I),      \
      vdJ_d2c_dydz(point_descriptor.I), vdJ_d2c_dzdz(point_descriptor.I)

struct LocalFirstDerivatives {
  CCTK_REAL da{0.0};
  CCTK_REAL db{0.0};
  CCTK_REAL dc{0.0};
};

struct LocalSecondDerivatives {
  CCTK_REAL dada{0.0};
  CCTK_REAL dadb{0.0};
  CCTK_REAL dadc{0.0};
  CCTK_REAL dbdb{0.0};
  CCTK_REAL dbdc{0.0};
  CCTK_REAL dcdc{0.0};
};

struct GlobalFirstDerivatives {
  CCTK_REAL dx{0.0};
  CCTK_REAL dy{0.0};
  CCTK_REAL dz{0.0};
};

struct GlobalSecondDerivatives {
  CCTK_REAL dxdx{0.0};
  CCTK_REAL dxdy{0.0};
  CCTK_REAL dxdz{0.0};
  CCTK_REAL dydy{0.0};
  CCTK_REAL dydz{0.0};
  CCTK_REAL dzdz{0.0};
};

struct Jacobians {
  CCTK_REAL da_dx{0.0};
  CCTK_REAL da_dy{0.0};
  CCTK_REAL da_dz{0.0};
  CCTK_REAL db_dx{0.0};
  CCTK_REAL db_dy{0.0};
  CCTK_REAL db_dz{0.0};
  CCTK_REAL dc_dx{0.0};
  CCTK_REAL dc_dy{0.0};
  CCTK_REAL dc_dz{0.0};
};

struct JacobianDerivatives {
  CCTK_REAL d2a_dxdx{0.0};
  CCTK_REAL d2a_dxdy{0.0};
  CCTK_REAL d2a_dxdz{0.0};
  CCTK_REAL d2a_dydy{0.0};
  CCTK_REAL d2a_dydz{0.0};
  CCTK_REAL d2a_dzdz{0.0};
  CCTK_REAL d2b_dxdx{0.0};
  CCTK_REAL d2b_dxdy{0.0};
  CCTK_REAL d2b_dxdz{0.0};
  CCTK_REAL d2b_dydy{0.0};
  CCTK_REAL d2b_dydz{0.0};
  CCTK_REAL d2b_dzdz{0.0};
  CCTK_REAL d2c_dxdx{0.0};
  CCTK_REAL d2c_dxdy{0.0};
  CCTK_REAL d2c_dxdz{0.0};
  CCTK_REAL d2c_dydy{0.0};
  CCTK_REAL d2c_dydz{0.0};
  CCTK_REAL d2c_dzdz{0.0};
};

static inline auto CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE
project_first(const LocalFirstDerivatives &lfd, const Jacobians &jac)
    -> GlobalFirstDerivatives {

  // Recover local derivatives
  const auto &[dgf_da, dgf_db, dgf_dc] = lfd;

  // Recover Jacobians
  const auto &[da_dx, da_dy, da_dz, db_dx, db_dy, db_dz, dc_dx, dc_dy, dc_dz] =
      jac;

  // Project
  return GlobalFirstDerivatives{
      dgf_db * db_dx + dgf_dc * dc_dx + da_dx * dgf_da,
      dgf_dc * dc_dy + db_dy * dgf_db + da_dy * dgf_da,
      dc_dz * dgf_dc + db_dz * dgf_db + da_dz * dgf_da,
  };
}

static inline auto CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE
project_second(const LocalFirstDerivatives &lfd,
               const LocalSecondDerivatives &lsd, const Jacobians &jac,
               const JacobianDerivatives &djac) -> GlobalSecondDerivatives {

  // Recover local derivatives
  const auto &[dgf_da, dgf_db, dgf_dc] = lfd;

  const auto &[d2gf_dada, d2gf_dadb, d2gf_dadc, d2gf_dbdb, d2gf_dbdc,
               d2gf_dcdc] = lsd;

  // Recover Jacobians
  const auto &[da_dx, da_dy, da_dz, db_dx, db_dy, db_dz, dc_dx, dc_dy, dc_dz] =
      jac;

  const auto &[d2a_dxdx, d2a_dxdy, d2a_dxdz, d2a_dydy, d2a_dydz, d2a_dzdz,
               d2b_dxdx, d2b_dxdy, d2b_dxdz, d2b_dydy, d2b_dydz, d2b_dzdz,
               d2c_dxdx, d2c_dxdy, d2c_dxdz, d2c_dydy, d2c_dydz, d2c_dzdz] =
      djac;

  // Project
  return GlobalSecondDerivatives{
      dc_dx * (d2gf_dbdc * db_dx + d2gf_dcdc * dc_dx + da_dx * d2gf_dadc) +
          db_dx * (d2gf_dbdb * db_dx + d2gf_dbdc * dc_dx + da_dx * d2gf_dadb) +
          dgf_da * d2a_dxdx + dgf_db * d2b_dxdx + dgf_dc * d2c_dxdx +
          da_dx * (dc_dx * d2gf_dadc + db_dx * d2gf_dadb + da_dx * d2gf_dada),

      dc_dx * (d2gf_dcdc * dc_dy + db_dy * d2gf_dbdc + da_dy * d2gf_dadc) +
          dgf_da * d2a_dxdy + dgf_db * d2b_dxdy + dgf_dc * d2c_dxdy +
          db_dx * (dc_dy * d2gf_dbdc + db_dy * d2gf_dbdb + da_dy * d2gf_dadb) +
          da_dx * (dc_dy * d2gf_dadc + db_dy * d2gf_dadb + da_dy * d2gf_dada),

      dgf_da * d2a_dxdz + dgf_db * d2b_dxdz + dgf_dc * d2c_dxdz +
          dc_dx * (dc_dz * d2gf_dcdc + db_dz * d2gf_dbdc + da_dz * d2gf_dadc) +
          db_dx * (dc_dz * d2gf_dbdc + db_dz * d2gf_dbdb + da_dz * d2gf_dadb) +
          da_dx * (dc_dz * d2gf_dadc + db_dz * d2gf_dadb + da_dz * d2gf_dada),

      dgf_db * d2b_dydy + dgf_dc * d2c_dydy + d2a_dydy * dgf_da +
          dc_dy * (d2gf_dcdc * dc_dy + db_dy * d2gf_dbdc + da_dy * d2gf_dadc) +
          db_dy * (dc_dy * d2gf_dbdc + db_dy * d2gf_dbdb + da_dy * d2gf_dadb) +
          da_dy * (dc_dy * d2gf_dadc + db_dy * d2gf_dadb + da_dy * d2gf_dada),

      dgf_db * d2b_dydz + dgf_dc * d2c_dydz + d2a_dydz * dgf_da +
          dc_dy * (dc_dz * d2gf_dcdc + db_dz * d2gf_dbdc + da_dz * d2gf_dadc) +
          db_dy * (dc_dz * d2gf_dbdc + db_dz * d2gf_dbdb + da_dz * d2gf_dadb) +
          da_dy * (dc_dz * d2gf_dadc + db_dz * d2gf_dadb + da_dz * d2gf_dada),

      dgf_dc * d2c_dzdz + d2b_dzdz * dgf_db + d2a_dzdz * dgf_da +
          dc_dz * (dc_dz * d2gf_dcdc + db_dz * d2gf_dbdc + da_dz * d2gf_dadc) +
          db_dz * (dc_dz * d2gf_dbdc + db_dz * d2gf_dbdb + da_dz * d2gf_dadb) +
          da_dz * (dc_dz * d2gf_dadc + db_dz * d2gf_dadb + da_dz * d2gf_dada)};
}

} // namespace CapyrX::MultiPatch::GlobalDerivatives

#endif // CAPYRX_GLOBAL_DERIVATIVES_HXX