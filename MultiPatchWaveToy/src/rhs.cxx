#include "mpwavetoy.hxx"

namespace MultiPatchWaveToy {

extern "C" void MultiPatchWaveToy_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_MultiPatchWaveToy_RHS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(boundary_condition, "CarpetX")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::pow;

          // Standard finite differencing for first order derivatives
          const CCTK_REAL duda{(-u(-p.DI[0] + p.I) + u(p.DI[0] + p.I)) /
                               (2. * p.DX[0])};

          const CCTK_REAL dudb{(-u(-p.DI[1] + p.I) + u(p.DI[1] + p.I)) /
                               (2. * p.DX[1])};

          const CCTK_REAL dudc{(-u(-p.DI[2] + p.I) + u(p.DI[2] + p.I)) /
                               (2. * p.DX[2])};

          // Standard finite differencing for second order derivatives
          const CCTK_REAL d2udada{
              (-2 * u(p.I) + u(-p.DI[0] + p.I) + u(p.DI[0] + p.I)) /
              pow(p.DX[0], 2)};

          const CCTK_REAL d2udadb{
              (u(-p.DI[0] - p.DI[1] + p.I) - u(p.DI[0] - p.DI[1] + p.I) -
               u(-p.DI[0] + p.DI[1] + p.I) + u(p.DI[0] + p.DI[1] + p.I)) /
              (4. * p.DX[0] * p.DX[1])};

          const CCTK_REAL d2udadc{
              (u(-p.DI[0] - p.DI[2] + p.I) - u(p.DI[0] - p.DI[2] + p.I) -
               u(-p.DI[0] + p.DI[2] + p.I) + u(p.DI[0] + p.DI[2] + p.I)) /
              (4. * p.DX[0] * p.DX[2])};

          const CCTK_REAL d2udbdb{
              (-2 * u(p.I) + u(-p.DI[1] + p.I) + u(p.DI[1] + p.I)) /
              pow(p.DX[1], 2)};

          const CCTK_REAL d2udbdc{
              (u(-p.DI[1] - p.DI[2] + p.I) - u(p.DI[1] - p.DI[2] + p.I) -
               u(-p.DI[1] + p.DI[2] + p.I) + u(p.DI[1] + p.DI[2] + p.I)) /
              (4. * p.DX[1] * p.DX[2])};

          const CCTK_REAL d2udcdc{
              (-2 * u(p.I) + u(-p.DI[2] + p.I) + u(p.DI[2] + p.I)) /
              pow(p.DX[2], 2)};

          // Projecting with Jacobians
          const CCTK_REAL dudxdx{d2udbdb * pow(vJ_db_dx(p.I), 2) +
                                 2 * d2udbdc * vJ_db_dx(p.I) * vJ_dc_dx(p.I) +
                                 d2udcdc * pow(vJ_dc_dx(p.I), 2) +
                                 dudc * (vJ_da_dx(p.I) * vdJ_d2a_dxdz(p.I) +
                                         vJ_db_dx(p.I) * vdJ_d2b_dxdz(p.I)) +
                                 dudc * vJ_dc_dx(p.I) * vdJ_d2c_dxdz(p.I) +
                                 2 * vJ_da_dx(p.I) * vJ_dc_dx(p.I) * d2udadc +
                                 dudb * vJ_da_dx(p.I) * vdJ_d2a_dxdy(p.I) +
                                 dudb * vJ_db_dx(p.I) * vdJ_d2b_dxdy(p.I) +
                                 dudb * vJ_dc_dx(p.I) * vdJ_d2c_dxdy(p.I) +
                                 2 * vJ_da_dx(p.I) * vJ_db_dx(p.I) * d2udadb +
                                 vJ_da_dx(p.I) * duda * vdJ_d2a_dxdx(p.I) +
                                 vJ_db_dx(p.I) * duda * vdJ_d2b_dxdx(p.I) +
                                 vJ_dc_dx(p.I) * duda * vdJ_d2c_dxdx(p.I) +
                                 pow(vJ_da_dx(p.I), 2) * d2udada};

          const CCTK_REAL dudydy{d2udcdc * pow(vJ_dc_dy(p.I), 2) +
                                 dudc * (vJ_da_dy(p.I) * vdJ_d2a_dydz(p.I) +
                                         vJ_db_dy(p.I) * vdJ_d2b_dydz(p.I) +
                                         vJ_dc_dy(p.I) * vdJ_d2c_dydz(p.I)) +
                                 2 * vJ_db_dy(p.I) * vJ_dc_dy(p.I) * d2udbdc +
                                 vJ_da_dy(p.I) * dudb * vdJ_d2a_dydy(p.I) +
                                 vJ_db_dy(p.I) * dudb * vdJ_d2b_dydy(p.I) +
                                 vJ_dc_dy(p.I) * dudb * vdJ_d2c_dydy(p.I) +
                                 pow(vJ_db_dy(p.I), 2) * d2udbdb +
                                 2 * vJ_da_dy(p.I) * vJ_dc_dy(p.I) * d2udadc +
                                 vJ_da_dy(p.I) * duda * vdJ_d2a_dxdy(p.I) +
                                 vJ_db_dy(p.I) * duda * vdJ_d2b_dxdy(p.I) +
                                 vJ_dc_dy(p.I) * duda * vdJ_d2c_dxdy(p.I) +
                                 2 * vJ_da_dy(p.I) * vJ_db_dy(p.I) * d2udadb +
                                 pow(vJ_da_dy(p.I), 2) * d2udada};

          const CCTK_REAL dudzdz{
              pow(vJ_db_dz(p.I), 2) * d2udbdb +
              vJ_db_dz(p.I) *
                  (dudc * vdJ_d2b_dzdz(p.I) + dudb * vdJ_d2b_dydz(p.I) +
                   2 * vJ_dc_dz(p.I) * d2udbdc + duda * vdJ_d2b_dxdz(p.I)) +
              vJ_dc_dz(p.I) *
                  (dudc * vdJ_d2c_dzdz(p.I) + vJ_dc_dz(p.I) * d2udcdc +
                   dudb * vdJ_d2c_dydz(p.I) + duda * vdJ_d2c_dxdz(p.I)) +
              vJ_da_dz(p.I) *
                  (dudc * vdJ_d2a_dzdz(p.I) + dudb * vdJ_d2a_dydz(p.I) +
                   duda * vdJ_d2a_dxdz(p.I) + 2 * vJ_dc_dz(p.I) * d2udadc +
                   2 * vJ_db_dz(p.I) * d2udadb) +
              pow(vJ_da_dz(p.I), 2) * d2udada};

          u_rhs(p.I) = rho(p.I);
          rho_rhs(p.I) = dudxdx + dudydy + dudzdz;
        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace MultiPatchWaveToy