#ifndef MULTIPATCH_INVERSE_CORE_DERIVS_HPP
#define MULTIPATCH_INVERSE_CORE_DERIVS_HPP

namespace MultiPatch {
namespace Cake {

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d_global_to_local_cake_core_da(const PatchTransformations &pt, CCTK_REAL a,
                               CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(b, 2);
  const CCTK_REAL v1 = pow(a, 2);
  const CCTK_REAL v2 = pow(var, 2);
  const CCTK_REAL v3 = -r1;
  const CCTK_REAL v4 = r0 + v3;
  const CCTK_REAL v5 = v0 + v1;
  const CCTK_REAL v6 =
      4 * v2 * v4 * (v3 + r0 * (1 + v5)) + pow(v5, 2) * pow(var, 4);

  return (2 * a * v2 *
          (2 * pow(r0, 2) - 2 * r0 * r1 + v0 * v2 + v1 * v2 + sqrt(v6))) /
         (pow(v4, 2) * sqrt(v6));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d_global_to_local_cake_core_db(const PatchTransformations &pt, CCTK_REAL a,
                               CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(b, 2);
  const CCTK_REAL v1 = pow(a, 2);
  const CCTK_REAL v2 = pow(var, 2);
  const CCTK_REAL v3 = -r1;
  const CCTK_REAL v4 = r0 + v3;
  const CCTK_REAL v5 = v0 + v1;
  const CCTK_REAL v6 =
      4 * v2 * v4 * (v3 + r0 * (1 + v5)) + pow(v5, 2) * pow(var, 4);

  return (2 * b * v2 *
          (2 * pow(r0, 2) - 2 * r0 * r1 + v0 * v2 + v1 * v2 + sqrt(v6))) /
         (pow(v4, 2) * sqrt(v6));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d_global_to_local_cake_core_dvar(const PatchTransformations &pt, CCTK_REAL a,
                                 CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(b, 2);
  const CCTK_REAL v1 = pow(a, 2);
  const CCTK_REAL v2 = -r1;
  const CCTK_REAL v3 = pow(var, 3);
  const CCTK_REAL v4 = pow(r0, 2);
  const CCTK_REAL v5 = r0 + v2;
  const CCTK_REAL v6 = v0 + v1;
  const CCTK_REAL v7 =
      4 * v5 * (v2 + r0 * (1 + v6)) * pow(var, 2) + pow(v6, 2) * pow(var, 4);
  const CCTK_REAL v8 = sqrt(v7);

  return (2 *
          (pow(a, 4) * v3 + pow(b, 4) * v3 + 2 * v0 * v1 * v3 -
           4 * r0 * r1 * var + 2 * pow(r1, 2) * var - 2 * r0 * r1 * v0 * var -
           2 * r0 * r1 * v1 * var + 2 * v4 * var + 2 * v0 * v4 * var +
           2 * v1 * v4 * var + v0 * v8 * var + v1 * v8 * var)) /
         (pow(v5, 2) * sqrt(v7));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d2_global_to_local_cake_core_da2(const PatchTransformations &pt, CCTK_REAL a,
                                 CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(var, 4);
  const CCTK_REAL v1 = pow(b, 2);
  const CCTK_REAL v2 = pow(a, 2);
  const CCTK_REAL v3 = pow(var, 2);
  const CCTK_REAL v4 = -r1;
  const CCTK_REAL v5 = pow(b, 4);
  const CCTK_REAL v6 = pow(a, 4);
  const CCTK_REAL v7 = pow(r1, 2);
  const CCTK_REAL v8 = pow(r0, 2);
  const CCTK_REAL v9 = pow(r0, 3);
  const CCTK_REAL v10 = pow(r0, 4);
  const CCTK_REAL v11 = r0 + v4;
  const CCTK_REAL v12 = v1 + v2;
  const CCTK_REAL v13 = v0 * pow(v12, 2) + 4 * v11 * v3 * (r0 * (1 + v12) + v4);
  const CCTK_REAL v14 = sqrt(v13);

  return (2 * v3 *
          (-8 * r0 * pow(r1, 3) + pow(a, 6) * v0 + pow(b, 6) * v0 + 8 * v10 +
           8 * v1 * v10 - 8 * r0 * r1 * v14 - 4 * r0 * r1 * v1 * v14 -
           4 * r0 * r1 * v14 * v2 - 8 * r0 * r1 * v1 * v3 -
           24 * r0 * r1 * v2 * v3 - 12 * r0 * r1 * v1 * v2 * v3 +
           2 * v1 * v14 * v2 * v3 + 3 * v0 * v2 * v5 - 6 * r0 * r1 * v3 * v5 +
           v14 * v3 * v5 + 3 * v0 * v1 * v6 - 6 * r0 * r1 * v3 * v6 +
           v14 * v3 * v6 + 4 * v14 * v7 + 4 * v1 * v3 * v7 + 12 * v2 * v3 * v7 +
           4 * v14 * v8 + 4 * v1 * v14 * v8 + 4 * v14 * v2 * v8 +
           4 * v1 * v3 * v8 + 12 * v2 * v3 * v8 + 12 * v1 * v2 * v3 * v8 +
           6 * v3 * v5 * v8 + 6 * v3 * v6 * v8 + 24 * v7 * v8 +
           8 * v1 * v7 * v8 - 24 * r1 * v9 - 16 * r1 * v1 * v9)) /
         (pow(v11, 2) * sqrt(v13) *
          (-8 * r0 * r1 - 4 * r0 * r1 * v1 - 4 * r0 * r1 * v2 +
           2 * v1 * v2 * v3 + v3 * v5 + v3 * v6 + 4 * v7 + 4 * v8 +
           4 * v1 * v8 + 4 * v2 * v8));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d2_global_to_local_cake_core_dadb(const PatchTransformations &pt, CCTK_REAL a,
                                  CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(var, 4);
  const CCTK_REAL v1 = pow(b, 2);
  const CCTK_REAL v2 = pow(a, 2);
  const CCTK_REAL v3 = pow(var, 2);
  const CCTK_REAL v4 = -r1;

  return (16 * a * b * v0 * (-pow(r0, 2) + v3)) /
         pow(v0 * pow(v1 + v2, 2) +
                 4 * v3 * (r0 + v4) * (r0 * (1 + v1 + v2) + v4),
             1.5);
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d2_global_to_local_cake_core_dadvar(const PatchTransformations &pt, CCTK_REAL a,
                                    CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(var, 4);
  const CCTK_REAL v1 = pow(b, 2);
  const CCTK_REAL v2 = pow(a, 2);
  const CCTK_REAL v3 = pow(var, 2);
  const CCTK_REAL v4 = -r1;
  const CCTK_REAL v5 = pow(b, 4);
  const CCTK_REAL v6 = pow(a, 3);
  const CCTK_REAL v7 = pow(a, 5);
  const CCTK_REAL v8 = pow(r1, 2);
  const CCTK_REAL v9 = pow(r0, 2);
  const CCTK_REAL v10 = pow(r0, 3);
  const CCTK_REAL v11 = pow(r0, 4);
  const CCTK_REAL v12 = r0 + v4;
  const CCTK_REAL v13 = v1 + v2;
  const CCTK_REAL v14 = v0 * pow(v13, 2) + 4 * v12 * v3 * (r0 * (1 + v13) + v4);
  const CCTK_REAL v15 = sqrt(v14);

  return (4 *
          (-4 * a * r0 * pow(r1, 3) + pow(a, 7) * v0 + a * pow(b, 6) * v0 -
           12 * a * r1 * v10 - 8 * a * r1 * v1 * v10 + 4 * a * v11 +
           4 * a * v1 * v11 - 8 * a * r0 * r1 * v15 -
           4 * a * r0 * r1 * v1 * v15 - 12 * a * r0 * r1 * v1 * v3 -
           6 * a * r0 * r1 * v3 * v5 + a * v15 * v3 * v5 - 8 * r1 * v10 * v6 +
           4 * v11 * v6 - 4 * r0 * r1 * v15 * v6 - 12 * r0 * r1 * v3 * v6 -
           12 * r0 * r1 * v1 * v3 * v6 + 2 * v1 * v15 * v3 * v6 +
           3 * v0 * v5 * v6 + 3 * v0 * v1 * v7 - 6 * r0 * r1 * v3 * v7 +
           v15 * v3 * v7 + 4 * a * v15 * v8 + 6 * a * v1 * v3 * v8 +
           6 * v3 * v6 * v8 + 4 * a * v15 * v9 + 4 * a * v1 * v15 * v9 +
           6 * a * v1 * v3 * v9 + 6 * a * v3 * v5 * v9 + 4 * v15 * v6 * v9 +
           6 * v3 * v6 * v9 + 12 * v1 * v3 * v6 * v9 + 6 * v3 * v7 * v9 +
           12 * a * v8 * v9 + 4 * a * v1 * v8 * v9 + 4 * v6 * v8 * v9) *
          var) /
         (pow(v12, 2) * sqrt(v14) *
          (-8 * r0 * r1 - 4 * r0 * r1 * v1 - 4 * r0 * r1 * v2 + pow(a, 4) * v3 +
           2 * v1 * v2 * v3 + v3 * v5 + 4 * v8 + 4 * v9 + 4 * v1 * v9 +
           4 * v2 * v9));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d2_global_to_local_cake_core_db2(const PatchTransformations &pt, CCTK_REAL a,
                                 CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(var, 4);
  const CCTK_REAL v1 = pow(b, 2);
  const CCTK_REAL v2 = pow(a, 2);
  const CCTK_REAL v3 = pow(var, 2);
  const CCTK_REAL v4 = -r1;
  const CCTK_REAL v5 = pow(b, 4);
  const CCTK_REAL v6 = pow(a, 4);
  const CCTK_REAL v7 = pow(r1, 2);
  const CCTK_REAL v8 = pow(r0, 2);
  const CCTK_REAL v9 = pow(r0, 3);
  const CCTK_REAL v10 = pow(r0, 4);
  const CCTK_REAL v11 = r0 + v4;
  const CCTK_REAL v12 = v1 + v2;
  const CCTK_REAL v13 = v0 * pow(v12, 2) + 4 * v11 * v3 * (r0 * (1 + v12) + v4);
  const CCTK_REAL v14 = sqrt(v13);

  return (2 * v3 *
          (-8 * r0 * pow(r1, 3) + pow(a, 6) * v0 + pow(b, 6) * v0 + 8 * v10 -
           8 * r0 * r1 * v14 - 4 * r0 * r1 * v1 * v14 + 8 * v10 * v2 -
           4 * r0 * r1 * v14 * v2 - 24 * r0 * r1 * v1 * v3 -
           8 * r0 * r1 * v2 * v3 - 12 * r0 * r1 * v1 * v2 * v3 +
           2 * v1 * v14 * v2 * v3 + 3 * v0 * v2 * v5 - 6 * r0 * r1 * v3 * v5 +
           v14 * v3 * v5 + 3 * v0 * v1 * v6 - 6 * r0 * r1 * v3 * v6 +
           v14 * v3 * v6 + 4 * v14 * v7 + 12 * v1 * v3 * v7 + 4 * v2 * v3 * v7 +
           4 * v14 * v8 + 4 * v1 * v14 * v8 + 4 * v14 * v2 * v8 +
           12 * v1 * v3 * v8 + 4 * v2 * v3 * v8 + 12 * v1 * v2 * v3 * v8 +
           6 * v3 * v5 * v8 + 6 * v3 * v6 * v8 + 24 * v7 * v8 +
           8 * v2 * v7 * v8 - 24 * r1 * v9 - 16 * r1 * v2 * v9)) /
         (pow(v11, 2) * sqrt(v13) *
          (-8 * r0 * r1 - 4 * r0 * r1 * v1 - 4 * r0 * r1 * v2 +
           2 * v1 * v2 * v3 + v3 * v5 + v3 * v6 + 4 * v7 + 4 * v8 +
           4 * v1 * v8 + 4 * v2 * v8));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d2_global_to_local_cake_core_dbdvar(const PatchTransformations &pt, CCTK_REAL a,
                                    CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(var, 4);
  const CCTK_REAL v1 = pow(b, 2);
  const CCTK_REAL v2 = pow(a, 2);
  const CCTK_REAL v3 = pow(var, 2);
  const CCTK_REAL v4 = -r1;
  const CCTK_REAL v5 = pow(b, 5);
  const CCTK_REAL v6 = pow(b, 3);
  const CCTK_REAL v7 = pow(a, 4);
  const CCTK_REAL v8 = pow(r1, 2);
  const CCTK_REAL v9 = pow(r0, 2);
  const CCTK_REAL v10 = pow(r0, 3);
  const CCTK_REAL v11 = pow(r0, 4);
  const CCTK_REAL v12 = r0 + v4;
  const CCTK_REAL v13 = v1 + v2;
  const CCTK_REAL v14 = v0 * pow(v13, 2) + 4 * v12 * v3 * (r0 * (1 + v13) + v4);
  const CCTK_REAL v15 = sqrt(v14);

  return (4 *
          (-4 * b * r0 * pow(r1, 3) + pow(a, 6) * b * v0 + pow(b, 7) * v0 -
           12 * b * r1 * v10 + 4 * b * v11 - 8 * b * r0 * r1 * v15 -
           8 * b * r1 * v10 * v2 + 4 * b * v11 * v2 -
           4 * b * r0 * r1 * v15 * v2 - 12 * b * r0 * r1 * v2 * v3 +
           3 * v0 * v2 * v5 - 6 * r0 * r1 * v3 * v5 + v15 * v3 * v5 -
           8 * r1 * v10 * v6 + 4 * v11 * v6 - 4 * r0 * r1 * v15 * v6 -
           12 * r0 * r1 * v3 * v6 - 12 * r0 * r1 * v2 * v3 * v6 +
           2 * v15 * v2 * v3 * v6 - 6 * b * r0 * r1 * v3 * v7 +
           b * v15 * v3 * v7 + 3 * v0 * v6 * v7 + 4 * b * v15 * v8 +
           6 * b * v2 * v3 * v8 + 6 * v3 * v6 * v8 + 4 * b * v15 * v9 +
           4 * b * v15 * v2 * v9 + 6 * b * v2 * v3 * v9 + 6 * v3 * v5 * v9 +
           4 * v15 * v6 * v9 + 6 * v3 * v6 * v9 + 12 * v2 * v3 * v6 * v9 +
           6 * b * v3 * v7 * v9 + 12 * b * v8 * v9 + 4 * b * v2 * v8 * v9 +
           4 * v6 * v8 * v9) *
          var) /
         (pow(v12, 2) * sqrt(v14) *
          (-8 * r0 * r1 - 4 * r0 * r1 * v1 - 4 * r0 * r1 * v2 + pow(b, 4) * v3 +
           2 * v1 * v2 * v3 + v3 * v7 + 4 * v8 + 4 * v9 + 4 * v1 * v9 +
           4 * v2 * v9));
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
d2_global_to_local_cake_core_dvar2(const PatchTransformations &pt, CCTK_REAL a,
                                   CCTK_REAL b, CCTK_REAL var) {

  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const CCTK_REAL v0 = pow(var, 4);
  const CCTK_REAL v1 = pow(b, 2);
  const CCTK_REAL v2 = pow(a, 2);
  const CCTK_REAL v3 = pow(var, 2);
  const CCTK_REAL v4 = -r1;
  const CCTK_REAL v5 = pow(b, 6);
  const CCTK_REAL v6 = pow(b, 4);
  const CCTK_REAL v7 = pow(a, 4);
  const CCTK_REAL v8 = pow(a, 6);
  const CCTK_REAL v9 = pow(r1, 2);
  const CCTK_REAL v10 = pow(r0, 2);
  const CCTK_REAL v11 = r0 + v4;
  const CCTK_REAL v12 = v1 + v2;
  const CCTK_REAL v13 = v0 * pow(v12, 2) + 4 * v11 * v3 * (r0 * (1 + v12) + v4);
  const CCTK_REAL v14 = sqrt(v13);

  return (2 *
          (pow(a, 8) * v0 + pow(b, 8) * v0 - 8 * r0 * r1 * v1 * v14 +
           4 * v1 * v10 * v14 - 8 * r0 * r1 * v14 * v2 -
           8 * r0 * r1 * v1 * v14 * v2 + 4 * v10 * v14 * v2 +
           8 * v1 * v10 * v14 * v2 - 24 * r0 * r1 * v1 * v2 * v3 +
           12 * v1 * v10 * v2 * v3 + 4 * v0 * v2 * v5 - 6 * r0 * r1 * v3 * v5 +
           6 * v10 * v3 * v5 + v14 * v3 * v5 - 4 * r0 * r1 * v14 * v6 +
           4 * v10 * v14 * v6 - 12 * r0 * r1 * v3 * v6 + 6 * v10 * v3 * v6 -
           18 * r0 * r1 * v2 * v3 * v6 + 18 * v10 * v2 * v3 * v6 +
           3 * v14 * v2 * v3 * v6 - 4 * r0 * r1 * v14 * v7 +
           4 * v10 * v14 * v7 - 12 * r0 * r1 * v3 * v7 -
           18 * r0 * r1 * v1 * v3 * v7 + 6 * v10 * v3 * v7 +
           18 * v1 * v10 * v3 * v7 + 3 * v1 * v14 * v3 * v7 + 6 * v0 * v6 * v7 +
           4 * v0 * v1 * v8 - 6 * r0 * r1 * v3 * v8 + 6 * v10 * v3 * v8 +
           v14 * v3 * v8 + 4 * v1 * v14 * v9 + 4 * v14 * v2 * v9 +
           12 * v1 * v2 * v3 * v9 + 6 * v3 * v6 * v9 + 6 * v3 * v7 * v9)) /
         (pow(v11, 2) * sqrt(v13) *
          (-8 * r0 * r1 - 4 * r0 * r1 * v1 + 4 * v10 + 4 * v1 * v10 -
           4 * r0 * r1 * v2 + 4 * v10 * v2 + 2 * v1 * v2 * v3 + v3 * v6 +
           v3 * v7 + 4 * v9));
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_INVERSE_CORE_DERIVS_HPP