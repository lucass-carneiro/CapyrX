#ifndef MULTIPATCH_CAKE_JACOBIANS_HXX
#define MULTIPATCH_CAKE_JACOBIANS_HXX

#include "multipatch.hxx"
#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

CCTK_DEVICE CCTK_HOST inline std_tuple<jac_t, djac_t>
cake_jacs(const PatchTransformations &pt, int patch, const svec &global_vars) {
  using std::pow;
  using std::sqrt;

  jac_t J{};
  djac_t dJ{};

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  switch (patch) {

  case static_cast<int>(patch_piece::cartesian):
    J(0)(0) = 1 / r0;
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)(0) = 0;
    J(1)(1) = 1 / r0;
    J(1)(2) = 0;
    J(2)(0) = 0;
    J(2)(1) = 0;
    J(2)(2) = 1 / r0;

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = 0;
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = 0;
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) = 0;
    dJ(2)(0, 1) = 0;
    dJ(2)(0, 2) = 0;
    dJ(2)(1, 0) = 0;
    dJ(2)(1, 1) = 0;
    dJ(2)(1, 2) = 0;
    dJ(2)(2, 0) = 0;
    dJ(2)(2, 1) = 0;
    dJ(2)(2, 2) = 0;
    break;

  case static_cast<int>(patch_piece::plus_x):
    J(0)(0) = -(z / pow(x, 2));
    J(0)(1) = 0;
    J(0)(2) = 1 / x;
    J(1)(0) = -(y / pow(x, 2));
    J(1)(1) = 1 / x;
    J(1)(2) = 0;
    J(2)
    (0) = (4 * x) /
          sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
               4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)
    (1) =
        (2 * y *
         (1 +
          (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);
    J(2)
    (2) =
        (2 * z *
         (1 +
          (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);

    dJ(0)(0, 0) = (2 * z) / pow(x, 3);
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = -pow(x, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = -pow(x, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = (2 * y) / pow(x, 3);
    dJ(1)(0, 1) = -pow(x, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = -pow(x, -2);
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        (4 * (pow(y, 2) + pow(z, 2)) *
         (4 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(0, 1) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(0, 2) =
        (-8 * x * z * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 0) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 1) =
        (2 -
         (4 * pow(y, 2) * pow(2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + 3 * pow(y, 2) + pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    dJ(2)(1, 2) =
        (16 * (-pow(r0, 2) + pow(x, 2)) * y * z) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 0) =
        (-8 * x * z * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 1) =
        (16 * (-pow(r0, 2) + pow(x, 2)) * y * z) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 2) =
        (2 -
         (4 * pow(z, 2) * pow(2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + pow(y, 2) + 3 * pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::minus_x):
    J(0)(0) = z / pow(x, 2);
    J(0)(1) = 0;
    J(0)(2) = -(1 / x);
    J(1)(0) = -(y / pow(x, 2));
    J(1)(1) = 1 / x;
    J(1)(2) = 0;
    J(2)
    (0) = (4 * x) /
          sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
               4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)
    (1) =
        (2 * y *
         (1 +
          (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);
    J(2)
    (2) =
        (2 * z *
         (1 +
          (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);

    dJ(0)(0, 0) = (-2 * z) / pow(x, 3);
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = pow(x, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = pow(x, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = (2 * y) / pow(x, 3);
    dJ(1)(0, 1) = -pow(x, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = -pow(x, -2);
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        (4 * (pow(y, 2) + pow(z, 2)) *
         (4 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(0, 1) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(0, 2) =
        (-8 * x * z * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 0) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 1) =
        (2 -
         (4 * pow(y, 2) * pow(2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + 3 * pow(y, 2) + pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    dJ(2)(1, 2) =
        (16 * (-pow(r0, 2) + pow(x, 2)) * y * z) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 0) =
        (-8 * x * z * (2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 1) =
        (16 * (-pow(r0, 2) + pow(x, 2)) * y * z) /
        pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 2) =
        (2 -
         (4 * pow(z, 2) * pow(2 * r0 * (r0 - r1) + pow(y, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + pow(y, 2) + 3 * pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::plus_y):
    J(0)(0) = 0;
    J(0)(1) = -(z / pow(y, 2));
    J(0)(2) = 1 / y;
    J(1)(0) = -(1 / y);
    J(1)(1) = x / pow(y, 2);
    J(1)(2) = 0;
    J(2)
    (0) =
        (2 * x *
         (1 +
          (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);
    J(2)
    (1) = (4 * y) /
          sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
               4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)));
    J(2)
    (2) =
        (2 * z *
         (1 +
          (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = (2 * z) / pow(y, 3);
    dJ(0)(1, 2) = -pow(y, -2);
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = -pow(y, -2);
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = pow(y, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = pow(y, -2);
    dJ(1)(1, 1) = (-2 * x) / pow(y, 3);
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        (2 -
         (4 * pow(x, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + 3 * pow(x, 2) + pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    dJ(2)(0, 1) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(0, 2) =
        (16 * x * (-pow(r0, 2) + pow(y, 2)) * z) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 0) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 1) =
        (4 * (pow(x, 2) + pow(z, 2)) *
         (4 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 2) =
        (-8 * y * z * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 0) =
        (16 * x * (-pow(r0, 2) + pow(y, 2)) * z) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 1) =
        (-8 * y * z * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 2) =
        (2 -
         (4 * pow(z, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + pow(x, 2) + 3 * pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::minus_y):
    J(0)(0) = 0;
    J(0)(1) = z / pow(y, 2);
    J(0)(2) = -(1 / y);
    J(1)(0) = -(1 / y);
    J(1)(1) = x / pow(y, 2);
    J(1)(2) = 0;
    J(2)
    (0) =
        (2 * x *
         (1 +
          (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);
    J(2)
    (1) = (4 * y) /
          sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
               4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
               4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)));
    J(2)
    (2) =
        (2 * z *
         (1 +
          (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2)) /
              sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                   4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                   4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2))))) /
        pow(r0 - r1, 2);

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = (-2 * z) / pow(y, 3);
    dJ(0)(1, 2) = pow(y, -2);
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = pow(y, -2);
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = pow(y, -2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = pow(y, -2);
    dJ(1)(1, 1) = (-2 * x) / pow(y, 3);
    dJ(1)(1, 2) = 0;
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = 0;
    dJ(1)(2, 2) = 0;
    dJ(2)(0, 0) =
        (2 -
         (4 * pow(x, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + 3 * pow(x, 2) + pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    dJ(2)(0, 1) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(0, 2) =
        (16 * x * (-pow(r0, 2) + pow(y, 2)) * z) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 0) =
        (-8 * x * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 1) =
        (4 * (pow(x, 2) + pow(z, 2)) *
         (4 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(1, 2) =
        (-8 * y * z * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 0) =
        (16 * x * (-pow(r0, 2) + pow(y, 2)) * z) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 1) =
        (-8 * y * z * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2))) /
        pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
            1.5);
    dJ(2)(2, 2) =
        (2 -
         (4 * pow(z, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(z, 2), 2)) /
             pow(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                     4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                     4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + pow(x, 2) + 3 * pow(z, 2))) /
             sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
                  4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
                  4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::plus_z):
    J(0)(0) = -(1 / z);
    J(0)(1) = 0;
    J(0)(2) = x / pow(z, 2);
    J(1)(0) = 0;
    J(1)(1) = 1 / z;
    J(1)(2) = -(y / pow(z, 2));
    J(2)
    (0) = (2 * x *
           (1 + (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) /
                    sqrt((pow(x, 2) + pow(y, 2)) *
                             (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                         4 * pow(r0 - r1, 2) * pow(z, 2)))) /
          pow(r0 - r1, 2);
    J(2)
    (1) = (2 * y *
           (1 + (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) /
                    sqrt((pow(x, 2) + pow(y, 2)) *
                             (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                         4 * pow(r0 - r1, 2) * pow(z, 2)))) /
          pow(r0 - r1, 2);
    J(2)
    (2) = (4 * z) / sqrt((pow(x, 2) + pow(y, 2)) *
                             (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                         4 * pow(r0 - r1, 2) * pow(z, 2));

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = pow(z, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = pow(z, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = (-2 * x) / pow(z, 3);
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = 0;
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = 0;
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = -pow(z, -2);
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = -pow(z, -2);
    dJ(1)(2, 2) = (2 * y) / pow(z, 3);
    dJ(2)(0, 0) =
        (2 -
         (4 * pow(x, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2), 2)) /
             pow((pow(x, 2) + pow(y, 2)) *
                         (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                     4 * pow(r0 - r1, 2) * pow(z, 2),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + 3 * pow(x, 2) + pow(y, 2))) /
             sqrt((pow(x, 2) + pow(y, 2)) *
                      (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                  4 * pow(r0 - r1, 2) * pow(z, 2))) /
        pow(r0 - r1, 2);
    dJ(2)(0, 1) = (16 * x * y * (-pow(r0, 2) + pow(z, 2))) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(0, 2) = (-8 * x * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(1, 0) = (16 * x * y * (-pow(r0, 2) + pow(z, 2))) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(1, 1) =
        (2 -
         (4 * pow(y, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2), 2)) /
             pow((pow(x, 2) + pow(y, 2)) *
                         (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                     4 * pow(r0 - r1, 2) * pow(z, 2),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + pow(x, 2) + 3 * pow(y, 2))) /
             sqrt((pow(x, 2) + pow(y, 2)) *
                      (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                  4 * pow(r0 - r1, 2) * pow(z, 2))) /
        pow(r0 - r1, 2);
    dJ(2)(1, 2) = (-8 * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(2, 0) = (-8 * x * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(2, 1) = (-8 * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(2, 2) = (4 * (pow(x, 2) + pow(y, 2)) *
                   (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2))) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    break;

  case static_cast<int>(patch_piece::minus_z):
    J(0)(0) = -(1 / z);
    J(0)(1) = 0;
    J(0)(2) = x / pow(z, 2);
    J(1)(0) = 0;
    J(1)(1) = -(1 / z);
    J(1)(2) = y / pow(z, 2);
    J(2)
    (0) = (2 * x *
           (1 + (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) /
                    sqrt((pow(x, 2) + pow(y, 2)) *
                             (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                         4 * pow(r0 - r1, 2) * pow(z, 2)))) /
          pow(r0 - r1, 2);
    J(2)
    (1) = (2 * y *
           (1 + (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) /
                    sqrt((pow(x, 2) + pow(y, 2)) *
                             (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                         4 * pow(r0 - r1, 2) * pow(z, 2)))) /
          pow(r0 - r1, 2);
    J(2)
    (2) = (4 * z) / sqrt((pow(x, 2) + pow(y, 2)) *
                             (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                         4 * pow(r0 - r1, 2) * pow(z, 2));

    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = pow(z, -2);
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = pow(z, -2);
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = (-2 * x) / pow(z, 3);
    dJ(1)(0, 0) = 0;
    dJ(1)(0, 1) = 0;
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = 0;
    dJ(1)(1, 1) = 0;
    dJ(1)(1, 2) = pow(z, -2);
    dJ(1)(2, 0) = 0;
    dJ(1)(2, 1) = pow(z, -2);
    dJ(1)(2, 2) = (-2 * y) / pow(z, 3);
    dJ(2)(0, 0) =
        (2 -
         (4 * pow(x, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2), 2)) /
             pow((pow(x, 2) + pow(y, 2)) *
                         (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                     4 * pow(r0 - r1, 2) * pow(z, 2),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + 3 * pow(x, 2) + pow(y, 2))) /
             sqrt((pow(x, 2) + pow(y, 2)) *
                      (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                  4 * pow(r0 - r1, 2) * pow(z, 2))) /
        pow(r0 - r1, 2);
    dJ(2)(0, 1) = (16 * x * y * (-pow(r0, 2) + pow(z, 2))) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(0, 2) = (-8 * x * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(1, 0) = (16 * x * y * (-pow(r0, 2) + pow(z, 2))) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(1, 1) =
        (2 -
         (4 * pow(y, 2) * pow(2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2), 2)) /
             pow((pow(x, 2) + pow(y, 2)) *
                         (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                     4 * pow(r0 - r1, 2) * pow(z, 2),
                 1.5) +
         (2 * (2 * r0 * (r0 - r1) + pow(x, 2) + 3 * pow(y, 2))) /
             sqrt((pow(x, 2) + pow(y, 2)) *
                      (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                  4 * pow(r0 - r1, 2) * pow(z, 2))) /
        pow(r0 - r1, 2);
    dJ(2)(1, 2) = (-8 * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(2, 0) = (-8 * x * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(2, 1) = (-8 * y * (2 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) * z) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    dJ(2)(2, 2) = (4 * (pow(x, 2) + pow(y, 2)) *
                   (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2))) /
                  pow((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2),
                      1.5);
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No jacobians available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_JACOBIANS_HXX