#include "cubed_sphere.hxx"

#include <cmath>
#include <random>

namespace CapyrX::MultiPatch::CubedSphere {

enum class PatchPiece : int {
  cartesian = 0,

  plus_x = 1,
  minus_x = 2,

  plus_y = 3,
  minus_y = 4,

  plus_z = 5,
  minus_z = 6,

  unknown = 7
};

static inline CCTK_HOST CCTK_DEVICE auto
get_owner_patch(const PatchParams &par,
                const svec_t &global_coords) -> PatchPiece {
  using std::distance;
  using std::fabs;
  using std::max_element;

  const auto x{global_coords(0)};
  const auto y{global_coords(1)};
  const auto z{global_coords(2)};

  const auto abs_x{fabs(x)};
  const auto abs_y{fabs(y)};
  const auto abs_z{fabs(z)};

  const auto r0{par.inner_boundary};

  if (abs_x <= r0 && abs_y <= r0 && abs_z <= r0) {
    return PatchPiece::cartesian;
  }

  std::array<CCTK_REAL, 3> abs_coords{abs_x, abs_y, abs_z};
  const auto max_coord_idx{distance(
      abs_coords.begin(), max_element(abs_coords.begin(), abs_coords.end()))};

  if (x > 0.0 && max_coord_idx == 0) {
    return PatchPiece::plus_x;
  }

  if (x < 0.0 && max_coord_idx == 0) {
    return PatchPiece::minus_x;
  }

  if (y > 0.0 && max_coord_idx == 1) {
    return PatchPiece::plus_y;
  }

  if (y < 0.0 && max_coord_idx == 1) {
    return PatchPiece::minus_y;
  }

  if (z > 0.0 && max_coord_idx == 2) {
    return PatchPiece::plus_z;
  }

  if (z < 0.0 && max_coord_idx == 2) {
    return PatchPiece::minus_z;
  }

  // We don't know where we are. This is unexpected
  CCTK_VINFO("Coordinate triplet (%.16f, %.16f, %.16f) cannot be located "
             "within the simulation domain",
             x, y, z);

  return PatchPiece::unknown;
}

CCTK_HOST CCTK_DEVICE auto
global2local(const PatchParams &par,
             const svec_t &global_coords) -> std_tuple<int, svec_t> {
  using std::pow;
  using std::sqrt;

  const auto r0{par.inner_boundary};
  const auto r1{par.outer_boundary};

  const auto x{global_coords(0)};
  const auto y{global_coords(1)};
  const auto z{global_coords(2)};

  const auto patch{get_owner_patch(par, global_coords)};

  svec_t local_coords{0.0, 0.0, 0.0};

  switch (patch) {

  case PatchPiece::cartesian:
    local_coords(0) = x / r0;
    local_coords(1) = y / r0;
    local_coords(2) = z / r0;
    break;

  case PatchPiece::plus_x:
    local_coords(0) = z / x;
    local_coords(1) = y / x;
    local_coords(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(y, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case PatchPiece::minus_x:
    local_coords(0) = -(z / x);
    local_coords(1) = y / x;
    local_coords(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(y, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case PatchPiece::plus_y:
    local_coords(0) = z / y;
    local_coords(1) = -(x / y);
    local_coords(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case PatchPiece::minus_y:
    local_coords(0) = -(z / y);
    local_coords(1) = -(x / y);
    local_coords(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case PatchPiece::plus_z:
    local_coords(0) = -(x / z);
    local_coords(1) = y / z;
    local_coords(2) = (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(y, 2) +
                       sqrt((pow(x, 2) + pow(y, 2)) *
                                (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                            4 * pow(r0 - r1, 2) * pow(z, 2))) /
                      pow(r0 - r1, 2);
    break;

  case PatchPiece::minus_z:
    local_coords(0) = -(x / z);
    local_coords(1) = -(y / z);
    local_coords(2) = (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(y, 2) +
                       sqrt((pow(x, 2) + pow(y, 2)) *
                                (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                            4 * pow(r0 - r1, 2) * pow(z, 2))) /
                      pow(r0 - r1, 2);
    break;

  default:
    CCTK_VERROR("Unable to compute global2local: Unknown patch piece");
    break;
  }

  return std_make_tuple(static_cast<int>(patch), local_coords);
}

CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &par, int patch,
                                        const svec_t &local_coords) -> svec_t {
  using std::pow;
  using std::sqrt;

  assert(0 <= patch && patch <= (static_cast<int>(PatchPiece::unknown) - 1));

  const auto r0{par.inner_boundary};
  const auto r1{par.outer_boundary};

  const auto a{local_coords(0)};
  const auto b{local_coords(1)};
  const auto c{local_coords(2)};

  svec_t global_coords = {0.0, 0.0, 0.0};

  switch (static_cast<PatchPiece>(patch)) {

  case PatchPiece::cartesian:
    global_coords(0) = a * r0;
    global_coords(1) = b * r0;
    global_coords(2) = c * r0;
    break;

  case PatchPiece::plus_x:
    global_coords(0) =
        (r0 - c * r0 + r1 + c * r1) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(1) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case PatchPiece::minus_x:
    global_coords(0) =
        -((r0 - c * r0 + r1 + c * r1) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_coords(1) =
        -((b * (r0 - c * r0 + r1 + c * r1)) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_coords(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case PatchPiece::plus_y:
    global_coords(0) =
        -((b * (r0 - c * r0 + r1 + c * r1)) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_coords(1) =
        (r0 - c * r0 + r1 + c * r1) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case PatchPiece::minus_y:
    global_coords(0) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(1) =
        -((r0 - c * r0 + r1 + c * r1) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_coords(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case PatchPiece::plus_z:
    global_coords(0) =
        -((a * (r0 - c * r0 + r1 + c * r1)) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_coords(1) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(2) =
        (r0 - c * r0 + r1 + c * r1) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case PatchPiece::minus_z:
    global_coords(0) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(1) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_coords(2) =
        -((r0 - c * r0 + r1 + c * r1) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    break;

  default:
    CCTK_VERROR("Unable to compute local2global: Unknown patch piece");
    break;
  }

  return global_coords;
}

static inline CCTK_HOST CCTK_DEVICE auto
cubed_sphere_jacs(const PatchParams &par, int patch,
                  const svec_t &global_coords) -> std_tuple<jac_t, djac_t> {
  using std::pow;
  using std::sqrt;

  assert(0 <= patch && patch <= (static_cast<int>(PatchPiece::unknown) - 1));

  jac_t J{};
  djac_t dJ{};

  const auto r0{par.inner_boundary};
  const auto r1{par.outer_boundary};

  const auto x{global_coords(0)};
  const auto y{global_coords(1)};
  const auto z{global_coords(2)};

  switch (static_cast<PatchPiece>(patch)) {

  case PatchPiece::cartesian:
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

  case PatchPiece::plus_x:
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

  case PatchPiece::minus_x:
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

  case PatchPiece::plus_y:
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

  case PatchPiece::minus_y:
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

  case PatchPiece::plus_z:
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

  case PatchPiece::minus_z:
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
    CCTK_VERROR("Unable to compute jacobians: Unknown patch piece");
    break;
  }

  return std_make_tuple(J, dJ);
}

CCTK_HOST CCTK_DEVICE auto
dlocal_dglobal(const PatchParams &par, int patch,
               const svec_t &local_coords) -> std_tuple<svec_t, jac_t> {
  const auto data{d2local_dglobal2(par, patch, local_coords)};
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &par, int patch,
                                            const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t> {
  const auto local_to_global_result{local2global(par, patch, local_coords)};
  const auto jacobian_results{
      cubed_sphere_jacs(par, patch, local_to_global_result)};

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
}

static inline auto make_patch(const PatchPiece &p,
                              const PatchParams &par) -> Patch {
  Patch patch{};

  switch (p) {

  case PatchPiece::cartesian:
    patch.name = "Cartesian";
    break;

  case PatchPiece::plus_x:
    patch.name = "Plus X";
    break;

  case PatchPiece::minus_x:
    patch.name = "Minus X";
    break;

  case PatchPiece::plus_y:
    patch.name = "Plus Y";
    break;

  case PatchPiece::minus_y:
    patch.name = "Minus Y";
    break;

  case PatchPiece::plus_z:
    patch.name = "Plus Z";
    break;

  case PatchPiece::minus_z:
    patch.name = "Minus Z";
    break;

  default:
    CCTK_VERROR("Unable to create patch. Unknown patch piece");
    break;
  }

  patch.ncells = {par.angular_cells, par.angular_cells, par.radial_cells};

  patch.xmin = {-1.0, -1.0, -1.0};
  patch.xmax = {1.0, 1.0, 1.0};

  patch.is_cartesian = false;

  PatchFace ob{true, -1};
  PatchFace co{false, static_cast<int>(PatchPiece::cartesian)};
  PatchFace px{false, static_cast<int>(PatchPiece::plus_x)};
  PatchFace mx{false, static_cast<int>(PatchPiece::minus_x)};
  PatchFace py{false, static_cast<int>(PatchPiece::plus_y)};
  PatchFace my{false, static_cast<int>(PatchPiece::minus_y)};
  PatchFace pz{false, static_cast<int>(PatchPiece::plus_z)};
  PatchFace mz{false, static_cast<int>(PatchPiece::minus_z)};

  if (p == PatchPiece::cartesian) {
    patch.ncells = {par.angular_cells, par.angular_cells, par.angular_cells};
    patch.is_cartesian = true;
    patch.faces = {{mx, my, mz}, {px, py, pz}};

  } else if (p == PatchPiece::plus_x) {
    patch.faces = {{mz, my, co}, {pz, py, ob}};

  } else if (p == PatchPiece::minus_x) {
    patch.faces = {{mz, py, co}, {pz, my, ob}};

  } else if (p == PatchPiece::plus_y) {
    patch.faces = {{mz, px, co}, {pz, mx, ob}};

  } else if (p == PatchPiece::minus_y) {
    patch.faces = {{mz, mx, co}, {pz, px, ob}};

  } else if (p == PatchPiece::plus_z) {
    patch.faces = {{px, my, co}, {mx, py, ob}};

  } else if (p == PatchPiece::minus_z) {
    patch.faces = {{mx, my, co}, {px, py, ob}};
  }

  return patch;
}

auto make_system(const PatchParams &par) -> PatchSystem {
  return PatchSystem{.name = "Cubed Sphere",
                     .id_tag = PatchSystems::cubed_spehre,
                     .patches = {make_patch(PatchPiece::cartesian, par),
                                 make_patch(PatchPiece::plus_x, par),
                                 make_patch(PatchPiece::minus_x, par),
                                 make_patch(PatchPiece::plus_y, par),
                                 make_patch(PatchPiece::minus_y, par),
                                 make_patch(PatchPiece::plus_z, par),
                                 make_patch(PatchPiece::minus_z, par)}};
}

template <typename fp_type>
static inline auto isapprox(fp_type x, fp_type y, fp_type atol = 0.0) -> bool {
  using std::abs;
  using std::max;
  using std::sqrt;

  const fp_type rtol{
      atol > 0.0 ? 0.0 : sqrt(std::numeric_limits<fp_type>::epsilon())};
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

auto unit_test(std::size_t repetitions, std::size_t seed,
               const PatchParams &par) -> bool {
  using std::cos;
  using std::sin;

  using real_dist = std::uniform_real_distribution<CCTK_REAL>;
  using int_dist = std::uniform_int_distribution<CCTK_INT>;

  std::mt19937 engine{static_cast<unsigned long>(seed)};

  real_dist r_dist{0.0, par.outer_boundary};
  real_dist theta_dist{0.0, M_PI};
  real_dist phi_dist{0.0, 2.0 * M_PI};

  real_dist local_dist{-1.0, 1.0};
  int_dist patch_dist{0, static_cast<CCTK_INT>(PatchPiece::unknown) - 1};

  bool all_pass{true};

  // local2global(global2local(global)) == global ?
  for (CCTK_INT i = 0; i < repetitions; i++) {
    const auto r{r_dist(engine)};
    const auto theta{theta_dist(engine)};
    const auto phi{phi_dist(engine)};

    const auto x{r * sin(theta) * cos(phi)};
    const auto y{r * sin(theta) * sin(phi)};
    const auto z{r * cos(theta)};

    const svec_t g_i{x, y, z};

    const auto l{global2local(par, g_i)};
    const auto g_f{local2global(par, std::get<0>(l), std::get<1>(l))};

    const auto passed{isapprox(g_i(0), g_f(0)) && isapprox(g_i(1), g_f(1)) &&
                      isapprox(g_i(2), g_f(2))};

    if (!passed) {
      CCTK_VINFO("local2global(global2local(global)) == global repetition %i "
                 "\033[1;31mFAILED\033[0m. Expected (%.16f, %.16f, %.16f) "
                 "but got (%.16f, %.16f, %.16f)",
                 i, g_i(0), g_i(1), g_i(2), g_f(0), g_f(1), g_f(2));
      all_pass = false;
    }
  }

  // global2local(local2global(local)) == global ?
  for (CCTK_INT i = 0; i < repetitions; i++) {
    const int p_i{patch_dist(engine)};
    const svec_t l_i{local_dist(engine), local_dist(engine),
                     local_dist(engine)};

    const auto g{local2global(par, p_i, l_i)};
    const auto l{global2local(par, g)};

    const auto &p_f{std::get<0>(l)};
    const auto &l_f{std::get<1>(l)};

    const auto passed{isapprox(l_i(0), l_f(0)) && isapprox(l_i(1), l_f(1)) &&
                      isapprox(l_i(2), l_f(2)) && p_i == p_f};

    if (!passed) {
      CCTK_VINFO(
          "global2local(local2global(local)) == global repetition %i "
          "\033[1;31mFAILED\033[0m. Expected (%.16f, %.16f, %.16f) patch %i "
          "but got (%.16f, %.16f, %.16f) patch %i",
          i, l_i(0), l_i(1), l_i(2), p_i, l_f(0), l_f(1), l_f(2), p_f);
      all_pass = false;
    }
  }

  return all_pass;
}

} // namespace CapyrX::MultiPatch::CubedSphere