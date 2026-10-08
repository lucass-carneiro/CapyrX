#include "llama.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>

namespace CapyrX::MultiPatch::Llama {

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

  // Ownership is split at the sphere r = R, reproducing real Llama's
  // global_to_local_Thornburg04 (Coordinates/src/thornburg04.cc: rp2 <
  // sphere_inner_radius^2 -> central cube). Only the inscribed sphere r < R is
  // cube-owned; the box-corner shell R < r < sqrt(3)R goes to the wedges, which
  // also hold interior cells there (overset). The cube patch grid is still the
  // full box [-R,R]^3 -- this classifies ownership, not grid extent. Compared
  // squared to avoid a sqrt.
  if (abs_x * abs_x + abs_y * abs_y + abs_z * abs_z < r0 * r0) {
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
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
  CCTK_VINFO("Coordinate triplet (%.16f, %.16f, %.16f) cannot be located "
             "within the simulation domain",
             x, y, z);
#else
  assert(false);
#endif

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
    local_coords(0) = x;
    local_coords(1) = y;
    local_coords(2) = z;
    break;

  case PatchPiece::plus_x:
    local_coords(0) = z / x;
    local_coords(1) = y / x;
    local_coords(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case PatchPiece::plus_y:
    local_coords(0) = z / y;
    local_coords(1) = -(x / y);
    local_coords(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case PatchPiece::minus_x:
    local_coords(0) = -(z / x);
    local_coords(1) = y / x;
    local_coords(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case PatchPiece::minus_y:
    local_coords(0) = -(z / y);
    local_coords(1) = -(x / y);
    local_coords(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case PatchPiece::plus_z:
    local_coords(0) = -(x / z);
    local_coords(1) = y / z;
    local_coords(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case PatchPiece::minus_z:
    local_coords(0) = -(x / z);
    local_coords(1) = -(y / z);
    local_coords(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to compute global2local: Unknown patch piece");
#else
    assert(false);
#endif
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
    global_coords(0) = a;
    global_coords(1) = b;
    global_coords(2) = c;
    break;

  case PatchPiece::plus_x:
    global_coords(0) =
        (r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(1) = (b * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case PatchPiece::plus_y:
    global_coords(0) = (b * ((-1 + c) * r0 - (1 + c) * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(1) =
        (r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case PatchPiece::minus_x:
    global_coords(0) =
        ((-1 + c) * r0 - (1 + c) * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(1) = (b * ((-1 + c) * r0 - (1 + c) * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case PatchPiece::minus_y:
    global_coords(0) = (b * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(1) =
        ((-1 + c) * r0 - (1 + c) * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case PatchPiece::plus_z:
    global_coords(0) = (a * ((-1 + c) * r0 - (1 + c) * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(1) = (b * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(2) =
        (r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case PatchPiece::minus_z:
    global_coords(0) = (a * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(1) = (b * (r0 - c * r0 + r1 + c * r1)) /
                       (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_coords(2) =
        ((-1 + c) * r0 - (1 + c) * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to compute local2global: Unknown patch piece");
#else
    assert(false);
#endif
    break;
  }

  return global_coords;
}

static inline CCTK_HOST CCTK_DEVICE auto
llama_jacs(const PatchParams &par, int patch, const svec_t &global_coords)
    -> std_tuple<jac_t, djac_t> {
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
    J(0)(0) = 1;
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)(0) = 0;
    J(1)(1) = 1;
    J(1)(2) = 0;
    J(2)(0) = 0;
    J(2)(1) = 0;
    J(2)(2) = 1;

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
    J(2)(0) = (-2 * x) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(1) = (-2 * y) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(2) = (-2 * z) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

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
    dJ(2)(0, 0) = (-2 * (pow(y, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 1) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 2) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 0) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 1) = (-2 * (pow(x, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 2) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 0) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 1) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 2) = (-2 * (pow(x, 2) + pow(y, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    break;

  case PatchPiece::plus_y:
    J(0)(0) = 0;
    J(0)(1) = -(z / pow(y, 2));
    J(0)(2) = 1 / y;
    J(1)(0) = -(1 / y);
    J(1)(1) = x / pow(y, 2);
    J(1)(2) = 0;
    J(2)(0) = (-2 * x) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(1) = (-2 * y) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(2) = (-2 * z) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

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
    dJ(2)(0, 0) = (-2 * (pow(y, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 1) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 2) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 0) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 1) = (-2 * (pow(x, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 2) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 0) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 1) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 2) = (-2 * (pow(x, 2) + pow(y, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    break;

  case PatchPiece::minus_x:
    J(0)(0) = z / pow(x, 2);
    J(0)(1) = 0;
    J(0)(2) = -(1 / x);
    J(1)(0) = -(y / pow(x, 2));
    J(1)(1) = 1 / x;
    J(1)(2) = 0;
    J(2)(0) = (-2 * x) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(1) = (-2 * y) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(2) = (-2 * z) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

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
    dJ(2)(0, 0) = (-2 * (pow(y, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 1) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 2) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 0) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 1) = (-2 * (pow(x, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 2) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 0) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 1) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 2) = (-2 * (pow(x, 2) + pow(y, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    break;

  case PatchPiece::minus_y:
    J(0)(0) = 0;
    J(0)(1) = z / pow(y, 2);
    J(0)(2) = -(1 / y);
    J(1)(0) = -(1 / y);
    J(1)(1) = x / pow(y, 2);
    J(1)(2) = 0;
    J(2)(0) = (-2 * x) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(1) = (-2 * y) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(2) = (-2 * z) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

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
    dJ(2)(0, 0) = (-2 * (pow(y, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 1) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 2) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 0) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 1) = (-2 * (pow(x, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 2) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 0) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 1) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 2) = (-2 * (pow(x, 2) + pow(y, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    break;

  case PatchPiece::plus_z:
    J(0)(0) = -(1 / z);
    J(0)(1) = 0;
    J(0)(2) = x / pow(z, 2);
    J(1)(0) = 0;
    J(1)(1) = 1 / z;
    J(1)(2) = -(y / pow(z, 2));
    J(2)(0) = (-2 * x) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(1) = (-2 * y) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(2) = (-2 * z) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

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
    dJ(2)(0, 0) = (-2 * (pow(y, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 1) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 2) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 0) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 1) = (-2 * (pow(x, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 2) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 0) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 1) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 2) = (-2 * (pow(x, 2) + pow(y, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    break;

  case PatchPiece::minus_z:
    J(0)(0) = -(1 / z);
    J(0)(1) = 0;
    J(0)(2) = x / pow(z, 2);
    J(1)(0) = 0;
    J(1)(1) = -(1 / z);
    J(1)(2) = y / pow(z, 2);
    J(2)(0) = (-2 * x) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(1) = (-2 * y) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    J(2)(2) = (-2 * z) / ((r0 - r1) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

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
    dJ(2)(0, 0) = (-2 * (pow(y, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 1) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(0, 2) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 0) =
        (2 * x * y) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 1) = (-2 * (pow(x, 2) + pow(z, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(1, 2) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 0) =
        (2 * x * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 1) =
        (2 * y * z) / ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    dJ(2)(2, 2) = (-2 * (pow(x, 2) + pow(y, 2))) /
                  ((r0 - r1) * pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to compute jacobians: Unknown patch piece");
#else
    assert(false);
#endif
    break;
  }

  return std_make_tuple(J, dJ);
}

CCTK_HOST CCTK_DEVICE auto dlocal_dglobal(const PatchParams &par, int patch,
                                          const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t> {
  const auto data{d2local_dglobal2(par, patch, local_coords)};
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &par, int patch,
                                            const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t> {
  const auto local_to_global_result{local2global(par, patch, local_coords)};
  const auto jacobian_results{llama_jacs(par, patch, local_to_global_result)};

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
}

static inline auto make_patch(const PatchPiece &p,
                              const PatchParams &par) -> Patch {

  const auto twice_overlap = 2 * par.patch_overlap;
  const CCTK_REAL angular_delta = 2.0 / par.angular_cells;
  const CCTK_REAL radial_delta = 2.0 / par.radial_cells;

  // Default: a wedge. Angular faces are interpatch (need overlap for donor
  // stencils); the radial direction is c. The inner radial face is fed by the
  // cube (co), the outer radial face is the physical outer boundary (ob).
  Patch patch{};

  patch.ncells = {par.angular_cells + twice_overlap,
                  par.angular_cells + twice_overlap,
                  par.radial_cells + par.patch_overlap};

  patch.xmin = {
      CCTK_REAL{-1.0} - par.patch_overlap * angular_delta,
      CCTK_REAL{-1.0} - par.patch_overlap * angular_delta,
      CCTK_REAL{-1.0} - par.patch_overlap * radial_delta,
  };

  patch.xmax = {
      CCTK_REAL{1.0} + par.patch_overlap * angular_delta,
      CCTK_REAL{1.0} + par.patch_overlap * angular_delta,
      CCTK_REAL{1.0},
  };

  patch.is_cartesian = false;
  patch.c_is_radial = true;

  PatchFace ob{true, -1};
  PatchFace co{false, static_cast<int>(PatchPiece::cartesian)};
  PatchFace px{false, static_cast<int>(PatchPiece::plus_x)};
  PatchFace mx{false, static_cast<int>(PatchPiece::minus_x)};
  PatchFace py{false, static_cast<int>(PatchPiece::plus_y)};
  PatchFace my{false, static_cast<int>(PatchPiece::minus_y)};
  PatchFace pz{false, static_cast<int>(PatchPiece::plus_z)};
  PatchFace mz{false, static_cast<int>(PatchPiece::minus_z)};

  switch (p) {

  case PatchPiece::cartesian: {
    patch.name = "Cartesian";

    // The cube is a true [-R,R]^3 Cartesian patch whose resolution is
    // independent of the angular grid (Llama's h_cartesian), so it must not
    // reuse angular_cells the way cubed_sphere does.
    //
    // All 6 faces are interpatch, so the cube carries the evolved overlap band
    // on every side (like the wedges): extend each extent by patch_overlap cells
    // and add 2*patch_overlap cells, keeping dx fixed. Without this the cube's
    // evolved interior stops at r = R and cannot donate a centered stencil for a
    // wedge's near-seam inner ghost (the R2 donor dropout).
    const CCTK_REAL cube_delta_i = 2.0 * par.inner_boundary / par.cube_ncells_i;
    const CCTK_REAL cube_delta_j = 2.0 * par.inner_boundary / par.cube_ncells_j;
    const CCTK_REAL cube_delta_k = 2.0 * par.inner_boundary / par.cube_ncells_k;

    patch.ncells = {par.cube_ncells_i + twice_overlap,
                    par.cube_ncells_j + twice_overlap,
                    par.cube_ncells_k + twice_overlap};

    patch.xmin = {
        -par.inner_boundary - par.patch_overlap * cube_delta_i,
        -par.inner_boundary - par.patch_overlap * cube_delta_j,
        -par.inner_boundary - par.patch_overlap * cube_delta_k,
    };

    patch.xmax = {
        par.inner_boundary + par.patch_overlap * cube_delta_i,
        par.inner_boundary + par.patch_overlap * cube_delta_j,
        par.inner_boundary + par.patch_overlap * cube_delta_k,
    };

    patch.faces = {{mx, my, mz}, {px, py, pz}};

    patch.is_cartesian = true;
    patch.c_is_radial = false;

    break;
  }

  case PatchPiece::plus_x:
    patch.name = "Plus X";
    patch.faces = {{mz, my, co}, {pz, py, ob}};
    break;

  case PatchPiece::minus_x:
    patch.name = "Minus X";
    patch.faces = {{mz, py, co}, {pz, my, ob}};
    break;

  case PatchPiece::plus_y:
    patch.name = "Plus Y";
    patch.faces = {{mz, px, co}, {pz, mx, ob}};
    break;

  case PatchPiece::minus_y:
    patch.name = "Minus Y";
    patch.faces = {{mz, mx, co}, {pz, px, ob}};
    break;

  case PatchPiece::plus_z:
    patch.name = "Plus Z";
    patch.faces = {{px, my, co}, {mx, py, ob}};
    break;

  case PatchPiece::minus_z:
    patch.name = "Minus Z";
    patch.faces = {{mx, my, co}, {px, py, ob}};
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to create patch. Unknown patch piece");
#else
    assert(false);
#endif
    break;
  }

  return patch;
}

auto make_system(const PatchParams &par) -> PatchSystem {
  return PatchSystem{.name = "Llama",
                     .id_tag = PatchSystems::llama,
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

  std::mt19937 engine{static_cast<std::mt19937>(seed)};

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

  // global2local(local2global(local)) round-trip. The patch/local identity only
  // holds in single-covered regions: a cube local point in the box-corner shell
  // (R<r<sqrt(3)R) is wedge-owned under the spherical rule (overset, design
  // S3.2). There the position must still round-trip, but ownership legitimately
  // switches to a wedge, so assert identity only when the generated point is
  // single-covered by the sampled patch.
  //
  // Both branches must actually be exercised, else the test passes vacuously:
  // count single-cover hits and overset corner-shell hits (a cube sample
  // reclassified to a wedge) and require each to be non-zero below.
  std::size_t single_cover_hits{0};
  std::size_t overset_shell_hits{0};
  for (CCTK_INT i = 0; i < repetitions; i++) {
    const int p_i{patch_dist(engine)};
    const svec_t l_i{local_dist(engine), local_dist(engine),
                     local_dist(engine)};

    const auto g{local2global(par, p_i, l_i)};
    const auto owner{get_owner_patch(par, g)};

    const auto l{global2local(par, g)};
    const auto &p_f{std::get<0>(l)};
    const auto &l_f{std::get<1>(l)};

    const auto g_rt{local2global(par, p_f, l_f)};

    bool passed{isapprox(g(0), g_rt(0)) && isapprox(g(1), g_rt(1)) &&
                isapprox(g(2), g_rt(2))};

    const bool single_covered{static_cast<int>(owner) == p_i};
    if (single_covered) {
      ++single_cover_hits;
      passed = passed && p_i == p_f && isapprox(l_i(0), l_f(0)) &&
               isapprox(l_i(1), l_f(1)) && isapprox(l_i(2), l_f(2));
    } else if (p_i == static_cast<int>(PatchPiece::cartesian) &&
               owner != PatchPiece::cartesian) {
      // A cube sample in the box-corner shell (R < r < sqrt(3)R) reclassified to
      // a wedge -- the overset region under the spherical ownership rule.
      ++overset_shell_hits;
    }

    if (!passed) {
      CCTK_VINFO(
          "global2local(local2global(local)) round-trip repetition %i "
          "\033[1;31mFAILED\033[0m. Patch %i local (%.16f, %.16f, %.16f) -> "
          "global (%.16f, %.16f, %.16f); got patch %i local "
          "(%.16f, %.16f, %.16f), owner patch %i",
          i, p_i, l_i(0), l_i(1), l_i(2), g(0), g(1), g(2), p_f, l_f(0), l_f(1),
          l_f(2), static_cast<int>(owner));
      all_pass = false;
    }
  }

  if (single_cover_hits == 0) {
    CCTK_VINFO("Round-trip test exercised no single-covered samples "
               "\033[1;31m(vacuous)\033[0m. Increase repetitions.");
    all_pass = false;
  }
  if (overset_shell_hits == 0) {
    CCTK_VINFO("Round-trip test exercised no overset corner-shell samples "
               "\033[1;31m(vacuous)\033[0m. Increase repetitions or check that "
               "the geometry has a box-corner shell (R < r < sqrt(3)R).");
    all_pass = false;
  }

  // Ownership spot-checks (design S3.2).
  // On-axis wedge midpoint: single-covered, owned by the +x wedge.
  {
    const svec_t global_coords{par.inner_boundary +
                                   (par.outer_boundary - par.inner_boundary) /
                                       CCTK_REAL{2.0},
                               CCTK_REAL{0.0}, CCTK_REAL{0.0}};
    const auto owner{get_owner_patch(par, global_coords)};

    if (owner != PatchPiece::plus_x) {
      CCTK_VINFO("Wedge-owner spot-check failed. Expected patch %i but got %i",
                 static_cast<int>(PatchPiece::plus_x),
                 static_cast<int>(owner));
      all_pass = false;
    }
  }

  // Cube interior (r < R): cube-owned.
  {
    const auto half{par.inner_boundary / CCTK_REAL{2.0}};
    const svec_t global_coords{half, half, half};
    const auto owner{get_owner_patch(par, global_coords)};

    if (owner != PatchPiece::cartesian) {
      CCTK_VINFO("Cube-interior spot-check failed. Expected patch %i but got %i",
                 static_cast<int>(PatchPiece::cartesian),
                 static_cast<int>(owner));
      all_pass = false;
    }
  }

  // Box-corner shell: (0.9R, 0.9R, 0.9R) has r = 1.56R > R but all |coord| < R.
  // Under the spherical ownership rule (matching Llama) this is wedge-owned
  // (+x by the x>=y>=z tie-break), NOT cube-owned. This pins the overset rule of
  // design S3.2 and guards against a regression back to the box classifier.
  {
    const auto s{CCTK_REAL{0.9} * par.inner_boundary};
    const svec_t global_coords{s, s, s};
    const auto owner{get_owner_patch(par, global_coords)};

    if (owner != PatchPiece::plus_x) {
      CCTK_VINFO("Box-corner-shell spot-check failed. Expected patch %i (+x "
                 "wedge) but got %i",
                 static_cast<int>(PatchPiece::plus_x),
                 static_cast<int>(owner));
      all_pass = false;
    }
  }

  // Analytic-vs-finite-difference Jacobian check (design S6). The round-trips
  // above validate only the coordinate maps; nothing exercises the reused
  // thornburg06 wedge Jacobians J and dJ. Verify them against central
  // differences of local2global using two identities that hold on every patch
  // and need no matrix inverse:
  //   (1) J . M = I, with M(k,m) = d global_k / d local_m (FD of local2global)
  //       and J(i)(k) = d local_i / d global_k (analytic).
  //   (2) d J(i)(j) / d local_m = sum_k dJ(i)(j,k) . M(k,m) (chain rule): the
  //       left side is a central difference of J, the right contracts analytic
  //       dJ with M.
  // d2local_dglobal2 takes an explicit patch, so no owner dispatch happens and
  // the overset shell is irrelevant here. Step size and tolerance scale with
  // CCTK_REAL precision so the check is meaningful in real64 and real32 builds.
  {
    using std::abs;
    using std::max;
    using std::pow;

    const CCTK_REAL eps{std::numeric_limits<CCTK_REAL>::epsilon()};
    const CCTK_REAL h{pow(eps, CCTK_REAL{1} / 3)};
    const CCTK_REAL tol{CCTK_REAL{1000} * pow(eps, CCTK_REAL{2} / 3)};

    const auto close{[&](CCTK_REAL a, CCTK_REAL b) {
      return abs(a - b) <= tol + tol * max(abs(a), abs(b));
    }};

    // Keep samples away from the local cube edges so the FD error stays small.
    real_dist interior_dist{-0.7, 0.7};
    std::size_t jac_checks{0};

    for (int p = 0; p <= static_cast<int>(PatchPiece::unknown) - 1; ++p) {
      for (CCTK_INT rep = 0; rep < repetitions; ++rep) {
        const svec_t l0{interior_dist(engine), interior_dist(engine),
                        interior_dist(engine)};

        const auto base{d2local_dglobal2(par, p, l0)};
        const auto &J0{std::get<1>(base)};
        const auto &dJ0{std::get<2>(base)};

        ++jac_checks;

        for (int m = 0; m < 3; ++m) {
          svec_t lp{l0}, lm{l0};
          lp(m) += h;
          lm(m) -= h;

          const auto gp{local2global(par, p, lp)};
          const auto gm{local2global(par, p, lm)};
          const auto base_p{d2local_dglobal2(par, p, lp)};
          const auto base_m{d2local_dglobal2(par, p, lm)};
          const auto &Jp{std::get<1>(base_p)};
          const auto &Jm{std::get<1>(base_m)};

          CCTK_REAL M[3];
          for (int k = 0; k < 3; ++k) {
            M[k] = (gp(k) - gm(k)) / (2 * h);
          }

          // Identity (1), column m of J.M.
          for (int i = 0; i < 3; ++i) {
            CCTK_REAL jm{0};
            for (int k = 0; k < 3; ++k) {
              jm += J0(i)(k) * M[k];
            }
            const CCTK_REAL expected{i == m ? CCTK_REAL{1} : CCTK_REAL{0}};
            if (!close(jm, expected)) {
              CCTK_VINFO("Jacobian J.M = I check \033[1;31mFAILED\033[0m on "
                         "patch %i, entry (%i,%i): got %.16e, expected %.16e",
                         p, i, m, jm, expected);
              all_pass = false;
            }
          }

          // Identity (2), column m of the chain rule.
          for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
              const CCTK_REAL L{(Jp(i)(j) - Jm(i)(j)) / (2 * h)};
              CCTK_REAL rhs{0};
              for (int k = 0; k < 3; ++k) {
                rhs += dJ0(i)(j, k) * M[k];
              }
              if (!close(L, rhs)) {
                CCTK_VINFO("Jacobian-derivative dJ check \033[1;31mFAILED\033[0m "
                           "on patch %i, entry (%i,%i,%i): FD %.16e, "
                           "analytic %.16e",
                           p, i, j, m, L, rhs);
                all_pass = false;
              }
            }
          }
        }
      }
    }

    if (jac_checks == 0) {
      CCTK_VINFO("Jacobian FD check exercised no samples "
                 "\033[1;31m(vacuous)\033[0m. Increase repetitions.");
      all_pass = false;
    }
  }

  return all_pass;
}

} // namespace CapyrX::MultiPatch::Llama
