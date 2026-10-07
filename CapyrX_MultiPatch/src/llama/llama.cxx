#include "llama.hxx"

#include <cmath>

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

} // namespace CapyrX::MultiPatch::Llama
