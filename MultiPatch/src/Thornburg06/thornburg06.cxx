#include "thornburg06.hxx"

namespace MultiPatch {
namespace Thornburg06 {

std::string piece_name(const patch_piece &p) {
  switch (p) {
  case plus_x:
    return "plus_x";
    break;
  case plus_y:
    return "plus_y";
    break;
  case minus_x:
    return "minus_x";
    break;
  case minus_y:
    return "minus_y";
    break;
  case plus_z:
    return "plus_z";
    break;
  case minus_z:
    return "minus_z";
    break;
  default:
    return "unknown";
    break;
  }
}

CCTK_DEVICE CCTK_HOST patch_piece
get_owner_patch(const PatchTransformations &pt, const svec &global_vars) {
  using std::distance;
  using std::fabs;
  using std::max_element;

  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  const auto abs_x{fabs(x)};
  const auto abs_y{fabs(y)};
  const auto abs_z{fabs(z)};

  std::array<CCTK_REAL, 3> abs_coords{abs_x, abs_y, abs_z};
  const auto max_coord_idx{distance(
      abs_coords.begin(), max_element(abs_coords.begin(), abs_coords.end()))};

  if (x > 0.0 && max_coord_idx == 0) {
    return patch_piece::plus_x;
  }

  if (x < 0.0 && max_coord_idx == 0) {
    return patch_piece::minus_x;
  }

  if (y > 0.0 && max_coord_idx == 1) {
    return patch_piece::plus_y;
  }

  if (y < 0.0 && max_coord_idx == 1) {
    return patch_piece::minus_y;
  }

  if (z > 0.0 && max_coord_idx == 2) {
    return patch_piece::plus_z;
  }

  if (z < 0.0 && max_coord_idx == 2) {
    return patch_piece::minus_z;
  }

// We don't know where we are. This is unexpected
#ifndef __CUDACC__
  CCTK_VERROR("Coordinate triplet (%f, %f, %f) cannot be located within "
              "the simulation domain",
              x, y, z);
#else
  assert(0);
#endif
  return patch_piece::unknown;
}

CCTK_DEVICE CCTK_HOST svec local2global(const PatchTransformations &pt,
                                        int patch, const svec &local_vars) {
  using std::pow;
  using std::sqrt;

  const auto r0{pt.thornburg06_inner_boundary_radius};
  const auto r1{pt.thornburg06_outer_boundary_radius};

  const auto a{local_vars(0)};
  const auto b{local_vars(1)};
  const auto c{local_vars(2)};

  svec global_vars = {0.0, 0.0, 0.0};

  switch (patch) {

  case patch_piece::plus_x:
    global_vars(0) =
        (r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(1) = (b * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case patch_piece::plus_y:
    global_vars(0) = -(b * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(1) =
        (r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case patch_piece::minus_x:
    global_vars(0) =
        -(r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(1) = -(b * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case patch_piece::minus_y:
    global_vars(0) = (b * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(1) =
        -(r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(2) = (a * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case patch_piece::plus_z:
    global_vars(0) = -(a * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(1) = (b * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(2) =
        (r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  case patch_piece::minus_z:
    global_vars(0) = (a * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(1) = (b * (r0 - c * r0 + r1 + c * r1)) /
                     (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    global_vars(2) =
        -(r0 - c * r0 + r1 + c * r1) / (2. * sqrt(1 + pow(a, 2) + pow(b, 2)));
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No local -> global transformations available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
    break;
  }

  return global_vars;
}

CCTK_DEVICE CCTK_HOST std_tuple<int, svec>
global2local(const PatchTransformations &pt, const svec &global_vars) {
  using std::pow;
  using std::sqrt;

  const auto r0{pt.thornburg06_inner_boundary_radius};
  const auto r1{pt.thornburg06_outer_boundary_radius};

  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  const auto patch{get_owner_patch(pt, global_vars)};

  svec local_vars{0.0, 0.0, 0.0};

  switch (patch) {

  case patch_piece::plus_x:
    local_vars(0) = z / x;
    local_vars(1) = y / x;
    local_vars(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case patch_piece::plus_y:
    local_vars(0) = z / y;
    local_vars(1) = -(x / y);
    local_vars(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case patch_piece::minus_x:
    local_vars(0) = -(z / x);
    local_vars(1) = y / x;
    local_vars(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case patch_piece::minus_y:
    local_vars(0) = -(z / y);
    local_vars(1) = -(x / y);
    local_vars(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case patch_piece::plus_z:
    local_vars(0) = -(x / z);
    local_vars(1) = y / z;
    local_vars(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  case patch_piece::minus_z:
    local_vars(0) = -(x / z);
    local_vars(1) = -(y / z);
    local_vars(2) =
        (r0 + r1 - 2 * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))) / (r0 - r1);
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("At point (%f, %f, %f): No global -> local transformations "
                "available for patch %s.",
                x, y, z, piece_name(patch).c_str());
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(static_cast<int>(patch), local_vars);
}

CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t, djac_t>
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const svec &local_vars) {
  using std::pow;

  jac_t J{};
  djac_t dJ{};

  const auto global_vars{local2global(pt, patch, local_vars)};

  const auto r0{pt.thornburg06_inner_boundary_radius};
  const auto r1{pt.thornburg06_outer_boundary_radius};

  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  // clang-format off
  switch (patch) {
  case patch_piece::plus_x:
    J(0)(0) = -(z/pow(x,2));
    J(0)(1) = 0;
    J(0)(2) = 1/x;
    J(1)(0) = -(y/pow(x,2));
    J(1)(1) = 1/x;
    J(1)(2) = 0;
    J(2)(0) = (-2*x)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(1) = (-2*y)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(2) = (-2*z)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

    dJ(0)(0,0) = (2*z)/pow(x,3);
    dJ(0)(0,1) = 0;
    dJ(0)(0,2) = -pow(x,-2);
    dJ(0)(1,0) = 0;
    dJ(0)(1,1) = 0;
    dJ(0)(1,2) = 0;
    dJ(0)(2,0) = -pow(x,-2);
    dJ(0)(2,1) = 0;
    dJ(0)(2,2) = 0;
    dJ(1)(0,0) = (2*y)/pow(x,3);
    dJ(1)(0,1) = -pow(x,-2);
    dJ(1)(0,2) = 0;
    dJ(1)(1,0) = -pow(x,-2);
    dJ(1)(1,1) = 0;
    dJ(1)(1,2) = 0;
    dJ(1)(2,0) = 0;
    dJ(1)(2,1) = 0;
    dJ(1)(2,2) = 0;
    dJ(2)(0,0) = (-2*(pow(y,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,1) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,2) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,0) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,1) = (-2*(pow(x,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,2) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,0) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,1) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,2) = (-2*(pow(x,2) + pow(y,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    break;

  case patch_piece::plus_y:
    J(0)(0) = 0;
    J(0)(1) = -(z/pow(y,2));
    J(0)(2) = 1/y;
    J(1)(0) = -(1/y);
    J(1)(1) = x/pow(y,2);
    J(1)(2) = 0;
    J(2)(0) = (-2*x)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(1) = (-2*y)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(2) = (-2*z)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

    dJ(0)(0,0) = 0;
    dJ(0)(0,1) = 0;
    dJ(0)(0,2) = 0;
    dJ(0)(1,0) = 0;
    dJ(0)(1,1) = (2*z)/pow(y,3);
    dJ(0)(1,2) = -pow(y,-2);
    dJ(0)(2,0) = 0;
    dJ(0)(2,1) = -pow(y,-2);
    dJ(0)(2,2) = 0;
    dJ(1)(0,0) = 0;
    dJ(1)(0,1) = pow(y,-2);
    dJ(1)(0,2) = 0;
    dJ(1)(1,0) = pow(y,-2);
    dJ(1)(1,1) = (-2*x)/pow(y,3);
    dJ(1)(1,2) = 0;
    dJ(1)(2,0) = 0;
    dJ(1)(2,1) = 0;
    dJ(1)(2,2) = 0;
    dJ(2)(0,0) = (-2*(pow(y,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,1) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,2) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,0) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,1) = (-2*(pow(x,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,2) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,0) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,1) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,2) = (-2*(pow(x,2) + pow(y,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    break;

  case patch_piece::minus_x:
    J(0)(0) = z/pow(x,2);
    J(0)(1) = 0;
    J(0)(2) = -(1/x);
    J(1)(0) = -(y/pow(x,2));
    J(1)(1) = 1/x;
    J(1)(2) = 0;
    J(2)(0) = (-2*x)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(1) = (-2*y)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(2) = (-2*z)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

    dJ(0)(0,0) = (-2*z)/pow(x,3);
    dJ(0)(0,1) = 0;
    dJ(0)(0,2) = pow(x,-2);
    dJ(0)(1,0) = 0;
    dJ(0)(1,1) = 0;
    dJ(0)(1,2) = 0;
    dJ(0)(2,0) = pow(x,-2);
    dJ(0)(2,1) = 0;
    dJ(0)(2,2) = 0;
    dJ(1)(0,0) = (2*y)/pow(x,3);
    dJ(1)(0,1) = -pow(x,-2);
    dJ(1)(0,2) = 0;
    dJ(1)(1,0) = -pow(x,-2);
    dJ(1)(1,1) = 0;
    dJ(1)(1,2) = 0;
    dJ(1)(2,0) = 0;
    dJ(1)(2,1) = 0;
    dJ(1)(2,2) = 0;
    dJ(2)(0,0) = (-2*(pow(y,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,1) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,2) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,0) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,1) = (-2*(pow(x,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,2) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,0) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,1) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,2) = (-2*(pow(x,2) + pow(y,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    break;

  case patch_piece::minus_y:
    J(0)(0) = 0;
    J(0)(1) = z/pow(y,2);
    J(0)(2) = -(1/y);
    J(1)(0) = -(1/y);
    J(1)(1) = x/pow(y,2);
    J(1)(2) = 0;
    J(2)(0) = (-2*x)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(1) = (-2*y)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(2) = (-2*z)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

    dJ(0)(0,0) = 0;
    dJ(0)(0,1) = 0;
    dJ(0)(0,2) = 0;
    dJ(0)(1,0) = 0;
    dJ(0)(1,1) = (-2*z)/pow(y,3);
    dJ(0)(1,2) = pow(y,-2);
    dJ(0)(2,0) = 0;
    dJ(0)(2,1) = pow(y,-2);
    dJ(0)(2,2) = 0;
    dJ(1)(0,0) = 0;
    dJ(1)(0,1) = pow(y,-2);
    dJ(1)(0,2) = 0;
    dJ(1)(1,0) = pow(y,-2);
    dJ(1)(1,1) = (-2*x)/pow(y,3);
    dJ(1)(1,2) = 0;
    dJ(1)(2,0) = 0;
    dJ(1)(2,1) = 0;
    dJ(1)(2,2) = 0;
    dJ(2)(0,0) = (-2*(pow(y,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,1) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,2) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,0) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,1) = (-2*(pow(x,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,2) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,0) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,1) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,2) = (-2*(pow(x,2) + pow(y,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    break;

  case patch_piece::plus_z:
    J(0)(0) = -(1/z);
    J(0)(1) = 0;
    J(0)(2) = x/pow(z,2);
    J(1)(0) = 0;
    J(1)(1) = 1/z;
    J(1)(2) = -(y/pow(z,2));
    J(2)(0) = (-2*x)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(1) = (-2*y)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(2) = (-2*z)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

    dJ(0)(0,0) = 0;
    dJ(0)(0,1) = 0;
    dJ(0)(0,2) = pow(z,-2);
    dJ(0)(1,0) = 0;
    dJ(0)(1,1) = 0;
    dJ(0)(1,2) = 0;
    dJ(0)(2,0) = pow(z,-2);
    dJ(0)(2,1) = 0;
    dJ(0)(2,2) = (-2*x)/pow(z,3);
    dJ(1)(0,0) = 0;
    dJ(1)(0,1) = 0;
    dJ(1)(0,2) = 0;
    dJ(1)(1,0) = 0;
    dJ(1)(1,1) = 0;
    dJ(1)(1,2) = -pow(z,-2);
    dJ(1)(2,0) = 0;
    dJ(1)(2,1) = -pow(z,-2);
    dJ(1)(2,2) = (2*y)/pow(z,3);
    dJ(2)(0,0) = (-2*(pow(y,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,1) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,2) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,0) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,1) = (-2*(pow(x,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,2) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,0) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,1) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,2) = (-2*(pow(x,2) + pow(y,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    break;

  case patch_piece::minus_z:
    J(0)(0) = -(1/z);
    J(0)(1) = 0;
    J(0)(2) = x/pow(z,2);
    J(1)(0) = 0;
    J(1)(1) = -(1/z);
    J(1)(2) = y/pow(z,2);
    J(2)(0) = (-2*x)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(1) = (-2*y)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    J(2)(2) = (-2*z)/((r0 - r1)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));

    dJ(0)(0,0) = 0;
    dJ(0)(0,1) = 0;
    dJ(0)(0,2) = pow(z,-2);
    dJ(0)(1,0) = 0;
    dJ(0)(1,1) = 0;
    dJ(0)(1,2) = 0;
    dJ(0)(2,0) = pow(z,-2);
    dJ(0)(2,1) = 0;
    dJ(0)(2,2) = (-2*x)/pow(z,3);
    dJ(1)(0,0) = 0;
    dJ(1)(0,1) = 0;
    dJ(1)(0,2) = 0;
    dJ(1)(1,0) = 0;
    dJ(1)(1,1) = 0;
    dJ(1)(1,2) = pow(z,-2);
    dJ(1)(2,0) = 0;
    dJ(1)(2,1) = pow(z,-2);
    dJ(1)(2,2) = (-2*y)/pow(z,3);
    dJ(2)(0,0) = (-2*(pow(y,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,1) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(0,2) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,0) = (2*x*y)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,1) = (-2*(pow(x,2) + pow(z,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(1,2) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,0) = (2*x*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,1) = (2*y*z)/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    dJ(2)(2,2) = (-2*(pow(x,2) + pow(y,2)))/((r0 - r1)*pow(pow(x,2) + pow(y,2) + pow(z,2),1.5));
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No jacobians available for patch %s", piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
    break;
  }
  // clang-format on

  return std_make_tuple(global_vars, J, dJ);
}

CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t>
dlocal_dglobal(const PatchTransformations &pt, int patch,
               const svec &local_vars) {
  const auto data = pt.d2local_dglobal2(pt, patch, local_vars);
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

template <patch_piece p> Patch make_patch(const PatchTransformations &pt) {
  Patch patch{};
  patch.name = piece_name(p);

  patch.ncells = {pt.thornburg06_angular_cells, pt.thornburg06_angular_cells,
                  pt.thornburg06_radial_cells};

  patch.xmin = {-1.0, -1.0, -1.0};
  patch.xmax = {1.0, 1.0, 1.0};

  patch.is_cartesian = false;

  const PatchFace px{false, patch_piece::plus_x};
  const PatchFace py{false, patch_piece::plus_y};
  const PatchFace mx{false, patch_piece::minus_x};
  const PatchFace my{false, patch_piece::minus_y};
  const PatchFace pz{false, patch_piece::plus_z};
  const PatchFace mz{false, patch_piece::minus_z};
  const PatchFace ob{true, -1};

  if constexpr (p == patch_piece::plus_x) {
    patch.faces = {{mz, my, ob}, {pz, py, ob}};
  } else if constexpr (p == patch_piece::plus_y) {
    patch.faces = {{mz, px, ob}, {pz, mx, ob}};
  } else if constexpr (p == patch_piece::minus_x) {
    patch.faces = {{mz, py, ob}, {pz, my, ob}};
  } else if constexpr (p == patch_piece::minus_y) {
    patch.faces = {{mz, mx, ob}, {pz, px, ob}};
  } else if constexpr (p == patch_piece::plus_z) {
    patch.faces = {{px, my, ob}, {mx, py, ob}};
  } else if constexpr (p == patch_piece::minus_z) {
    patch.faces = {{mx, my, ob}, {px, py, ob}};
  }

  return patch;
}

} // namespace Thornburg06

PatchSystem SetupThornburg06() {
  PatchTransformations pt;
  pt.global2local = &Thornburg06::global2local;
  pt.local2global = &Thornburg06::local2global;
  pt.dlocal_dglobal = &Thornburg06::dlocal_dglobal;
  pt.d2local_dglobal2 = &Thornburg06::d2local_dglobal2;
  pt.global2local_device = &Thornburg06::global2local;
  pt.local2global_device = &Thornburg06::local2global;
  pt.dlocal_dglobal_device = &Thornburg06::dlocal_dglobal;
  pt.d2local_dglobal2_device = &Thornburg06::d2local_dglobal2;

  const auto patches = std::vector<Patch>{
      Thornburg06::make_patch<Thornburg06::patch_piece::plus_x>(pt),
      Thornburg06::make_patch<Thornburg06::patch_piece::plus_y>(pt),
      Thornburg06::make_patch<Thornburg06::patch_piece::minus_x>(pt),
      Thornburg06::make_patch<Thornburg06::patch_piece::minus_y>(pt),
      Thornburg06::make_patch<Thornburg06::patch_piece::plus_z>(pt),
      Thornburg06::make_patch<Thornburg06::patch_piece::minus_z>(pt)};

  return PatchSystem("Thornburg06", std::move(patches), std::move(pt));
}

} // namespace MultiPatch