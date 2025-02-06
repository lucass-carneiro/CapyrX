#include "multipatch.hxx"
#include "tests.hxx"
#include "cake.hxx"
#include "cake_jacobians.hxx"

#include <cassert>
#include <cmath>
#include <string>

namespace MultiPatch {
namespace Cake {

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

  const auto r0 = pt.cake_inner_boundary_radius;

  if (abs_x <= r0 && abs_y <= r0 && abs_z <= r0) {
    return patch_piece::cartesian;
  }

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

CCTK_DEVICE CCTK_HOST svec local2global_impl(const PatchTransformations &pt,
                                             int patch,
                                             const svec &local_vars) {
  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const auto a{local_vars(0)};
  const auto b{local_vars(1)};
  const auto c{local_vars(2)};

  svec global_vars = {0.0, 0.0, 0.0};

  switch (patch) {

  case static_cast<int>(patch_piece::cartesian):
    global_vars(0) = a * r0;
    global_vars(1) = b * r0;
    global_vars(2) = c * r0;
    break;

  case static_cast<int>(patch_piece::plus_x):
    global_vars(0) =
        (r0 - c * r0 + r1 + c * r1) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(1) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case static_cast<int>(patch_piece::minus_x):
    global_vars(0) =
        -((r0 - c * r0 + r1 + c * r1) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_vars(1) =
        -((b * (r0 - c * r0 + r1 + c * r1)) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_vars(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case static_cast<int>(patch_piece::plus_y):
    global_vars(0) =
        -((b * (r0 - c * r0 + r1 + c * r1)) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_vars(1) =
        (r0 - c * r0 + r1 + c * r1) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case static_cast<int>(patch_piece::minus_y):
    global_vars(0) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(1) =
        -((r0 - c * r0 + r1 + c * r1) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_vars(2) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case static_cast<int>(patch_piece::plus_z):
    global_vars(0) =
        -((a * (r0 - c * r0 + r1 + c * r1)) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
    global_vars(1) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(2) =
        (r0 - c * r0 + r1 + c * r1) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    break;

  case static_cast<int>(patch_piece::minus_z):
    global_vars(0) =
        (a * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(1) =
        (b * (r0 - c * r0 + r1 + c * r1)) /
        (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c)));
    global_vars(2) =
        -((r0 - c * r0 + r1 + c * r1) /
          (sqrt(2) * sqrt(2 + pow(a, 2) * (1 + c) + pow(b, 2) * (1 + c))));
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
global2local_impl(const PatchTransformations &pt, const svec &global_vars) {
  using std::pow;
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  const auto piece{get_owner_patch(pt, global_vars)};

  svec local_vars{0.0, 0.0, 0.0};

  switch (static_cast<int>(piece)) {

  case static_cast<int>(patch_piece::cartesian):
    local_vars(0) = x / r0;
    local_vars(1) = y / r0;
    local_vars(2) = z / r0;
    break;

  case static_cast<int>(patch_piece::plus_x):
    local_vars(0) = z / x;
    local_vars(1) = y / x;
    local_vars(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(y, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::minus_x):
    local_vars(0) = -(z / x);
    local_vars(1) = y / x;
    local_vars(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(y, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(x, 2) + pow(pow(y, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (2 * pow(x, 2) + pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::plus_y):
    local_vars(0) = z / y;
    local_vars(1) = -(x / y);
    local_vars(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::minus_y):
    local_vars(0) = -(z / y);
    local_vars(1) = -(x / y);
    local_vars(2) =
        (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(z, 2) +
         sqrt(4 * pow(r1, 2) * pow(y, 2) + pow(pow(x, 2) + pow(z, 2), 2) +
              4 * pow(r0, 2) * (pow(x, 2) + pow(y, 2) + pow(z, 2)) -
              4 * r0 * r1 * (pow(x, 2) + 2 * pow(y, 2) + pow(z, 2)))) /
        pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::plus_z):
    local_vars(0) = -(x / z);
    local_vars(1) = y / z;
    local_vars(2) = (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(y, 2) +
                     sqrt((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2))) /
                    pow(r0 - r1, 2);
    break;

  case static_cast<int>(patch_piece::minus_z):
    local_vars(0) = -(x / z);
    local_vars(1) = -(y / z);
    local_vars(2) = (pow(r0, 2) - pow(r1, 2) + pow(x, 2) + pow(y, 2) +
                     sqrt((pow(x, 2) + pow(y, 2)) *
                              (4 * r0 * (r0 - r1) + pow(x, 2) + pow(y, 2)) +
                          4 * pow(r0 - r1, 2) * pow(z, 2))) /
                    pow(r0 - r1, 2);
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("At point (%f, %f, %f): No global -> local transformations "
                "available for patch %s.",
                x, y, z, piece_name(piece).c_str());
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(static_cast<int>(piece), local_vars);
}

CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t, djac_t>
d2local_dglobal2_impl(const PatchTransformations &pt, int patch,
                      const svec &local_vars) {

  const auto local_to_global_result = pt.local2global(pt, patch, local_vars);
  const auto jacobian_results = cake_jacs(pt, patch, local_to_global_result);

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
}

CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t>
dlocal_dglobal_impl(const PatchTransformations &pt, int patch,
                    const svec &local_vars) {

  const auto data = pt.d2local_dglobal2(pt, patch, local_vars);
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

/**
 * @brief Creates a cake patch piece
 *
 * @tparam p The piece of the patch to make.
 * @param pt The patch transformation object with patch data.
 * @return The constructed patch piece.
 */
template <patch_piece p> Patch make_patch(const PatchTransformations &pt) {
  Patch patch;

  patch.name = piece_name(p);

  // Basic configuration for a thornburg patch piece
  patch.ncells = {pt.cake_angular_cells, pt.cake_angular_cells,
                  pt.cake_radial_cells};

  patch.xmin = {-1.0, -1.0, -1.0};
  patch.xmax = {1.0, 1.0, 1.0};

  patch.is_cartesian = false;

  PatchFace ob{true, -1};
  PatchFace co{false, static_cast<int>(patch_piece::cartesian)};
  PatchFace px{false, static_cast<int>(patch_piece::plus_x)};
  PatchFace mx{false, static_cast<int>(patch_piece::minus_x)};
  PatchFace py{false, static_cast<int>(patch_piece::plus_y)};
  PatchFace my{false, static_cast<int>(patch_piece::minus_y)};
  PatchFace pz{false, static_cast<int>(patch_piece::plus_z)};
  PatchFace mz{false, static_cast<int>(patch_piece::minus_z)};

  if constexpr (p == patch_piece::cartesian) {
    patch.ncells = {pt.cake_cartesian_ncells_i, pt.cake_cartesian_ncells_j,
                    pt.cake_cartesian_ncells_k};

    patch.is_cartesian = true;

    patch.faces = {{mx, my, mz}, {px, py, pz}};

  } else if constexpr (p == patch_piece::plus_x) {
    patch.faces = {{mz, my, co}, {pz, py, ob}};

  } else if constexpr (p == patch_piece::minus_x) {
    patch.faces = {{mz, py, co}, {pz, my, ob}};

  } else if constexpr (p == patch_piece::plus_y) {
    patch.faces = {{mz, px, co}, {pz, mx, ob}};

  } else if constexpr (p == patch_piece::minus_y) {
    patch.faces = {{mz, mx, co}, {pz, px, ob}};

  } else if constexpr (p == patch_piece::plus_z) {
    patch.faces = {{px, my, co}, {mx, py, ob}};

  } else if constexpr (p == patch_piece::minus_z) {
    patch.faces = {{mx, my, co}, {px, py, ob}};
  }

  return patch;
}

// Host functions
svec local2global(const PatchTransformations &pt, int patch,
                  const svec &local_vars) {
  return local2global_impl(pt, patch, local_vars);
}

std_tuple<int, svec> global2local(const PatchTransformations &pt,
                                  const svec &global_vars) {
  return global2local_impl(pt, global_vars);
}

std_tuple<svec, jac_t> dlocal_dglobal(const PatchTransformations &pt, int patch,
                                      const svec &local_vars) {
  return dlocal_dglobal_impl(pt, patch, local_vars);
}

std_tuple<svec, jac_t, djac_t> d2local_dglobal2(const PatchTransformations &pt,
                                                int patch,
                                                const svec &local_vars) {
  return d2local_dglobal2_impl(pt, patch, local_vars);
}

// Device functions
CCTK_DEVICE svec local2global_device(const PatchTransformations &pt, int patch,
                                     const svec &local_vars) {
  return local2global_impl(pt, patch, local_vars);
}

CCTK_DEVICE std_tuple<int, svec>
global2local_device(const PatchTransformations &pt, const svec &global_vars) {
  return global2local_impl(pt, global_vars);
}

CCTK_DEVICE std_tuple<svec, jac_t>
dlocal_dglobal_device(const PatchTransformations &pt, int patch,
                      const svec &local_vars) {
  return dlocal_dglobal_impl(pt, patch, local_vars);
}

CCTK_DEVICE std_tuple<svec, jac_t, djac_t>
d2local_dglobal2_device(const PatchTransformations &pt, int patch,
                        const svec &local_vars) {
  return d2local_dglobal2_impl(pt, patch, local_vars);
}

} // namespace Cake

PatchSystem SetupCake() {
  PatchTransformations pt;
  pt.global2local = &Cake::global2local;
  pt.local2global = &Cake::local2global;
  pt.dlocal_dglobal = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2 = &Cake::d2local_dglobal2;
  pt.global2local_device = &Cake::global2local_device;
  pt.local2global_device = &Cake::local2global_device;
  pt.dlocal_dglobal_device = &Cake::dlocal_dglobal_device;
  pt.d2local_dglobal2_device = &Cake::d2local_dglobal2_device;

  const auto patches =
      std::vector<Patch>{Cake::make_patch<Cake::patch_piece::cartesian>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_x>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_x>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_y>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_y>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_z>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_z>(pt)};

  return PatchSystem("Cake", std::move(patches), std::move(pt));
}

} // namespace MultiPatch
