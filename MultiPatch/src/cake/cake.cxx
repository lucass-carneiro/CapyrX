#include "multipatch.hxx"
#include "tests.hxx"
#include "cake.hxx"
#include "cake_jacobians.hxx"

#include <cassert>
#include <cmath>
#include <string>

namespace MultiPatch {
namespace Cake {

/**
 * @brief The local to global coordinate transformation implementation.
 * This is where the actual coordinate transformations are performed but it is
 * also suitable for passing to CarpetX.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A spatial vector containing the coordinate transformations.
 */
CCTK_DEVICE CCTK_HOST svec local2global(const PatchTransformations &pt,
                                        int patch, const svec &local_vars) {
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
    global_vars = local_vars;
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

/**
 * @brief The global to local coordinate transformation implementation.
 * This is where the actual coordinate transformations are performed but it is
 * also suitable to be passed to CarpetX without any wrappers.
 *
 * @param pt The patch data
 * @param global_vars The values of the local global (x, y, z)
 * @return A tuple containing the patch piece and the global coordinate triplet.
 */
CCTK_DEVICE CCTK_HOST std_tuple<int, svec>
global2local(const PatchTransformations &pt, const svec &global_vars) {
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
    local_vars = global_vars;
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

/**
 * @brief This function computes the local to global coordinate
 * transformation, the jacobian and it's derivative. It can be passed directly
 * to CarpetX.
 *
 * @note The Jacobians are defined as
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 * dJ(i)(j,k) = $dJ^{i}_{j k} = \frac{d^2 a^i}{d x^j d x^k}
 * \right)$.
 *
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation,
 * the local to global jacobian matrix and it's derivative.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t, djac_t>
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const svec &local_vars) {

  const auto local_to_global_result = pt.local2global(pt, patch, local_vars);
  const auto jacobian_results = cake_jacs(pt, patch, local_to_global_result);

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
} // namespace Cake

/**
 * @brief This function computes the local to global coordinate transformation
 * and the jacobian. It can be passed directly to CarpetX.
 *
 * @note NThe Jacobians is defined as:
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation
 * and the local to global jacobian matrix.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t>
dlocal_dglobal(const PatchTransformations &pt, int patch,
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

  PatchFace co{false, static_cast<int>(patch_piece::cartesian)};
  PatchFace ex{true, static_cast<int>(patch_piece::exterior)};
  PatchFace px{false, static_cast<int>(patch_piece::plus_x)};
  PatchFace mx{false, static_cast<int>(patch_piece::minus_x)};
  PatchFace py{false, static_cast<int>(patch_piece::plus_y)};
  PatchFace my{false, static_cast<int>(patch_piece::minus_y)};
  PatchFace pz{false, static_cast<int>(patch_piece::plus_z)};
  PatchFace mz{false, static_cast<int>(patch_piece::minus_z)};

  if constexpr (p == patch_piece::cartesian) {
    patch.ncells = {pt.cake_cartesian_ncells_i, pt.cake_cartesian_ncells_j,
                    pt.cake_cartesian_ncells_k};

    patch.xmin = {-pt.cake_inner_boundary_radius,
                  -pt.cake_inner_boundary_radius,
                  -pt.cake_inner_boundary_radius};

    patch.xmax = {pt.cake_inner_boundary_radius, pt.cake_inner_boundary_radius,
                  pt.cake_inner_boundary_radius};

    patch.is_cartesian = true;

    patch.faces = {{mx, my, mz}, {px, py, pz}};

  } else if constexpr (p == patch_piece::plus_x) {
    patch.faces = {{mz, my, co}, {pz, py, ex}};

  } else if constexpr (p == patch_piece::minus_x) {
    patch.faces = {{mz, py, co}, {pz, my, ex}};

  } else if constexpr (p == patch_piece::plus_y) {
    patch.faces = {{mz, px, co}, {pz, mx, ex}};

  } else if constexpr (p == patch_piece::minus_y) {
    patch.faces = {{mz, mx, co}, {pz, px, ex}};

  } else if constexpr (p == patch_piece::plus_z) {
    patch.faces = {{px, my, co}, {mx, py, ex}};

  } else if constexpr (p == patch_piece::minus_z) {
    patch.faces = {{mx, my, co}, {px, py, ex}};
  }

  return patch;
}

} // namespace Cake

/**
 * @brief Creates the Cake patch system
 *
 * @return PatchSystem object with Cake data and functions
 */
PatchSystem SetupCake() {
  PatchTransformations pt;
  pt.global2local = &Cake::global2local;
  pt.local2global = &Cake::local2global;
  pt.dlocal_dglobal = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2 = &Cake::d2local_dglobal2;
  pt.global2local_device = &Cake::global2local;
  pt.local2global_device = &Cake::local2global;
  pt.dlocal_dglobal_device = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2_device = &Cake::d2local_dglobal2;

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
