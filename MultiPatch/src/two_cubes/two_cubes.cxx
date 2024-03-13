#include "multipatch.hxx"

namespace MultiPatch {
namespace TwoCubes {

using svec = vec<CCTK_REAL, dim>;
using smat = smat<CCTK_REAL, dim>;
using jac_t = vec<vec<CCTK_REAL, dim>, dim>;
using djac_t = vec<smat, dim>;

enum class patch_piece : int {
  left_cube = 0,
  right_cube = 1,
  exterior = 2,
  unknown = 3
};

inline const char *piece_name(const patch_piece &p) {
  switch (static_cast<int>(p)) {
  case static_cast<int>(patch_piece::left_cube):
    return "left cube";
  case static_cast<int>(patch_piece::right_cube):
    return "right cube";
  case static_cast<int>(patch_piece::exterior):
    return "exterior";
  default:
    return "unknown";
  }
}

CCTK_DEVICE CCTK_HOST svec local2global(const PatchTransformations &pt,
                                        int patch, const svec &local_vars) {

  svec global_vars = {0.0, 0.0, 0.0};

  const auto a{local_vars(0)};
  const auto b{local_vars(1)};
  const auto c{local_vars(2)};

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};

  const auto ymin{pt.two_cubes_ymin};
  const auto ymax{pt.two_cubes_ymax};

  const auto zmin{pt.two_cubes_zmin};
  const auto zmax{pt.two_cubes_zmax};

  switch (patch) {
  case static_cast<int>(patch_piece::left_cube):
    global_vars = {0.25 * (1.0 + a) * (xmax - xmin) + xmin,
                   0.5 * ((1.0 + b) * ymax + (1.0 - b) * ymin),
                   0.5 * ((1.0 + c) * zmax + (1.0 - c) * zmin)};
    break;

  case static_cast<int>(patch_piece::right_cube):
    global_vars = {0.25 * ((3.0 + a) * xmax + (1.0 - a) * xmin),
                   0.5 * ((1.0 + b) * ymax + (1.0 - b) * ymin),
                   0.5 * ((1.0 + c) * zmax + (1.0 - c) * zmin)};
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No local -> global transformations available for patch %s",
                piece_name(static_cast<patch_piece>(patch)));
#else
    assert(0);
#endif
    break;
  }

  return global_vars;
}

CCTK_DEVICE CCTK_HOST std_tuple<int, svec>
global2local(const PatchTransformations &pt, const svec &global_vars) {
  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};

  const auto ymin{pt.two_cubes_ymin};
  const auto ymax{pt.two_cubes_ymax};

  const auto zmin{pt.two_cubes_zmin};
  const auto zmax{pt.two_cubes_zmax};

  // TODO: Implement this
  const auto piece = get_owner_patch(pt, global_vars);

  CCTK_REAL a{0}, b{0}, c{0};

  switch (static_cast<int>(piece)) {
  case static_cast<int>(patch_piece::left_cube):
    a = 3.0 + (4 * (x - xmax)) / (xmax - xmin);
    b = -(-2 * y + ymax + ymin) / (ymax - ymin);
    c = -(-2 * z + zmax + zmin) / (zmax - zmin);
    break;

  case static_cast<int>(patch_piece::right_cube):
    a = 1.0 + (4 * (x - xmax)) / (xmax - xmin);
    b = -(-2 * y + ymax + ymin) / (ymax - ymin);
    c = -(-2 * z + zmax + zmin) / (zmax - zmin);
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

  return std_make_tuple(static_cast<int>(piece), svec{a, b, c});
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
} // namespace TwoCubes

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

} // namespace TwoCubes

/**
 * @brief Creates the TwoCubes patch system
 *
 * @return PatchSystem object with TwoCubes data and functions
 */
PatchSystem SetupCake() {
  PatchTransformations pt;
  pt.global2local = &TwoCubes::global2local;
  pt.local2global = &TwoCubes::local2global;
  pt.dlocal_dglobal = &TwoCubes::dlocal_dglobal;
  pt.d2local_dglobal2 = &TwoCubes::d2local_dglobal2;
  pt.global2local_device = &TwoCubes::global2local;
  pt.local2global_device = &TwoCubes::local2global;
  pt.dlocal_dglobal_device = &TwoCubes::dlocal_dglobal;
  pt.d2local_dglobal2_device = &TwoCubes::d2local_dglobal2;

  const auto patches = std::vector<Patch>{
      TwoCubes::make_patch<TwoCubes::patch_piece::cartesian>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::minus_x>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::plus_x>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::minus_y>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::plus_y>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::minus_z>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::plus_z>(pt)};

  return PatchSystem("TwoCubes", std::move(patches), std::move(pt));
}

} // namespace MultiPatch
