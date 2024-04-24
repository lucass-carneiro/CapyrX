#include "cake.hxx"

#include <array>
#include <algorithm>
#include <cassert>

/**
 * @brief Get the patch piece that owns a global coordinate point.
 *
 * @param pt The PatchTransformations structure describing the patch system.
 * @param global_vars The global coordinate triplet to locate the owner for.
 * @return The patch piece owning the global coordinates.
 */
CCTK_DEVICE CCTK_HOST MultiPatch::Cake::patch_piece
MultiPatch::Cake::get_owner_patch(const PatchTransformations &pt,
                                  const svec &global_vars) {
  using std::abs;
  using std::distance;
  using std::max_element;

  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  const auto r0 = pt.cake_inner_boundary_radius;

  const bool x_within_r0{-r0 < x && x < r0};
  const bool y_within_r0{-r0 < y && y < r0};
  const bool z_within_r0{-r0 < z && z < r0};

  // Are we on the Cartesian core?
  if (x_within_r0 && y_within_r0 && z_within_r0) {
    return patch_piece::cartesian;
  }

  /* We are on one of the 6 Thornburg patches. The patch direction is the
   * determined by the coordinate with the largest absolute value
   */
  const std::array<CCTK_REAL, 3> abs_global_coords{abs(x), abs(y), abs(z)};
  const auto max_abs_global_coord_idx{distance(
      abs_global_coords.begin(),
      max_element(abs_global_coords.begin(), abs_global_coords.end()))};

  /* After the direction is know, we must dtermine the sense, which is
   * represented by the sign of the coordinate with the largest absolute
   * value.
   */
  const int max_coord_sign{global_vars(max_abs_global_coord_idx) > 0.0 ? 1
                                                                       : -1};

  if (max_coord_sign < 0 && max_abs_global_coord_idx == 0) {
    return patch_piece::minus_x;
  } else if (max_coord_sign > 0 && max_abs_global_coord_idx == 0) {
    return patch_piece::plus_x;
  } else if (max_coord_sign < 0 && max_abs_global_coord_idx == 1) {
    return patch_piece::minus_y;
  } else if (max_coord_sign > 0 && max_abs_global_coord_idx == 1) {
    return patch_piece::plus_y;
  } else if (max_coord_sign < 0 && max_abs_global_coord_idx == 2) {
    return patch_piece::minus_z;
  } else if (max_coord_sign > 0 && max_abs_global_coord_idx == 2) {
    return patch_piece::plus_z;
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
