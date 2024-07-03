#include "cake.hxx"

#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>

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
