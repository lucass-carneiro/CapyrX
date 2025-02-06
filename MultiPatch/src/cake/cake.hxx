#ifndef MULTIPATCH_CAKE_HXX
#define MULTIPATCH_CAKE_HXX

#include "multipatch.hxx"
#include "tests.hxx"

#include <string>
#include <type_traits>

namespace MultiPatch {
namespace Cake {

/**
 * @brief Stores a 3-vector
 */
using svec = vec<CCTK_REAL, dim>;

/**
 * @brief Stores a 3-matrix
 */
using smat = smat<CCTK_REAL, dim>;

/**
 * @brief Stores a Jacobian matrix of 3 dimensions
 */
using jac_t = vec<vec<CCTK_REAL, dim>, dim>;

/**
 * @brief Stores a the derivatives of a Jacobian matrix of 3 dimensions
 */
using djac_t = vec<smat, dim>;

/**
 * @brief Tags for each patch piece in the cake
 */
enum class patch_piece : int {
  cartesian = 0,

  plus_x = 1,
  minus_x = 2,

  plus_y = 3,
  minus_y = 4,

  plus_z = 5,
  minus_z = 6,

  unknown = 7
};

/**
 * @brief Gets name from patch piece
 *
 * @param p A patch piece type.
 * @return A string representing the name of the piece.
 */
inline const std::string piece_name(const patch_piece &p) {
  switch (static_cast<int>(p)) {
  case static_cast<int>(patch_piece::cartesian):
    return "cartesian";
  case static_cast<int>(patch_piece::plus_x):
    return "plus x";
  case static_cast<int>(patch_piece::minus_x):
    return "minus x";
  case static_cast<int>(patch_piece::plus_y):
    return "plus y";
  case static_cast<int>(patch_piece::minus_y):
    return "minus y";
  case static_cast<int>(patch_piece::plus_z):
    return "plus z";
  case static_cast<int>(patch_piece::minus_z):
    return "minus z";
  default:
    return "unknown";
  }
}

/**
 * @brief Get the patch piece that owns a global coordinate point.
 *
 * @param pt The PatchTransformations structure describing the patch system.
 * @param global_vars The global coordinate triplet to locate the owner for.
 * @return The patch piece owning the global coordinates.
 */
CCTK_DEVICE CCTK_HOST patch_piece
get_owner_patch(const PatchTransformations &pt, const svec &global_vars);

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
CCTK_DEVICE CCTK_HOST svec local2global_impl(const PatchTransformations &pt,
                                             int patch, const svec &local_vars);

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
global2local_impl(const PatchTransformations &pt, const svec &global_vars);

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
d2local_dglobal2_impl(const PatchTransformations &pt, int patch,
                      const svec &local_vars);

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
dlocal_dglobal_impl(const PatchTransformations &pt, int patch,
                    const svec &local_vars);

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_HXX
