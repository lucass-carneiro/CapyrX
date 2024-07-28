#ifndef MULTIPATCH_THORNBURG06_HXX
#define MULTIPATCH_THORNBURG06_HXX

#include "multipatch.hxx"

#include <string>

namespace MultiPatch::Thornburg06 {

/* Some type aliases to help conserve finger grease */
using svec = vec<CCTK_REAL, dim>;            // A 3-vector
using smat = smat<CCTK_REAL, dim>;           // A 3x3 matrix
using jac_t = vec<vec<CCTK_REAL, dim>, dim>; // A Jacobian
using djac_t = vec<smat, dim>;               // A Jacobian derivative

/* An enum with patch names, so we can refer to names instead of indices */
enum patch_piece : int {
  plus_x = 0,
  plus_y = 1,
  minus_x = 2,
  minus_y = 3,
  plus_z = 4,
  minus_z = 5,
  unknown = 6
};

/**
 * @brief Gets name from a patch piece
 *
 * @param p A patch piece type.
 * @return A string representing the name of the piece.
 */
std::string piece_name(const patch_piece &p);

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
CCTK_DEVICE CCTK_HOST svec local2global(const PatchTransformations &pt,
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
global2local(const PatchTransformations &pt, const svec &global_vars);

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
dlocal_dglobal(const PatchTransformations &pt, int patch,
               const svec &local_vars);

/**
 * @brief Creates patch piece
 *
 * @tparam p The piece of the patch to make.
 * @param pt The patch transformation object with patch data.
 * @return The constructed patch piece.
 */
template <patch_piece p> Patch make_patch(const PatchTransformations &pt);

} // namespace MultiPatch::Thornburg06

#endif // MULTIPATCH_THORNBURG06_HXX