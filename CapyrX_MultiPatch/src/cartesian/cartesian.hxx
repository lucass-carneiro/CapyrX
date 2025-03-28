#ifndef CAPYRX_PATCH_CARTESIAN_HXX
#define CAPYRX_PATCH_CARTESIAN_HXX

#include "multipatch.hxx"

#include <loop_device.hxx>
#include <tuple.hxx>

namespace CapyrX::MultiPatch::Cartesian {

/**
 * Parameters that define a Cartesian patch
 */
struct PatchParams {
  CCTK_INT ncells_i{10};
  CCTK_INT ncells_j{10};
  CCTK_INT ncells_k{10};

  CCTK_REAL xmin{1};
  CCTK_REAL ymin{1};
  CCTK_REAL zmin{1};

  CCTK_REAL xmax{11};
  CCTK_REAL ymax{11};
  CCTK_REAL zmax{11};
};

/**
 * Transforms global coordinates to patch local coordinates.
 *
 * Inputs:
 * 1. `const PatchParams &p`: Patch parameters.
 * 2. `const svec_t &global_vars`: Global coordinates.
 *
 * Returns:
 *  A 2-element tuple containing:
 *  1. The index of the patch owning the local coordinates.
 *  2. The local coordinates.
 */
CCTK_HOST CCTK_DEVICE auto
global2local(const PatchParams &p,
             const svec_t &global_coords) -> std_tuple<int, svec_t>;

/**
 * Transforms patch local coordinates to global coordinates
 *
 * Inputs:
 * 1. `const PatchTransformations &pt`: Reference to a PatchTransformations
 * structure from which patch data will be obtained.
 * 2. `int patch`: The index of the patch that owns the local coordinates to
 * be transformed.
 * 3. `const svec_t &local_vars`: Reference to a tuple of local coordinates
 * to be transformed to global coordiantes.
 *
 * Returns:
 * An svec_t object containing the transformd global coordinates.
 */
CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &p, int patch,
                                        const svec_t &local_coords) -> svec_t;

/**
 * Coordinate transformation's Jacobian components.
 *
 * Inputs:
 * 1. `const PatchParams &p`: Reference to a PatchTransformations
 * structure from which patch data will be obtained.
 * 2. `int patch`: The index representing the patch for which Jacobians should
 * be calculated.
 * 3. `const svec_t &local_vars`: Reference to a tuple of local coordinates
 * where the Jacobian will be evaluated.
 *
 * Returns:
 * A 2 element tuple containing:
 * 1. The global coordinates obtained by calling `local2global` on the
 * `local_vars` input.
 * 2. The components of the local -> global coordinate transformation
 * Jacobian.
 *
 * @note Mathematically, Jacobians are defined as
 * \f[J^{i}_{\phantom{i}j} = \frac{\mathrm{d} a^i}{\mathrm{d} x^j}\f]
 * In code, \f$J^{i}_{\phantom{i}j}\f$ is acessed by calling `J(i)(j)` on a
 * `jac_t` type object.
 */
CCTK_HOST CCTK_DEVICE auto
dlocal_dglobal(const PatchParams &p, int patch,
               const svec_t &local_coords) -> std_tuple<svec_t, jac_t>;

/**
 * Computes Jacobian derivative components.
 *
 * @par Inputs:
 * 1. `const PatchParams &p`: Reference to a PatchTransformations
 * structure from which patch data will be obtained.
 * 2. `int patch`: The index representing the patch for which Jacobian
 * derivatives should be calculated.
 * 3. `const svec_t &local_vars`: Reference to a tuple of local coordinates
 * where the Jacobian will be evaluated.
 *
 * @par Returns:
 * A 3 element tuple containing:
 * 1. The global coordinates obtained by calling `local2global` on the
 * `local_vars` input.
 * 2. The components of the local -> global coordinate transformation
 * Jacobian, obtained by calling `dlocal_dglobal` on the input data.
 * 3. The components of the local -> global coordinate transformation Jacobian
 * derivatives.
 *
 * @note Mathematically, Jacobians are defined as
 * \f[J^{i}_{\phantom{i}j} = \frac{\mathrm{d} a^i}{\mathrm{d} x^j}\f]
 * and their derivatives as
 * \f[J^{i}_{\phantom{i}jk} = \frac{\mathrm{d}^2 a^i}{\mathrm{d} x^j
 * \mathrm{d}x^k}\f] In code, \f$J^{i}_{\phantom{i}jk}\f$ is acessed by
 * calling `J(i)(j,k)` on a `djac_t` type object.
 */
CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &p, int patch, const svec_t &local_coords) -> std_tuple<svec_t, jac_t, djac_t>;

/**
 * Creates a cartesian PatchSystem
 */
auto make_system(const PatchParams &p) -> PatchSystem;

} // namespace CapyrX::MultiPatch::Cartesian

#endif // CAPYRX_PATCH_CARTESIAN_HXX
