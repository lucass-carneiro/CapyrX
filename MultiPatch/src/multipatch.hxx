#ifndef CAPYRX_MULTIPATCH_MULTIPATCH_HXX
#define CAPYRX_MULTIPATCH_MULTIPATCH_HXX

#include <cctk.h>

#include <loop.hxx>
#include <mat.hxx>
#include <tuple.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <memory>
#include <string>
#include <utility>
#include <vector>

/**
 * The namespace where all MultiPatch related infrastructure is located
 */
namespace MultiPatch {

using namespace Arith;

/**
 * The number of spatial dimentions each patch will cover
 */
constexpr int dim = 3;

////////////////////////////////////////////////////////////////////////////////

/**
 * Describes the connections of one patch with another patch.
 *
 * The interface between two patches is called a patch face. This struct
 * describes such interfaces wherever two patches meet.
 */
struct PatchFace {
  /**
   * True if the face connects with the exterior of the simulation domain.
   */
  bool is_outer_boundary;

  /**
   * The index of the neighbouring patch, or -1 if the face is the outer
   *boundary.
   **/
  int other_patch;
};

/**
 * Describes a patch in the patch system
 */
struct Patch {
  /**
   * String representation of the name of the patch.
   */
  std::string name;

  /**
   * The number of cells in the logical x, y and z dimentions that the patch wil
   * contain.
   */
  vect<int, dim> ncells;

  /**
   * Lower coordinate boundary for each logical x, y, z dimention in the patch
   */
  vect<CCTK_REAL, dim> xmin;

  /**
   * Upper coordinate boundary for each logical x, y, z dimention in the patch
   */
  vect<CCTK_REAL, dim> xmax;

  /**
   * Wether or not a patch is cartesian. If true, this flag is used to save time
   * during Jacobian calculations.
   */
  bool is_cartesian;

  /**
   * Store the 6 faces of a patch, repreenting its connection with other patches
   */
  vect<vect<PatchFace, dim>, 2> faces;
};

////////////////////////////////////////////////////////////////////////////////

/**
 * Stores patch transformation functions and patch data.
 *
 * %Patch transformation functions refer to the local -> global and global ->
 * local coordinate transformations, as well as coordinate transformation
 * Jacobian and their respective derivatives. These functions are patch system
 * specific, and so they are stored as funcion pointers within the struct. All
 * function pointers are duplicate, since one set of functions is intended to be
 * executed on hosts (CPUs) and another set is intended to be executed on
 * devices (GPUs).
 *
 * %Patch data refers to each patch specific data that is required to setup a
 * patch, such as the number of cells in the patch, global extensions and other
 * parameters that may be required.
 */
struct PatchTransformations {
  /**
   * Construct a new Patch Transformations object.
   *
   * The constructor reads cactus parameters and assigns them to the internal
   * data storage during construction.
   */
  PatchTransformations();

  /**
   * Construct a new Patch Transformations object by copying another object of
   * the same type
   */
  PatchTransformations(const PatchTransformations &) = default;

  /**
   * Construct a new Patch Transformations object by moving from another object
   * of the same type
   */
  PatchTransformations(PatchTransformations &&) = default;

  /**
   * Construct a new Patch Transformations object using the assignment operator
   * by copying another object of the same type
   */
  PatchTransformations &operator=(const PatchTransformations &) = default;

  /**
   * Construct a new Patch Transformations object using the assignment operator
   * by moving from another object of the same type
   */
  PatchTransformations &operator=(PatchTransformations &&) = default;

  /**
   * Type alias for representing spatial vectors. Used for representing
   * coordinate tuples.
   */
  using svec_t = vec<CCTK_REAL, dim>;

  /**
   * Type alias for representing Jacobian components
   */
  using jac_t = vec<vec<CCTK_REAL, dim>, dim>;

  /**
   * Type alias for representing Jacobian derivative components
   */
  using djac_t = vec<smat<CCTK_REAL, dim>, dim>;

  /**
   * Type alias for a function pointer that points to a `global2local` like
   * function, responsible for transforming coordinates from global coordinates
   * to patch local coordinates.
   *
   * @par Inputs:
   * 1. `const PatchTransformations &pt`: Reference to a PatchTransformations
   * structure from which patch data will be obtained.
   * 2. `const svec_t &global_vars`: Reference to a tuple of global coordinates
   * to be transformed to local coordiantes.
   *
   * @par Returns:
   *  A 2-element tuple containing:
   *  1. The index of the patch owning the local coordinates.\n
   *  2. The local coordinates.
   */
  using global2local_func = std_tuple<int, svec_t> (*)(
      const PatchTransformations &pt, const svec_t &global_vars);

  /**
   * Type alias for a function pointer that points to a `local2global` like
   * function, responsible for transforming coordinates from patch local
   * coordinates to global coordinates.
   *
   * @par Inputs:
   * 1. `const PatchTransformations &pt`: Reference to a PatchTransformations
   * structure from which patch data will be obtained.
   * 2. `int patch`: The index of the patch that owns the local coordinates to
   * be transformed.
   * 3. `const svec_t &local_vars`: Reference to a tuple of local coordinates
   * to be transformed to global coordiantes.
   *
   * @par Returns:
   * An svec_t object containing the transformd global coordinates.
   */
  using local2global_func = svec_t (*)(const PatchTransformations &pt,
                                       int patch, const svec_t &local_vars);

  /**
   * Type alias for a function pointer that points to a `dlocal_dglobal` like
   * function, responsible for computing the coordinate transformation's
   * Jacobian components.
   *
   * @par Inputs:
   * 1. `const PatchTransformations &pt`: Reference to a PatchTransformations
   * structure from which patch data will be obtained.
   * 2. `int patch`: The index representing the patch for which Jacobians should
   * be calculated.
   * 3. `const svec_t &local_vars`: Reference to a tuple of local coordinates
   * where the Jacobian will be evaluated.
   *
   * @par Returns:
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
  using dlocal_dglobal_fun = std_tuple<svec_t, jac_t> (*)(
      const PatchTransformations &pt, int patch, const svec_t &local_vars);

  /**
   * Type alias for a function pointer that points to a `d2local_dglobal2_fun`
   * like function, responsible for computing the coordinate transformation's
   * Jacobian derivative components.
   *
   * @par Inputs:
   * 1. `const PatchTransformations &pt`: Reference to a PatchTransformations
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
  using d2local_dglobal2_fun = std_tuple<svec_t, jac_t, djac_t> (*)(
      const PatchTransformations &pt, int patch, const svec_t &local_vars);

  /**
   * Pointer to a function with signature as described in global2local_func,
   * performing global -> local coordinate transformations. This function is
   * called on the host (CPU).
   */
  global2local_func global2local{nullptr};

  /**
   * Pointer to a function with signature as described in local2global_func,
   * performing local -> global coordinate transformations. This function is
   * called on the host (CPU).
   */
  local2global_func local2global{nullptr};

  /**
   * Pointer to a function with signature as described in dlocal_dglobal_fun,
   * computing local -> global coordinate transformations Jacobian components.
   * This function is called on the host (CPU).
   */
  dlocal_dglobal_fun dlocal_dglobal{nullptr};

  /**
   * Pointer to a function with signature as described in d2local_dglobal2_fun,
   * computing local -> global coordinate transformations. Jacobian derivative
   * components. This function is called on the host (CPU).
   */
  d2local_dglobal2_fun d2local_dglobal2{nullptr};

  /**
   * Pointer to a function with signature as described in global2local_func,
   * performing global -> local coordinate transformations. This function is
   * called on the device (GPU).
   */
  global2local_func global2local_device{nullptr};

  /**
   * Pointer to a function with signature as described in local2global_func,
   * performing local -> global coordinate transformations. This function is
   * called on the device (GPU).
   */
  local2global_func local2global_device{nullptr};

  /**
   * Pointer to a function with signature as described in dlocal_dglobal_fun,
   * computing local -> global coordinate transformations Jacobian components.
   * This function is called on the device (GPU).
   */
  dlocal_dglobal_fun dlocal_dglobal_device{nullptr};

  /**
   * Pointer to a function with signature as described in d2local_dglobal2_fun,
   * computing local -> global coordinate transformations. Jacobian derivative
   * components. This function is called on the device (GPU).
   */
  d2local_dglobal2_fun d2local_dglobal2_device{nullptr};

  /**
   * Cartesian %Patch: Upper x coordinate boundary
   */
  CCTK_REAL cartesian_xmax;

  /**
   * Cartesian %Patch: Lower x coordinate boundary
   */
  CCTK_REAL cartesian_xmin;

  /**
   * Cartesian %Patch: Upper y coordinate boundary
   */
  CCTK_REAL cartesian_ymax;

  /**
   * Cartesian %Patch: Lower y coordinate boundary
   */
  CCTK_REAL cartesian_ymin;

  /**
   * Cartesian %Patch: Upper z coordinate boundary
   */
  CCTK_REAL cartesian_zmax;

  /**
   * Cartesian %Patch: Lower z coordinate boundary
   */
  CCTK_REAL cartesian_zmin;

  /**
   * Cartesian %Patch: Number of cells in the x direction
   */
  CCTK_INT cartesian_ncells_i;

  /**
   * Cartesian %Patch: Number of cells in the y direction
   */
  CCTK_INT cartesian_ncells_j;

  /**
   * Cartesian %Patch: Number of cells in the z direction
   */
  CCTK_INT cartesian_ncells_k;

  /**
   * Cubed Sphere %Patch: Inner radius
   */
  CCTK_REAL cubed_sphere_rmin;

  /**
   * Cubed Sphere %Patch: Outer radius
   */
  CCTK_REAL cubed_sphere_rmax;

  /**
   * Swirl %Patch: Number of cells in the x direction
   */
  CCTK_INT swirl_ncells_i;

  /**
   * Swirl %Patch: Number of cells in the y direction
   */
  CCTK_INT swirl_ncells_j;

  /**
   * Swirl %Patch: Number of cells in the z direction
   */
  CCTK_INT swirl_ncells_k;

  /**
   * Cake %Patch: Radius of the outer boundary
   */
  CCTK_REAL cake_outer_boundary_radius;

  /**
   * Cake %Patch: Half the coordinate length of the central cartesian cube's
   * face
   */
  CCTK_REAL cake_inner_boundary_radius;

  /**
   * Cake %Patch: The number of cells in the x direction of the central
   * cartesian cube.
   */
  CCTK_INT cake_cartesian_ncells_i;

  /**
   * Cake %Patch: The number of cells in the y direction of the central
   * cartesian cube.
   */
  CCTK_INT cake_cartesian_ncells_j;

  /**
   * Cake %Patch: The number of cells in the z direction of the central
   * cartesian cube.
   */
  CCTK_INT cake_cartesian_ncells_k;

  /**
   * Cake %Patch: The number of cells in the angular direction of spherical
   * patches.
   */
  CCTK_INT cake_angular_cells;

  /**
   * Cake %Patch: The number of cells in the radial direction of spherical
   * patches.
   */
  CCTK_INT cake_radial_cells;

  /**
   * Two Cubes %Patch: The starting x value of the coordinate system.
   */
  CCTK_REAL two_cubes_xmin;

  /**
   * Two Cubes %Patch: The final x value of the coordinate system.
   */
  CCTK_REAL two_cubes_xmax;

  /**
   * Two Cubes %Patch: The starting y value of the coordinate system.
   */
  CCTK_REAL two_cubes_ymin;

  /**
   * Two Cubes %Patch: The final y value of the coordinate system.
   */
  CCTK_REAL two_cubes_ymax;

  /**
   * Two Cubes %Patch: The starting z value of the coordinate system.
   */
  CCTK_REAL two_cubes_zmin;

  /**
   * Two Cubes %Patch: The final z value of the coordinate system.
   */
  CCTK_REAL two_cubes_zmax;

  /**
   * Two Cubes %Patch:The shift of the patch interface in the y axis.
   */
  CCTK_REAL two_cubes_delta_y;

  /**
   * Two Cubes %Patch: Number of cells in the left cube.
   */
  CCTK_INT two_cubes_ncells_left;

  /**
   * Two Cubes %Patch: Number of cells in the right cube.
   */
  CCTK_INT two_cubes_ncells_right;

  /**
   * Two Cubes %Patch: Number of cells in the y direction.
   */
  CCTK_INT two_cubes_ncells_y;

  /**
   * Two Cubes %Patch: Number of cells in the z direction.
   */
  CCTK_INT two_cubes_ncells_z;

  /**
   * Thornburg06 patch outer boundary radius.
   */
  CCTK_REAL thornburg06_outer_boundary_radius;

  /**
   * Thornburg06 patch inner boundary radius.
   */
  CCTK_REAL thornburg06_inner_boundary_radius;

  /**
   * Thornburg06 patch angular cells.
   */
  CCTK_INT thornburg06_angular_cells;

  /**
   * Thornburg06 patch radial cells.
   */
  CCTK_INT thornburg06_radial_cells;
};

////////////////////////////////////////////////////////////////////////////////

/**
 * Describes a patch system, that is, a collection of patches and patch
 * transformations
 */
struct PatchSystem {

  /**
   * String repreentation of the name of the patch system
   */
  std::string name;

  /**
   * A vector containing all patches in the system.
   */
  std::vector<Patch> patches;

  /**
   * Returns the current number of patches in the system.
   *
   * @par Returns:
   * an `int` containing the number of patches in the system.
   */
  int num_patches() const { return patches.size(); }

  /**
   * A PatchTransformations containing patch data and coordinate transfomation
   * functions.
   */
  PatchTransformations transformations;

  /**
   * Constructs a new empty PatchSystem object
   */
  PatchSystem() {}

  /**
   * Construct a new Patch System object
   *
   * @par Inputs:
   * 1. `std::string name`: The name of the patch system.
   * 2. `std::vector<Patch> patches`: The patches comprising the patch system.
   * 3. `PatchTransformations transformations`: The patch system coordinate
   * transformations and data object.
   */
  PatchSystem(std::string name, std::vector<Patch> patches,
              PatchTransformations transformations)
      : name(std::move(name)), patches(std::move(patches)),
        transformations(std::move(transformations)) {}
};

////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Creates a Cartesian patch system
 *
 * @par Returns:
 * A PatchSystem object with Cartesian coordinates.
 */
PatchSystem SetupCartesian();

/**
 * @brief Creates a Cubed Sphere sphere patch system
 *
 * @par Returns:
 * A PatchSystem object with Cubed Sphere coordinates.
 */
PatchSystem SetupCubedSphere();

/**
 * @brief Creates a Swirl patch system
 *
 * @par Returns:
 * A PatchSystem object with Swirl coordinates.
 */
PatchSystem SetupSwirl();

/**
 * @brief Creates a Cake patch system
 *
 * @par Returns:
 * A PatchSystem object with Cake coordinates.
 */
PatchSystem SetupCake();

/**
 * @brief Creates a Two Cubes patch system
 *
 * @par Returns:
 * A PatchSystem object with Two Cubes coordinates.
 */
PatchSystem SetupTwoCubes();

/**
 * @brief Creates a Thornburg06 patch system
 *
 * @par Returns:
 * A PatchSystem object with Thornburg06 coordinates.
 */
PatchSystem SetupThornburg06();

/**
 * Pointer to the global PatchSystem objct use during a simulation
 */
extern std::unique_ptr<PatchSystem> the_patch_system;

} // namespace MultiPatch

#endif // #ifndef CAPYRX_MULTIPATCH_MULTIPATCH_HXX
