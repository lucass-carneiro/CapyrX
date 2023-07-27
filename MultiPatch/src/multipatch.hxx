#ifndef MULTIPATCH_HXX
#define MULTIPATCH_HXX

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

namespace MultiPatch {
using namespace Arith;

/**
 * The number of spatial dimentions each patch will cover
 */
constexpr int dim = 3;

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
   * String repreentation of the name of the patch.
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
   * TODO: Request Erik's help
   */
  vect<vect<PatchFace, dim>, 2> faces;
};

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
   * Cartesian %Patch: Lower x coordinate boundary
   */
  CCTK_REAL cartesian_xmin;

  /**
   * Cartesian %Patch: Upper x coordinate boundary
   */
  CCTK_REAL cartesian_xmax;

  /**
   * Cartesian %Patch: Lower y coordinate boundary
   */
  CCTK_REAL cartesian_ymin;

  /**
   * Cartesian %Patch: Upper y coordinate boundary
   */
  CCTK_REAL cartesian_ymax;

  /**
   * Cartesian %Patch: Lower z coordinate boundary
   */
  CCTK_REAL cartesian_zmin;

  /**
   * Cartesian %Patch: Upper z coordinate boundary
   */
  CCTK_REAL cartesian_zmax;

  /**
   * Cartesian %Patch: Number of cells in the x direction
   */
  int cartesian_ncells_i;

  /**
   * Cartesian %Patch: Number of cells in the y direction
   */
  int cartesian_ncells_j;

  /**
   * Cartesian %Patch: Number of cells in the z direction
   */
  int cartesian_ncells_k;

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
  int swirl_ncells_i;

  /**
   * Swirl %Patch: Number of cells in the y direction
   */
  int swirl_ncells_j;

  /**
   * Swirl %Patch: Number of cells in the z direction
   */
  int swirl_ncells_k;

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
  int cake_cartesian_ncells_i;

  /**
   * Cake %Patch: The number of cells in the y direction of the central
   * cartesian cube.
   */
  int cake_cartesian_ncells_j;

  /**
   * Cake %Patch: The number of cells in the z direction of the central
   * cartesian cube.
   */
  int cake_cartesian_ncells_k;

  /**
   * Cake %Patch: The number of cells in the angular direction of spherical
   * patches.
   */
  int cake_angular_cells;

  /**
   * Cake %Patch: The number of cells in the radial direction of spherical
   * patches.
   */
  int cake_radial_cells;

  /**
   * Constructs a PatchTransformations object by filling its patch data members
   * using information from `params.ccl`
   */
  PatchTransformations();

  PatchTransformations(const PatchTransformations &) = default;
  PatchTransformations(PatchTransformations &&) = default;
  PatchTransformations &operator=(const PatchTransformations &) = default;
  PatchTransformations &operator=(PatchTransformations &&) = default;

  std_tuple<int, vec<CCTK_REAL, dim> > (*global2local)(
      const PatchTransformations &pt, const vec<CCTK_REAL, dim> &x) = 0;

  vec<CCTK_REAL, dim> (*local2global)(const PatchTransformations &pt, int patch,
                                      const vec<CCTK_REAL, dim> &a) = 0;

  // Calculating global derivatives d/dx from local derivatives d/da
  // requires the Jacobian da/dx, and also its derivative d^2/dx^2 for
  // second derivatives

  // da/dx[i,j] = da[i] / dx[j]
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> > (
      *dlocal_dglobal)(const PatchTransformations &pt, int patch,
                       const vec<CCTK_REAL, dim> &a) = 0;

  // d^2a/dx^2[i,j,k] = d^2a[i] / dx^2[j,k]
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
            vec<smat<CCTK_REAL, dim>, dim> > (*d2local_dglobal2)(
      const PatchTransformations &pt, int patch,
      const vec<CCTK_REAL, dim> &a) = 0;

  // Device functions mirroring the function above
  std_tuple<int, vec<CCTK_REAL, dim> > (*global2local_device)(
      const PatchTransformations &pt, const vec<CCTK_REAL, dim> &x) = 0;
  vec<CCTK_REAL, dim> (*local2global_device)(const PatchTransformations &pt,
                                             int patch,
                                             const vec<CCTK_REAL, dim> &a) = 0;
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> > (
      *dlocal_dglobal_device)(const PatchTransformations &pt, int patch,
                              const vec<CCTK_REAL, dim> &a) = 0;
  std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
            vec<smat<CCTK_REAL, dim>, dim> > (*d2local_dglobal2_device)(
      const PatchTransformations &pt, int patch,
      const vec<CCTK_REAL, dim> &a) = 0;
};

struct PatchSystem {
  std::string name;
  std::vector<Patch> patches;
  int num_patches() const { return patches.size(); }

  PatchTransformations transformations;

  PatchSystem() {}
  PatchSystem(std::string name, std::vector<Patch> patches,
              PatchTransformations transformations)
      : name(std::move(name)), patches(std::move(patches)),
        transformations(std::move(transformations)) {}
};

////////////////////////////////////////////////////////////////////////////////

PatchSystem SetupCartesian();
PatchSystem SetupCubedSphere();
PatchSystem SetupSwirl();
PatchSystem SetupCake();

extern std::unique_ptr<PatchSystem> the_patch_system;

} // namespace MultiPatch

#endif // #ifndef MULTIPATCH_HXX
