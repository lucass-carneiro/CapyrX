#ifndef CAPYRX_MULTIPATCH_MULTIPATCH_HXX
#define CAPYRX_MULTIPATCH_MULTIPATCH_HXX

#include "type_aliases.hxx"
#include "patch_systems.hxx"

#include <vector>
#include <memory>

namespace CapyrX::MultiPatch {

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
  bool is_outer_boundary{false};

  /**
   * The index of the neighbouring patch, or -1 if the face is the outer
   *boundary.
   **/
  int other_patch{0};
};

/**
 * Describes a patch in the patch system
 */
struct Patch {
  /**
   * String representation of the name of the patch.
   */
  const char *name{"Patch"};

  /**
   * The number of cells in the logical x, y and z dimentions that the patch wil
   * contain.
   */
  ivec_t ncells{0};

  /**
   * Lower coordinate boundary for each logical x, y, z dimention in the patch
   */
  rvec_t xmin{0, 0, 0};

  /**
   * Upper coordinate boundary for each logical x, y, z dimention in the patch
   */
  rvec_t xmax{0, 0, 0};

  /**
   * Wether or not a patch is cartesian. If true, this flag is used to save time
   * during Jacobian calculations.
   */
  bool is_cartesian{false};

  /**
   * Store the 6 faces of a patch, repreenting its connection with other patches
   */
  faces_t faces{};
};

/**
 * Describes a patch system, that is, a collection of patches and patch
 * transformations
 */
struct PatchSystem {
  /**
   * String repreentation of the name of the patch system
   */
  const char *name{"No System"};

  /**
   * ID tag of the patch system
   */
  PatchSystems id_tag{PatchSystems::none};

  /**
   * A vector containing all patches in the system.
   */
  std::vector<Patch> patches{};
};

extern std::unique_ptr<PatchSystem> g_patch_system;

} // namespace CapyrX::MultiPatch

#endif // #ifndef CAPYRX_MULTIPATCH_MULTIPATCH_HXX
