#ifndef CAPYRX_MULTIPATCH_MULTIPATCH_HXX
#define CAPYRX_MULTIPATCH_MULTIPATCH_HXX

#include "type_aliases.hxx"

#include <vector>

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
class PatchSystem {
private:
  /**
   * String repreentation of the name of the patch system
   */
  const char *name{"Patch system"};

  /**
   * A vector containing all patches in the system.
   */
  std::vector<Patch> patches{};

public:
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
  PatchSystem(const char *nm, std::vector<Patch> ptc)
      : name{nm}, patches{ptc} {}

  /**
   * Returns the current number of patches in the system.
   *
   * @par Returns:
   * An `int` containing the number of patches in the system.
   */
  [[nodiscard]] auto get_num_patches() const -> int { return patches.size(); }

  /**
   * Returns the name of the system.
   *
   * @par Returns:
   * A C-style string containing the number of patches in the system.
   */
  [[nodiscard]] auto get_name() const -> const char * { return name; }
};

} // namespace CapyrX::MultiPatch

#endif // #ifndef CAPYRX_MULTIPATCH_MULTIPATCH_HXX
