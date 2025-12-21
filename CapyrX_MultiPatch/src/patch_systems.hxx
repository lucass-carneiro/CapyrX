#ifndef CAPYRX_MULTIPATCH_PATCH_SYSTEMS_HXX
#define CAPYRX_MULTIPATCH_PATCH_SYSTEMS_HXX

namespace CapyrX::MultiPatch {

/**
 * Tag for identifying patch systems.
 */
enum class PatchSystems : int {
  none = 0,
  cartesian,
  cubed_spehre,
  thornburg06
};

} // namespace CapyrX::MultiPatch

#endif // CAPYRX_MULTIPATCH_PATCH_SYSTEMS_HXX