#include "multipatch.hxx"

#include <cctk.h>

namespace CapyrX::MultiPatch {

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
  // TODO
  return;
}

} // namespace CapyrX::MultiPatch
