#include "multipatch.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "cartesian.hxx"

namespace CapyrX::MultiPatch {

extern "C" CCTK_INT
MultiPatch1_GetSystemSpecification(CCTK_INT *restrict const npatches) {
  // TODO
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetPatchSpecification(
    const CCTK_INT ipatch, CCTK_INT *restrict const is_cartesian,
    const CCTK_INT size, CCTK_INT *restrict const ncells,
    CCTK_REAL *restrict const xmin, CCTK_REAL *restrict const xmax) {
  // TODO
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetBoundarySpecification2(
    const CCTK_INT ipatch, const CCTK_INT size,
    CCTK_INT *restrict const is_interpatch_boundary) {
  // TODO
  return 0;
}

extern "C" void MultiPatch1_GlobalToLocal2(
    const CCTK_INT npoints, const CCTK_REAL *restrict const globalsx,
    const CCTK_REAL *restrict const globalsy,
    const CCTK_REAL *restrict const globalsz, CCTK_INT *restrict const patches,
    CCTK_REAL *restrict const localsx, CCTK_REAL *restrict const localsy,
    CCTK_REAL *restrict const localsz) {
  // TODO
}

extern "C" void MultiPatch1_LocalToGlobal2(
    const CCTK_INT npoints, const CCTK_INT *restrict const patches,
    const CCTK_REAL *restrict const localsx,
    const CCTK_REAL *restrict const localsy,
    const CCTK_REAL *restrict const localsz, CCTK_REAL *restrict const globalsx,
    CCTK_REAL *restrict const globalsy, CCTK_REAL *restrict const globalsz) {
  // TODO
}

extern "C" int CapyrX_MultiPatch_Setup() {
  DECLARE_CCTK_PARAMETERS;
  return 0;
}

extern "C" void CapyrX_MultiPatch_Coordinates_Setup(CCTK_ARGUMENTS) {
  // No idea why the following lines are broken, but they are
  // DECLARE_CCTK_ARGUMENTSX_CapyrX_MultiPatch_Coordinates_Setup;
  // DECLARE_CCTK_PARAMETERS;
}

extern "C" void CapyrX_MultiPatch_Check_Parameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CapyrX_MultiPatch_Check_Parameters;
  DECLARE_CCTK_PARAMETERS;
}

} // namespace CapyrX::MultiPatch
