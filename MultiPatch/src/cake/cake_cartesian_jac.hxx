
#ifndef MULTIPATCH_CAKE_CARTESIAN_JAC_HXX
#define MULTIPATCH_CAKE_CARTESIAN_JAC_HXX

#include "cake.hxx"

namespace MultiPatch {
namespace Cake {

CCTK_DEVICE CCTK_HOST inline std_tuple<jac_t, djac_t>
cake_cartesian_jac(const PatchTransformations &, const svec_u &) {
  jac_t J = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  djac_t dJ{};

  dJ(0)(0, 0) = 0;
  dJ(0)(0, 1) = 0;
  dJ(0)(0, 2) = 0;
  dJ(0)(1, 0) = 0;
  dJ(0)(1, 1) = 0;
  dJ(0)(1, 2) = 0;
  dJ(0)(2, 0) = 0;
  dJ(0)(2, 1) = 0;
  dJ(0)(2, 2) = 0;

  dJ(1)(0, 0) = 0;
  dJ(1)(0, 1) = 0;
  dJ(1)(0, 2) = 0;
  dJ(1)(1, 0) = 0;
  dJ(1)(1, 1) = 0;
  dJ(1)(1, 2) = 0;
  dJ(1)(2, 0) = 0;
  dJ(1)(2, 1) = 0;
  dJ(1)(2, 2) = 0;

  dJ(2)(0, 0) = 0;
  dJ(2)(0, 1) = 0;
  dJ(2)(0, 2) = 0;
  dJ(2)(1, 0) = 0;
  dJ(2)(1, 1) = 0;
  dJ(2)(1, 2) = 0;
  dJ(2)(2, 0) = 0;
  dJ(2)(2, 1) = 0;
  dJ(2)(2, 2) = 0;

  return std_make_tuple(J, dJ);
}

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_CARTESIAN_JAC_HXX
