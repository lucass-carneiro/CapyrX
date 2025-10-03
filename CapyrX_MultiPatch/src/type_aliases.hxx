#ifndef CAPYRX_MULTIPATCH_TYPE_ALIASES_HXX
#define CAPYRX_MULTIPATCH_TYPE_ALIASES_HXX

#include <cctk.h>

#include <vec.hxx>
#include <vect.hxx>
#include <mat.hxx>

namespace CapyrX::MultiPatch {

constexpr std::size_t dim{3};

struct PatchFace;

/**
 * Type alias for representing spatial vectors. Used for representing
 * coordinate tuples.
 */
using svec_t = Arith::vec<CCTK_REAL, dim>;

/**
 * Type alias for representing Jacobian components
 */
using jac_t = Arith::vec<Arith::vec<CCTK_REAL, dim>, dim>;

/**
 * Type alias for representing Jacobian derivative components
 */
using djac_t = Arith::vec<Arith::smat<CCTK_REAL, dim>, dim>;

/**
 * Type alias for representing an array of integers
 */
using ivec_t = Arith::vect<int, dim>;

/**
 * Type alias for representing an array of reals
 */
using rvec_t = Arith::vect<CCTK_REAL, dim>;

/**
 * Type alias for representing an array of patch faces
 */
using faces_t = Arith::vect<Arith::vect<PatchFace, dim>, 2>;

} // namespace CapyrX::MultiPatch

#endif // CAPYRX_MULTIPATCH_TYPE_ALIASES_HXX