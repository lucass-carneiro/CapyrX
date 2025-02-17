#include "multipatch.hxx"
#include "tests.hxx"

#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <string>
#include <sstream>
#include <random>

namespace MultiPatch {
namespace TwoCubes {

using svec = vec<CCTK_REAL, dim>;
using smat = smat<CCTK_REAL, dim>;
using jac_t = vec<vec<CCTK_REAL, dim>, dim>;
using djac_t = vec<smat, dim>;

enum class patch_piece : int { left_cube = 0, right_cube = 1, unknown = 2 };

inline const char *piece_name(const patch_piece &p) {
  switch (static_cast<int>(p)) {
  case static_cast<int>(patch_piece::left_cube):
    return "left cube";
  case static_cast<int>(patch_piece::right_cube):
    return "right cube";
  default:
    return "unknown";
  }
}

CCTK_DEVICE CCTK_HOST svec local2global(const PatchTransformations &pt,
                                        int patch, const svec &local_vars) {

  svec global_vars = {0.0, 0.0, 0.0};

  const auto a{local_vars(0)};
  const auto b{local_vars(1)};
  const auto c{local_vars(2)};

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};

  const auto ymin{pt.two_cubes_ymin};
  const auto ymax{pt.two_cubes_ymax};

  const auto zmin{pt.two_cubes_zmin};
  const auto zmax{pt.two_cubes_zmax};

  const auto dy{pt.two_cubes_delta_y};

  switch (patch) {
  case static_cast<int>(patch_piece::left_cube):
    global_vars = {((1 + a) * (xmax - xmin)) / 4. + xmin,
                   ymax - ((-1 + b) * (dy + a * dy - 2 * ymax + 2 * ymin)) / 4.,
                   (zmax + c * zmax + zmin - c * zmin) / 2.};
    break;

  case static_cast<int>(patch_piece::right_cube):
    global_vars = {
        ((3 + a) * xmax + xmin - a * xmin) / 4.,
        ((-1 + a) * (-1 + b) * dy + 2 * (ymax + b * ymax + ymin - b * ymin)) /
            4.,
        (zmax + c * zmax + zmin - c * zmin) / 2.};
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("No local -> global transformations available for patch %s",
                piece_name(static_cast<patch_piece>(patch)));
#else
    assert(0);
#endif
    break;
  }

  return global_vars;
}

static inline patch_piece get_owner_patch(const PatchTransformations &pt,
                                          const svec &global_vars) {
  using MultiPatchTests::isapprox;

  const auto x = global_vars(0);
  const auto y = global_vars(0);
  const auto z = global_vars(0);

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};

  const auto x_bnd{xmin + (xmax - xmin) / 2.0};

  /* Note: Here when x = x_bnd, we return the left cube patch as the owner.
   * This is an arbitrary choice. The x = x_bnd point is ambiguous but it is
   * never used in interpatch interpolation, so we can freely choose
   * what to return. This return is used later when unit tests are being
   * performed.
   */
  const bool x_left{x < x_bnd || isapprox(x, x_bnd)};
  const bool x_right{x_bnd < x};

  if (x_left) {
    return patch_piece::left_cube;
  } else if (x_right) {
    return patch_piece::right_cube;
  } else {
#ifndef __CUDACC__
    CCTK_VERROR("Coordinate triplet (%f, %f, %f) cannot be located within "
                "the simulation domain",
                x, y, z);
#else
    assert(0);
#endif
    return patch_piece::unknown;
  }
}

CCTK_DEVICE CCTK_HOST std_tuple<int, svec>
global2local(const PatchTransformations &pt, const svec &global_vars) {
  const auto x{global_vars(0)};
  const auto y{global_vars(1)};
  const auto z{global_vars(2)};

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};

  const auto ymin{pt.two_cubes_ymin};
  const auto ymax{pt.two_cubes_ymax};

  const auto zmin{pt.two_cubes_zmin};
  const auto zmax{pt.two_cubes_zmax};

  const auto dy{pt.two_cubes_delta_y};

  const auto piece = get_owner_patch(pt, global_vars);

  CCTK_REAL a{0}, b{0}, c{0};

  switch (static_cast<int>(piece)) {
  case static_cast<int>(patch_piece::left_cube):
    a = 3 + (4 * (x - xmax)) / (xmax - xmin);
    b = 1 + (2 * (xmax - xmin) * (y - ymax)) /
                (2 * dy * (-x + xmin) + (xmax - xmin) * (ymax - ymin));
    c = (-2 * z + zmax + zmin) / (-zmax + zmin);
    break;

  case static_cast<int>(patch_piece::right_cube):
    a = 1 + (4 * (x - xmax)) / (xmax - xmin);
    b = 1 + (2 * (xmax - xmin) * (y - ymax)) /
                (2 * dy * (x - xmax) + (xmax - xmin) * (ymax - ymin));
    c = (-2 * z + zmax + zmin) / (-zmax + zmin);
    break;

  default:
#ifndef __CUDACC__
    CCTK_VERROR("At point (%f, %f, %f): No global -> local transformations "
                "available for patch %s.",
                x, y, z, piece_name(piece));
#else
    assert(0);
#endif
    break;
  }

  return std_make_tuple(static_cast<int>(piece), svec{a, b, c});
}

/**
 * @brief Computes the jacobians and Jacobian derivatives of the Two Cubes
 * patch
 *
 * @param pt The patch transformation structures for the patch
 * @param patch The index of the patch to find the Jacobians for
 * @param global_vars The global point to evaluate the jacobian in.
 * @return The jacobian and jacobian derivative data arrays
 */
static inline std_tuple<jac_t, djac_t>
two_cubes_jacs(const PatchTransformations &pt, int patch,
               const svec &global_vars) {

  using std::pow;
  using std::sqrt;

  jac_t J{};
  djac_t dJ{};

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};

  const auto ymin{pt.two_cubes_ymin};
  const auto ymax{pt.two_cubes_ymax};

  const auto zmin{pt.two_cubes_zmin};
  const auto zmax{pt.two_cubes_zmax};

  const auto dy{pt.two_cubes_delta_y};

  const auto x{global_vars(0)};
  const auto y{global_vars(1)};

  if (patch != static_cast<int>(patch_piece::left_cube) &&
      patch != static_cast<int>(patch_piece::right_cube)) {
#ifndef __CUDACC__
    CCTK_VERROR("No jacobians available for patch %s",
                piece_name(static_cast<patch_piece>(patch)));
#else
    assert(0);
#endif
  }

  if (patch == static_cast<int>(patch_piece::left_cube)) {
    J(0)(0) = 4 / (xmax - xmin);
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)
    (0) = (4 * dy * (xmax - xmin) * (y - ymax)) /
          pow(2 * dy * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 2);
    J(1)
    (1) = (2 * (xmax - xmin)) /
          (2 * dy * (-x + xmin) + (xmax - xmin) * (ymax - ymin));
    J(1)(2) = 0;
    J(2)(0) = 0;
    J(2)(1) = 0;
    J(2)(2) = 2 / (zmax - zmin);
    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = (16 * pow(dy, 2) * (xmax - xmin) * (y - ymax)) /
                  pow(2 * dy * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 3);
    dJ(1)(0, 1) = (4 * dy * (xmax - xmin)) /
                  pow(2 * dy * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = (4 * dy * (xmax - xmin)) /
                  pow(2 * dy * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 2);
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

  } else if (patch == static_cast<int>(patch_piece::right_cube)) {
    J(0)(0) = 4 / (xmax - xmin);
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)
    (0) = (-4 * dy * (xmax - xmin) * (y - ymax)) /
          pow(2 * dy * (x - xmax) + (xmax - xmin) * (ymax - ymin), 2);
    J(1)
    (1) = (2 * (xmax - xmin)) /
          (2 * dy * (x - xmax) + (xmax - xmin) * (ymax - ymin));
    J(1)(2) = 0;
    J(2)(0) = 0;
    J(2)(1) = 0;
    J(2)(2) = 2 / (zmax - zmin);
    dJ(0)(0, 0) = 0;
    dJ(0)(0, 1) = 0;
    dJ(0)(0, 2) = 0;
    dJ(0)(1, 0) = 0;
    dJ(0)(1, 1) = 0;
    dJ(0)(1, 2) = 0;
    dJ(0)(2, 0) = 0;
    dJ(0)(2, 1) = 0;
    dJ(0)(2, 2) = 0;
    dJ(1)(0, 0) = (16 * pow(dy, 2) * (xmax - xmin) * (y - ymax)) /
                  pow(2 * dy * (x - xmax) + (xmax - xmin) * (ymax - ymin), 3);
    dJ(1)(0, 1) = (4 * dy * (-xmax + xmin)) /
                  pow(2 * dy * (x - xmax) + (xmax - xmin) * (ymax - ymin), 2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = (4 * dy * (-xmax + xmin)) /
                  pow(2 * dy * (x - xmax) + (xmax - xmin) * (ymax - ymin), 2);
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
  }

  return std_make_tuple(J, dJ);
}

/**
 * @brief This function computes the local to global coordinate
 * transformation, the jacobian and it's derivative. It can be passed
 * directly to CarpetX.
 *
 * @note The Jacobians are defined as
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 * dJ(i)(j,k) = $dJ^{i}_{j k} = \frac{d^2 a^i}{d x^j d x^k}
 * \right)$.
 *
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation,
 * the local to global jacobian matrix and it's derivative.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t, djac_t>
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const svec &local_vars) {

  const auto local_to_global_result{pt.local2global(pt, patch, local_vars)};
  const auto jacobian_results{
      two_cubes_jacs(pt, patch, local_to_global_result)};

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
} // namespace TwoCubes

/**
 * @brief This function computes the local to global coordinate
 * transformation and the jacobian. It can be passed directly to CarpetX.
 *
 * @note NThe Jacobians is defined as:
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation
 * and the local to global jacobian matrix.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t>
dlocal_dglobal(const PatchTransformations &pt, int patch,
               const svec &local_vars) {

  const auto data = pt.d2local_dglobal2(pt, patch, local_vars);
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

/**
 * @brief Creates a Two Cubes patch piece
 *
 * @tparam p The piece of the patch to make.
 * @param pt The patch transformation object with patch data.
 * @return The constructed patch piece.
 */
template <patch_piece p> Patch make_patch(const PatchTransformations &pt) {
  Patch patch;

  patch.name = piece_name(p);
  patch.xmin = {-1.0, -1.0, -1.0};
  patch.xmax = {1.0, 1.0, 1.0};
  patch.is_cartesian = false;

  PatchFace ob{true, -1};
  PatchFace lc{false, static_cast<int>(patch_piece::left_cube)};
  PatchFace rc{false, static_cast<int>(patch_piece::right_cube)};

  if constexpr (p == patch_piece::left_cube) {
    patch.ncells = {pt.two_cubes_ncells_left, pt.two_cubes_ncells_y,
                    pt.two_cubes_ncells_z};

    patch.faces = {{ob, ob, ob}, {rc, ob, ob}};

  } else if constexpr (p == patch_piece::right_cube) {
    patch.ncells = {pt.two_cubes_ncells_right, pt.two_cubes_ncells_y,
                    pt.two_cubes_ncells_z};

    patch.faces = {{lc, ob, ob}, {ob, ob, ob}};
  }

  return patch;
}

} // namespace TwoCubes

/**
 * @brief Creates the TwoCubes patch system
 *
 * @return PatchSystem object with TwoCubes data and functions
 */
PatchSystem SetupTwoCubes() {
  PatchTransformations pt;
  pt.global2local = &TwoCubes::global2local;
  pt.local2global = &TwoCubes::local2global;
  pt.dlocal_dglobal = &TwoCubes::dlocal_dglobal;
  pt.d2local_dglobal2 = &TwoCubes::d2local_dglobal2;
  pt.global2local_device = &TwoCubes::global2local;
  pt.local2global_device = &TwoCubes::local2global;
  pt.dlocal_dglobal_device = &TwoCubes::dlocal_dglobal;
  pt.d2local_dglobal2_device = &TwoCubes::d2local_dglobal2;

  const auto patches = std::vector<Patch>{
      TwoCubes::make_patch<TwoCubes::patch_piece::left_cube>(pt),
      TwoCubes::make_patch<TwoCubes::patch_piece::right_cube>(pt)};

  return PatchSystem("Two Cubes", std::move(patches), std::move(pt));
}

namespace TwoCubesTests {

/**
 * Tests the get_owner_patch function.
 *
 * @param pt The patch transformations structure
 * @param x A global point to test.
 * @param expected The expected result of the get_owner_patch execution.
 * @return A string indicating the test status.
 */
std::string patch_owner_test(const PatchTransformations &pt,
                             const TwoCubes::svec &x,
                             TwoCubes::patch_piece expected) {

  using namespace MultiPatch::TwoCubes;
  using namespace MultiPatchTests;

  std::ostringstream msg;
  msg << "has ";

  const auto owner_patch = get_owner_patch(pt, x);

  if (owner_patch == expected) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: Expected to get patch "
        << piece_name(expected) << " and got " << piece_name(owner_patch);
  }

  msg << ".";
  return msg.str();
}

/**
 * Tests if local2global(global2local(global)) == global
 *
 * @param pt The patch transformations structure
 * @param global_vars A global point to test.
 * @return A string indicating the test status.
 */
std::string
global_identity_test(const PatchTransformations &pt,
                     const MultiPatch::TwoCubes::svec &global_vars) {

  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::ostringstream msg;
  msg << "has ";

  const auto g2l = pt.global2local(pt, global_vars);

  const auto owner_patch_idx = std::get<0>(g2l);
  const auto local_vars = std::get<1>(g2l);

  const auto l2g = pt.local2global(pt, owner_patch_idx, local_vars);

  const auto test1 = isapprox(l2g(0), global_vars(0));
  const auto test2 = isapprox(l2g(1), global_vars(1));
  const auto test3 = isapprox(l2g(2), global_vars(2));

  if (test1 && test2 && test3) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: ";

    if (!test1) {
      msg << l2g(0) << " =/= " << global_vars(0) << ". ";
    }

    if (!test2) {
      msg << l2g(1) << " =/= " << global_vars(1) << ". ";
    }

    if (!test3) {
      msg << l2g(2) << " =/= " << global_vars(2) << ". ";
    }
  }

  return msg.str();
}

/**
 * Tests if global2local(local2global(local, patch)) == (local, patch)
 *
 * @param pt The patch transformations structure
 * @param patch The patch index to test.
 * @param local_point A local point to test.
 * @return A string indicating the test status.
 */
std::string local_identity_test(const PatchTransformations &pt, int patch,
                                const MultiPatch::TwoCubes::svec &local_point) {

  using MultiPatch::TwoCubes::patch_piece;
  using MultiPatch::TwoCubes::piece_name;
  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::ostringstream msg;
  msg << "has ";

  const auto l2g = pt.local2global(pt, patch, local_point);
  const auto g2l = pt.global2local(pt, l2g);

  const auto computed_patch_idx = std::get<0>(g2l);
  const auto computed_local_point = std::get<1>(g2l);

  const bool test1 = computed_patch_idx == patch;
  const bool test2 = isapprox(computed_local_point(0), local_point(0));
  const bool test3 = isapprox(computed_local_point(1), local_point(1));
  const bool test4 = isapprox(computed_local_point(2), local_point(2));

  if (test1 && test2 && test3 && test4) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: ";

    if (!test1) {
      msg << "the computed patch is "
          << piece_name(static_cast<patch_piece>(computed_patch_idx)) << ". ";
    }

    if (!test2) {
      msg << "the computed local coordinate a is " << computed_local_point(0)
          << ". ";
    }

    if (!test3) {
      msg << "the computed local coordinate a is " << computed_local_point(1)
          << ". ";
    }

    if (!test4) {
      msg << "the computed local coordinate a is " << computed_local_point(2)
          << ". ";
    }
  }

  return msg.str();
}

} // namespace TwoCubesTests
} // namespace MultiPatch

extern "C" void MultiPatch_run_two_cubes_tests(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using MultiPatch::SetupTwoCubes;
  using MultiPatch::TwoCubes::patch_piece;
  using MultiPatch::TwoCubes::svec;

  using MultiPatch::TwoCubesTests::global_identity_test;
  using MultiPatch::TwoCubesTests::local_identity_test;
  using MultiPatch::TwoCubesTests::patch_owner_test;

  using MultiPatchTests::random_seed;
  using std::mt19937_64;
  using std::uniform_int_distribution;
  using std::uniform_real_distribution;

  const auto ps{SetupTwoCubes()};
  const auto &pt{ps.transformations};

  const auto xmin{pt.two_cubes_xmin};
  const auto xmax{pt.two_cubes_xmax};
  const auto xmid{xmin + (xmax - xmin) / 2.0};

  const auto ymin{pt.two_cubes_ymin};
  const auto ymax{pt.two_cubes_ymax};
  const auto ymid{ymin + (ymax - ymin) / 2};

  const auto zmin{pt.two_cubes_zmin};
  const auto zmax{pt.two_cubes_zmax};
  const auto zmid{zmin + (zmax - zmin) / 2};

  const auto xcent_l{xmin + (xmax - xmin) / 4.0};
  const auto xcent_r{xmax - (xmax - xmin) / 4.0};

  CCTK_INFO("Running Two Cubes patch tests:");

  /*
   * Fixed point tests
   */
  const std::array<std::pair<svec, patch_piece>, 35> owner_test_data = {
      // X variation
      std::make_pair(svec{xmin, ymid, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymid, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xmid, ymid, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xcent_r, ymid, zmid}, patch_piece::right_cube),
      std::make_pair(svec{xmax, ymid, zmid}, patch_piece::right_cube),

      // Y variation left
      std::make_pair(svec{xcent_l, ymin, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymid, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymax, zmid}, patch_piece::left_cube),

      // Y variation right
      std::make_pair(svec{xcent_r, ymin, zmid}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymid, zmid}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymax, zmid}, patch_piece::right_cube),

      // Z variation left
      std::make_pair(svec{xcent_l, ymid, zmin}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymid, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymid, zmax}, patch_piece::left_cube),

      // Z variation right
      std::make_pair(svec{xcent_r, ymid, zmin}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymid, zmid}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymid, zmax}, patch_piece::right_cube),

      // Corners
      std::make_pair(svec{xmax, ymax, zmax}, patch_piece::right_cube),
      std::make_pair(svec{xmax, ymin, zmax}, patch_piece::right_cube),
      std::make_pair(svec{xmin, ymin, zmax}, patch_piece::left_cube),
      std::make_pair(svec{xmin, ymax, zmax}, patch_piece::left_cube),

      std::make_pair(svec{xmax, ymax, zmin}, patch_piece::right_cube),
      std::make_pair(svec{xmax, ymin, zmin}, patch_piece::right_cube),
      std::make_pair(svec{xmin, ymin, zmin}, patch_piece::left_cube),
      std::make_pair(svec{xmin, ymax, zmin}, patch_piece::left_cube),

      // Exteriors
      std::make_pair(svec{xcent_l, ymid, zmax + 1.0}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymid, zmin - 1.0}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymax + 1, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xcent_l, ymin - 1, zmid}, patch_piece::left_cube),
      std::make_pair(svec{xmin - 1, ymid, zmid}, patch_piece::left_cube),

      std::make_pair(svec{xcent_r, ymid, zmax + 1.0}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymid, zmin - 1.0}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymax + 1, zmid}, patch_piece::right_cube),
      std::make_pair(svec{xcent_r, ymin - 1, zmid}, patch_piece::right_cube),
      std::make_pair(svec{xmax + 1, ymid, zmid}, patch_piece::right_cube)};

  // Tests if the patch owner is correct.
  for (const auto &data : owner_test_data) {
    CCTK_VINFO("  Patch owner test at point (%f, %f, %f) %s", data.first(0),
               data.first(1), data.first(2),
               patch_owner_test(pt, data.first, data.second).c_str());
  }

  /*
   * Random point tests
   */
  mt19937_64 engine(random_seed);

  uniform_real_distribution<CCTK_REAL> x_distrib(xmin, xmax);
  uniform_real_distribution<CCTK_REAL> y_distrib(ymin, ymax);
  uniform_real_distribution<CCTK_REAL> z_distrib(zmin, zmax);
  uniform_real_distribution<CCTK_REAL> local_distrib(-1, 1);
  uniform_int_distribution<int> patch_distrib(0, 1);

  auto global_point = svec{0, 0, 0};
  auto local_point = svec{0, 0, 0};
  int patch = 0;

  // Tests if local2global(global2local(global)) == global
  for (int i = 0; i < test_repetitions; i++) {
    global_point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};

    CCTK_VINFO("  local2global(global2local(global)) transformation test at "
               "point (%f, %f, %f) %s",
               global_point(0), global_point(1), global_point(2),
               global_identity_test(pt, global_point).c_str());
  }

  // Tests if global2local(local2global(local, patch)) == (local, patch)
  for (int i = 0; i < test_repetitions; i++) {
    local_point = {local_distrib(engine), local_distrib(engine),
                   local_distrib(engine)};

    patch = patch_distrib(engine);

    CCTK_VINFO("  global2local(local2global(local, patch)) transformation test "
               "at point (%f, %f, %f) patch %s %s",
               local_point(0), local_point(1), local_point(2),
               piece_name(static_cast<patch_piece>(patch)),
               local_identity_test(pt, patch, local_point).c_str());
  }

  // TODO: Implement jacobian tests
}
