#include "two_cubes.hxx"

#include <cmath>
#include <random>

namespace CapyrX::MultiPatch::TwoCubes {

enum class PatchPiece : int { left = 0, right = 1 };

static inline CCTK_HOST CCTK_DEVICE auto
get_owner_patch(const PatchParams &par, const svec_t &global_coords)
    -> PatchPiece {
  const auto x{global_coords(0)};
  const auto interface{(par.xmax + par.xmin) / 2.0};

  if (x <= interface) {
    return PatchPiece::left;
  } else {
    return PatchPiece::right;
  }
}

CCTK_HOST CCTK_DEVICE auto global2local(const PatchParams &par,
                                        const svec_t &global_coords)
    -> std_tuple<int, svec_t> {
  const auto xmin{par.xmin};
  const auto xmax{par.xmax};

  const auto ymin{par.ymin};
  const auto ymax{par.ymax};

  const auto zmin{par.zmin};
  const auto zmax{par.zmax};

  const auto bend{par.bend};

  const auto x{global_coords(0)};
  const auto y{global_coords(1)};
  const auto z{global_coords(2)};

  const auto patch{get_owner_patch(par, global_coords)};

  svec_t local_coords{0.0, 0.0, 0.0};

  switch (patch) {

  case PatchPiece::left:
    local_coords(0) = 3 + (4 * (x - xmax)) / (xmax - xmin);
    local_coords(1) =
        1 + (2 * (xmax - xmin) * (y - ymax)) /
                (2 * bend * (-x + xmin) + (xmax - xmin) * (ymax - ymin));
    local_coords(2) = (-2 * z + zmax + zmin) / (-zmax + zmin);
    break;

  case PatchPiece::right:
    local_coords(0) = 1 + (4 * (x - xmax)) / (xmax - xmin);
    local_coords(1) =
        1 + (2 * (xmax - xmin) * (y - ymax)) /
                (2 * bend * (x - xmax) + (xmax - xmin) * (ymax - ymin));
    local_coords(2) = (-2 * z + zmax + zmin) / (-zmax + zmin);
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to compute global2local: Unknown patch piece");
#else
    assert(false);
#endif
    break;
  }

  return std_make_tuple(static_cast<int>(patch), local_coords);
}

CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &par, int patch,
                                        const svec_t &local_coords) -> svec_t {
  const auto xmin{par.xmin};
  const auto xmax{par.xmax};

  const auto ymin{par.ymin};
  const auto ymax{par.ymax};

  const auto zmin{par.zmin};
  const auto zmax{par.zmax};

  const auto bend{par.bend};

  const auto a{local_coords(0)};
  const auto b{local_coords(1)};
  const auto c{local_coords(2)};

  svec_t global_coords = {0.0, 0.0, 0.0};

  switch (static_cast<PatchPiece>(patch)) {

  case PatchPiece::left:
    global_coords(0) = ((1 + a) * (xmax - xmin)) / 4. + xmin;
    global_coords(1) =
        ymax - ((-1 + b) * (bend + a * bend - 2 * ymax + 2 * ymin)) / 4.;
    global_coords(2) = (zmax + c * zmax + zmin - c * zmin) / 2.;
    break;

  case PatchPiece::right:
    global_coords(0) = ((3 + a) * xmax + xmin - a * xmin) / 4.;
    global_coords(1) =
        ((-1 + a) * (-1 + b) * bend + 2 * (ymax + b * ymax + ymin - b * ymin)) /
        4.;
    global_coords(2) = (zmax + c * zmax + zmin - c * zmin) / 2.;
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to compute local2global: Unknown patch piece");
#else
    assert(false);
#endif
    break;
  }

  return global_coords;
}

static inline CCTK_HOST CCTK_DEVICE auto
get_jacs(const PatchParams &par, int patch, const svec_t &global_coords)
    -> std_tuple<jac_t, djac_t> {
  using std::pow;

  jac_t J{};
  djac_t dJ{};

  const auto xmin{par.xmin};
  const auto xmax{par.xmax};

  const auto ymin{par.ymin};
  const auto ymax{par.ymax};

  const auto zmin{par.zmin};
  const auto zmax{par.zmax};

  const auto bend{par.bend};

  const auto x{global_coords(0)};
  const auto y{global_coords(1)};

  switch (static_cast<PatchPiece>(patch)) {

  case PatchPiece::left:
    J(0)(0) = 4 / (xmax - xmin);
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)(0) = (4 * bend * (xmax - xmin) * (y - ymax)) /
              pow(2 * bend * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 2);
    J(1)(1) = (2 * (xmax - xmin)) /
              (2 * bend * (-x + xmin) + (xmax - xmin) * (ymax - ymin));
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
    dJ(1)(0, 0) =
        (16 * pow(bend, 2) * (xmax - xmin) * (y - ymax)) /
        pow(2 * bend * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 3);
    dJ(1)(0, 1) =
        (4 * bend * (xmax - xmin)) /
        pow(2 * bend * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) =
        (4 * bend * (xmax - xmin)) /
        pow(2 * bend * (-x + xmin) + (xmax - xmin) * (ymax - ymin), 2);
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
    break;

  case PatchPiece::right:
    J(0)(0) = 4 / (xmax - xmin);
    J(0)(1) = 0;
    J(0)(2) = 0;
    J(1)(0) = (-4 * bend * (xmax - xmin) * (y - ymax)) /
              pow(2 * bend * (x - xmax) + (xmax - xmin) * (ymax - ymin), 2);
    J(1)(1) = (2 * (xmax - xmin)) /
              (2 * bend * (x - xmax) + (xmax - xmin) * (ymax - ymin));
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
    dJ(1)(0, 0) = (16 * pow(bend, 2) * (xmax - xmin) * (y - ymax)) /
                  pow(2 * bend * (x - xmax) + (xmax - xmin) * (ymax - ymin), 3);
    dJ(1)(0, 1) = (4 * bend * (-xmax + xmin)) /
                  pow(2 * bend * (x - xmax) + (xmax - xmin) * (ymax - ymin), 2);
    dJ(1)(0, 2) = 0;
    dJ(1)(1, 0) = (4 * bend * (-xmax + xmin)) /
                  pow(2 * bend * (x - xmax) + (xmax - xmin) * (ymax - ymin), 2);
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
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to compute jacobians: Unknown patch piece");
#else
    assert(false);
#endif
    break;
  }

  return std_make_tuple(J, dJ);
}

CCTK_HOST CCTK_DEVICE auto dlocal_dglobal(const PatchParams &par, int patch,
                                          const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t> {
  const auto data{d2local_dglobal2(par, patch, local_coords)};
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &par, int patch,
                                            const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t> {
  const auto local_to_global_result{local2global(par, patch, local_coords)};
  const auto jacobian_results{get_jacs(par, patch, local_to_global_result)};

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
}

static inline auto make_patch(const PatchPiece &p, const PatchParams &par)
    -> Patch {

  const auto twice_overlap = 2 * par.patch_overlap;
  const CCTK_REAL dx{(par.xmax - par.xmin) / par.ncells_x};
  const CCTK_REAL dy{(par.ymax - par.ymin) / par.ncells_y};
  const CCTK_REAL dz{(par.zmax - par.zmin) / par.ncells_z};

  PatchFace ob{true, -1};
  PatchFace lc{false, static_cast<int>(PatchPiece::left)};
  PatchFace rc{false, static_cast<int>(PatchPiece::right)};

  // Common patch data
  Patch patch{};

  patch.ncells = {par.ncells_x + twice_overlap, par.ncells_y + twice_overlap,
                  par.ncells_z + twice_overlap};

  patch.xmin = {
      CCTK_REAL{-1.0} - par.patch_overlap * dx,
      CCTK_REAL{-1.0} - par.patch_overlap * dy,
      CCTK_REAL{-1.0} - par.patch_overlap * dz,
  };

  patch.xmax = {CCTK_REAL{1.0} + par.patch_overlap * dx,
                CCTK_REAL{1.0} + par.patch_overlap * dy,
                CCTK_REAL{1.0} + par.patch_overlap * dz};

  patch.is_cartesian = false;

  // Specific patch data
  switch (p) {

  case PatchPiece::left:
    patch.name = "Left";
    patch.faces = {{ob, ob, ob}, {rc, ob, ob}};
    break;

  case PatchPiece::right:
    patch.name = "Right";
    patch.faces = {{lc, ob, ob}, {ob, ob, ob}};
    break;

  default:
#if !defined(__CUDACC__) && !defined(__HIP_PLATFORM_AMD__) &&                  \
    !defined(__HIP_PLATFORM_HCC__) && !defined(__INTEL_LLVM_COMPILER)
    CCTK_VERROR("Unable to create patch. Unknown patch piece");
#else
    assert(false);
#endif
    break;
  } // namespace CapyrX::MultiPatch::TwoCubes

  return patch;
}

auto make_system(const PatchParams &par) -> PatchSystem {
  return PatchSystem{.name = "Two Cubes",
                     .id_tag = PatchSystems::two_cubes,
                     .patches = {make_patch(PatchPiece::left, par),
                                 make_patch(PatchPiece::right, par)}};
}

template <typename fp_type>
static inline auto isapprox(fp_type x, fp_type y, fp_type atol = 0.0) -> bool {
  using std::abs;
  using std::max;
  using std::sqrt;

  const fp_type rtol{
      atol > 0.0 ? 0.0 : sqrt(std::numeric_limits<fp_type>::epsilon())};
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

auto unit_test(std::size_t repetitions, std::size_t seed,
               const PatchParams &par) -> bool {
  using real_dist = std::uniform_real_distribution<CCTK_REAL>;
  using int_dist = std::uniform_int_distribution<CCTK_INT>;

  std::mt19937 engine{static_cast<unsigned long>(seed)};

  real_dist x_dist{par.xmin, par.xmax};
  real_dist z_dist{par.ymin, par.ymax};
  real_dist y_dist{par.zmin, par.zmax};
  real_dist local_dist{-1.0, 1.0};
  int_dist patch_dist{0, 1};

  bool all_pass{true};

  // local2global(global2local(global)) == global ?
  for (CCTK_INT i = 0; i < repetitions; i++) {
    const auto x{x_dist(engine)};
    const auto y{y_dist(engine)};
    const auto z{z_dist(engine)};

    const svec_t g_i{x, y, z};

    const auto l{global2local(par, g_i)};
    const auto g_f{local2global(par, std::get<0>(l), std::get<1>(l))};

    const auto passed{isapprox(g_i(0), g_f(0)) && isapprox(g_i(1), g_f(1)) &&
                      isapprox(g_i(2), g_f(2))};

    if (!passed) {
      CCTK_VINFO("local2global(global2local(global)) == global repetition %i "
                 "\033[1;31mFAILED\033[0m. Expected (%.16f, %.16f, %.16f) "
                 "but got (%.16f, %.16f, %.16f)",
                 i, g_i(0), g_i(1), g_i(2), g_f(0), g_f(1), g_f(2));
      all_pass = false;
    }
  }

  // global2local(local2global(local)) == global ?
  for (CCTK_INT i = 0; i < repetitions; i++) {
    const int p_i{patch_dist(engine)};
    const svec_t l_i{local_dist(engine), local_dist(engine),
                     local_dist(engine)};

    const auto g{local2global(par, p_i, l_i)};
    const auto l{global2local(par, g)};

    const auto &p_f{std::get<0>(l)};
    const auto &l_f{std::get<1>(l)};

    const auto passed{isapprox(l_i(0), l_f(0)) && isapprox(l_i(1), l_f(1)) &&
                      isapprox(l_i(2), l_f(2)) && p_i == p_f};

    if (!passed) {
      CCTK_VINFO(
          "global2local(local2global(local)) == global repetition %i "
          "\033[1;31mFAILED\033[0m. Expected (%.16f, %.16f, %.16f) patch %i "
          "but got (%.16f, %.16f, %.16f) patch %i",
          i, l_i(0), l_i(1), l_i(2), p_i, l_f(0), l_f(1), l_f(2), p_f);
      all_pass = false;
    }
  }

  return all_pass;
}

} // namespace CapyrX::MultiPatch::TwoCubes
