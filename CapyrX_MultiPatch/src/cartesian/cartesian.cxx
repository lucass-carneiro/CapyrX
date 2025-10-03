#include "cartesian.hxx"
#include "multipatch.hxx"
#include "patch_systems.hxx"
#include "type_aliases.hxx"

#include <random>
#include <cmath>

namespace CapyrX::MultiPatch::Cartesian {

CCTK_HOST CCTK_DEVICE auto global2local(const PatchParams &,
                                        const svec_t &global_coords)
    -> std_tuple<int, svec_t> {
  return std_make_tuple(0, global_coords);
}

CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &p, int patch,
                                        const svec_t &local_coords) -> svec_t {
  const auto x_dx{dlocal_dglobal(p, patch, local_coords)};
  return std::get<0>(x_dx);
}

CCTK_HOST CCTK_DEVICE auto dlocal_dglobal(const PatchParams &p, int patch,
                                          const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t> {
  const auto x_dx_ddx{d2local_dglobal2(p, patch, local_coords)};
  return std_make_tuple(std::get<0>(x_dx_ddx), std::get<1>(x_dx_ddx));
}

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &, int,
                                            const svec_t &local_coords)
    -> std_tuple<svec_t, jac_t, djac_t> {
  using namespace Arith;
  return std_make_tuple(local_coords, zero<vec<vec<CCTK_REAL, dim>, dim> >()(),
                        zero<vec<smat<CCTK_REAL, dim>, dim> >()());
}

auto make_system(const PatchParams &p) -> PatchSystem {
  const PatchFace ob{.is_outer_boundary = true, .other_patch = -1};

  const Patch patch{.name = "Cartesian",
                    .ncells = {p.ncells_i, p.ncells_j, p.ncells_k},
                    .xmin = {p.xmin, p.ymin, p.zmin},
                    .xmax = {p.xmax, p.ymax, p.zmax},
                    .is_cartesian = true,
                    .faces = {{ob, ob, ob}, {ob, ob, ob}}};

  return PatchSystem{.name = "Cartesian",
                     .id_tag = PatchSystems::cartesian,
                     .patches = {patch}};
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
  real_dist y_dist{par.ymin, par.ymax};
  real_dist z_dist{par.zmin, par.zmax};

  bool all_pass{true};

  // local2global(global2local(global)) == global ?
  for (CCTK_INT i = 0; i < repetitions; i++) {
    const svec_t g_i{x_dist(engine), y_dist(engine), z_dist(engine)};
    const auto l{Cartesian::global2local(par, g_i)};
    const auto g_f{
        Cartesian::local2global(par, std::get<0>(l), std::get<1>(l))};

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
    const svec_t l_i{x_dist(engine), y_dist(engine), z_dist(engine)};

    const auto g{Cartesian::local2global(par, 0, l_i)};
    const auto l{Cartesian::global2local(par, g)};

    const auto &p{std::get<0>(l)};
    const auto &l_f{std::get<1>(l)};

    const auto passed{isapprox(l_i(0), l_f(0)) && isapprox(l_i(1), l_f(1)) &&
                      isapprox(l_i(2), l_f(2)) && p == 0};

    if (!passed) {
      CCTK_VINFO("global2local(local2global(local)) == global repetition %i "
                 "\033[1;31mFAILED\033[0m. Expected (%.16f, %.16f, %.16f) "
                 "but got (%.16f, %.16f, %.16f)",
                 i, l_i(0), l_i(1), l_i(2), l_f(0), l_f(1), l_f(2));
      all_pass = false;
    }
  }

  return all_pass;
}

} // namespace CapyrX::MultiPatch::Cartesian
