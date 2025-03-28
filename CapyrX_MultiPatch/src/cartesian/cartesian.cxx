#include "cartesian.hxx"
#include "multipatch.hxx"
#include "patch_systems.hxx"
#include "type_aliases.hxx"

namespace CapyrX::MultiPatch::Cartesian {

CCTK_HOST CCTK_DEVICE auto
global2local(const PatchParams &,
             const svec_t &global_coords) -> std_tuple<int, svec_t> {
  return std_make_tuple(0, global_coords);
}

CCTK_HOST CCTK_DEVICE auto local2global(const PatchParams &p, int patch,
                                        const svec_t &local_coords) -> svec_t {
  const auto x_dx{dlocal_dglobal(p, patch, local_coords)};
  return std::get<0>(x_dx);
}

CCTK_HOST CCTK_DEVICE auto
dlocal_dglobal(const PatchParams &p, int patch,
               const svec_t &local_coords) -> std_tuple<svec_t, jac_t> {
  const auto x_dx_ddx{d2local_dglobal2(p, patch, local_coords)};
  return std_make_tuple(std::get<0>(x_dx_ddx), std::get<1>(x_dx_ddx));
}

CCTK_HOST CCTK_DEVICE auto d2local_dglobal2(const PatchParams &, int, const svec_t &local_coords) -> std_tuple<svec_t, jac_t, djac_t> {
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

} // namespace CapyrX::MultiPatch::Cartesian
