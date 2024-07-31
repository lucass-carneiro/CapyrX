#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"

#include <cstddef>
#include <cstdio>
#include <loop.hxx>

#include <cctk.h>

#include <array>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace {
// <https://stackoverflow.com/questions/2590677/how-do-i-combine-hash-values-in-c0x>
constexpr inline std::size_t hash_combine(std::size_t h1, std::size_t h2) {
  return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}
} // namespace

namespace MultiPatch {

struct Location {
  int patch;
  int level;
  int index;
  int component;
};

#ifdef CCTK_DEBUG

struct PointDebugData {
  CCTK_REAL coord{};
  std::array<PatchFace, 2> patch_faces{};
  Loop::PointDesc p;
};

using PointListDebugData = std::array<std::vector<PointDebugData>, dim>;
using PointListDebugMap = std::map<Location, PointListDebugData>;
using InterpResultDebugMap =
    std::map<Location, std::vector<std::vector<CCTK_REAL> > >;

#endif

struct point_location {
  int patch;
  int level;
  int component;
  int i;
  int j;
  int k;
};

struct source_points {
  std::vector<CCTK_REAL> global_xs{};
  std::vector<CCTK_REAL> global_ys{};
  std::vector<CCTK_REAL> global_zs{};

  std::vector<CCTK_REAL> local_as{};
  std::vector<CCTK_REAL> local_bs{};
  std::vector<CCTK_REAL> local_cs{};

  std::vector<point_location> locations{};
};

using PointList = std::array<std::vector<CCTK_REAL>, dim>;

} // namespace MultiPatch

namespace std {
template <> struct equal_to<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::equal_to<std::array<int, 4> >()(
        std::array<int, 4>{x.patch, x.level, x.index, x.component},
        std::array<int, 4>{y.patch, y.level, y.index, y.component});
  }
};
template <> struct less<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::less<std::array<int, 4> >()(
        std::array<int, 4>{x.patch, x.level, x.index, x.component},
        std::array<int, 4>{y.patch, y.level, y.index, y.component});
  }
};
template <> struct hash<MultiPatch::Location> {
  std::size_t operator()(const MultiPatch::Location &x) const {
    return hash_combine(hash_combine(hash_combine(std::hash<int>()(x.patch),
                                                  std::hash<int>()(x.level)),
                                     std::hash<int>()(x.index)),
                        std::hash<int>()(x.component));
  }
};
} // namespace std

namespace MultiPatch {

#ifdef CCTK_DEBUG

static inline void dump_src_pts(const std::vector<CCTK_INT> &varinds,
                                const cGH *cctkGH,
                                const source_points &src_pts) {
  using std::fopen, std::fclose;

  for (const auto &varind : varinds) {
    const auto var_name{std::string(CCTK_FullVarName(varind))};
    const auto iter_idx{std::to_string(cctkGH->cctk_iteration)};
    const auto out_file_name{var_name + "_src_pts_iter_" + iter_idx +
                             std::string(".txt")};

    auto file{fopen(out_file_name.c_str(), "w")};

    fprintf(file, "#1:patch\t"
                  "2:level\t"
                  "3:component\t"
                  "4:i\t"
                  "5:j\t"
                  "6:k\t"
                  "7:a\t"
                  "8:b\t"
                  "9:c\t"
                  "10:vcoordx\t"
                  "11:vcoordy\t"
                  "12:vcoordz\n");

    for (std::size_t I = 0; I < src_pts.global_xs.size(); I++) {
      const auto patch{src_pts.locations[I].patch};
      const auto level{src_pts.locations[I].level};
      const auto component{src_pts.locations[I].component};

      const auto i{src_pts.locations[I].i};
      const auto j{src_pts.locations[I].j};
      const auto k{src_pts.locations[I].k};

      const auto gx{src_pts.global_xs[I]};
      const auto gy{src_pts.global_ys[I]};
      const auto gz{src_pts.global_zs[I]};

      const auto la{src_pts.local_as[I]};
      const auto lb{src_pts.local_bs[I]};
      const auto lc{src_pts.local_cs[I]};

      fprintf(file,
              "%i\t"     // patch
              "%i\t"     // level
              "%i\t"     // component
              "%i\t"     // i
              "%i\t"     // j
              "%i\t"     // k
              "%.16f\t"  // a
              "%.16f\t"  // b
              "%.16f\t"  // c
              "%.16f\t"  // vcoordx
              "%.16f\t"  // vcoordy
              "%.16f\n", // vcoordz
              patch, level, component, i, j, k, la, lb, lc, gx, gy, gz);
    }

    fclose(file);
  }
}

static inline void
dump_interp_results(int varind, const cGH *cctkGH, const source_points &src_pts,
                    const std::vector<CCTK_REAL> &interp_results) {
  using std::fopen, std::fclose;

  const auto var_name{std::string(CCTK_FullVarName(varind))};
  const auto iter_idx{std::to_string(cctkGH->cctk_iteration)};
  const auto out_file_name{var_name + "_interp_result_iter_" + iter_idx +
                           std::string(".txt")};

  auto file{fopen(out_file_name.c_str(), "w")};

  fprintf(file, "#1:patch\t"
                "2:level\t"
                "3:component\t"
                "4:i\t"
                "5:j\t"
                "6:k\t"
                "7:a\t"
                "8:b\t"
                "9:c\t"
                "10:vcoordx\t"
                "11:vcoordy\t"
                "12:vcoordz\t"
                "13:interp_result\n");

  for (std::size_t I = 0; I < src_pts.global_xs.size(); I++) {
    const auto patch{src_pts.locations[I].patch};
    const auto level{src_pts.locations[I].level};
    const auto component{src_pts.locations[I].component};

    const auto i{src_pts.locations[I].i};
    const auto j{src_pts.locations[I].j};
    const auto k{src_pts.locations[I].k};

    const auto gx{src_pts.global_xs[I]};
    const auto gy{src_pts.global_ys[I]};
    const auto gz{src_pts.global_zs[I]};

    const auto la{src_pts.local_as[I]};
    const auto lb{src_pts.local_bs[I]};
    const auto lc{src_pts.local_cs[I]};

    const auto interp_result{interp_results[I]};

    fprintf(file,
            "%i\t"     // patch
            "%i\t"     // level
            "%i\t"     // component
            "%i\t"     // i
            "%i\t"     // j
            "%i\t"     // k
            "%.16f\t"  // a
            "%.16f\t"  // b
            "%.16f\t"  // c
            "%.16f\t"  // vcoordx
            "%.16f\t"  // vcoordy
            "%.16f\t"  // vcoordz
            "%.16f\n", // interp result
            patch, level, component, i, j, k, la, lb, lc, gx, gy, gz,
            interp_result);
  }

  fclose(file);
}

static inline void
dump_point_map(const std::vector<CCTK_INT> &varinds, const cGH *cctkGH,
               const PointListDebugMap &source_mapping_debug_data) noexcept {
  using std::fopen, std::fclose, std::fprintf;

  CCTK_VINFO("Dumping source points data");

  for (const auto &varind : varinds) {
    const auto var_name{std::string(CCTK_FullVarName(varind))};
    const auto iter_idx{std::to_string(cctkGH->cctk_iteration)};
    const auto out_file_name{var_name + "_src_pts_iter_" + iter_idx +
                             std::string(".txt")};

    auto file{fopen(out_file_name.c_str(), "w")};

    fprintf(file, "#1:patch\t"
                  "2:level\t"
                  "3:index\t"
                  "4:component\t"
                  "5:direction\t"
                  "6:NI[direction]\t"
                  "7:patch_faces[0].is_outer_boundary\t"
                  "8:patch_faces[0].other_patch\t"
                  "9:patch_faces[1].is_outer_boundary\t"
                  "10:patch_faces[1].other_patch\t"
                  "11:i\t"
                  "12:j\t"
                  "13:k\t"
                  "14:x\t"
                  "15:y\t"
                  "16:z\t"
                  "17:coord\n");

    for (const auto &[location, point_list] : source_mapping_debug_data) {
      for (int d = 0; d < dim; ++d) {
        for (const auto &point_data : point_list[d]) {
          fprintf(file,
                  "%i\t"     // patch
                  "%i\t"     // level
                  "%i\t"     // index
                  "%i\t"     // component
                  "%i\t"     // direction
                  "%i\t"     // NI[direction]
                  "%c\t"     // patch_faces[0].is_outer_boundary
                  "%i\t"     // patch_faces[0].other_patch
                  "%c\t"     // patch_faces[1].is_outer_boundary
                  "%i\t"     // patch_faces[1].other_patch
                  "%i\t"     // i
                  "%i\t"     // j
                  "%i\t"     // k
                  "%.16f\t"  // x
                  "%.16f\t"  // y
                  "%.16f\t"  // z
                  "%.16f\n", // coord
                  location.patch, location.level, location.index,
                  location.component, d, point_data.p.NI[d],
                  point_data.patch_faces[0].is_outer_boundary ? 'y' : 'n',
                  point_data.patch_faces[0].other_patch,
                  point_data.patch_faces[1].is_outer_boundary ? 'y' : 'n',
                  point_data.patch_faces[1].other_patch, point_data.p.i,
                  point_data.p.j, point_data.p.k, point_data.p.x,
                  point_data.p.y, point_data.p.z, point_data.coord);
        }
      }
    }
    fclose(file);
  }
}

static inline void
dump_interp_result(const std::vector<CCTK_INT> &varinds, const cGH *cctkGH,
                   const InterpResultDebugMap &interp_result_dbg_map) noexcept {
  using std::fopen, std::fclose, std::fprintf;

  CCTK_VINFO("Dumping interp results");

  for (const auto &[location, results] : interp_result_dbg_map) {
    for (std::size_t i = 0; i < results.size(); i++) {
      const auto varind{varinds[i]};
      const auto var_name{std::string(CCTK_FullVarName(varind))};
      const auto iter_idx{std::to_string(cctkGH->cctk_iteration)};
      const auto out_file_name{var_name + "_interp_result_iter_" + iter_idx +
                               std::string(".txt")};

      auto file{fopen(out_file_name.c_str(), "a")};

      for (const auto &result : results[i]) {
        fprintf(file, "%i\t%i\t%i\t%i\t%.16f\n", location.patch, location.level,
                location.index, location.component, result);
      }

      fclose(file);
    }
  }
}

#endif

static inline void interpolate_batch(const CCTK_POINTER_TO_CONST cctkGH_,
                                     const CCTK_INT nvars_,
                                     const CCTK_INT *restrict const varinds_) {
  // Step 0: Check input

  // Function Input checking
  if (cctkGH_ == nullptr) {
    CCTK_VERROR("The cctkGH_ pointer is null. Unable to continue");
  }

  if (nvars_ < 0) {
    CCTK_VERROR("The nvars_ variable is negatie. Unable to continue");
  }

  if (varinds_ == nullptr) {
    CCTK_VERROR("The varinds_ pointer is null. Unable to continue");
  }

  // Cast GH and wrap varinds
  const auto cctkGH{static_cast<const cGH *>(cctkGH_)};
  const std::vector<CCTK_INT> varinds(varinds_, varinds_ + nvars_);

  // Check input varinds validity
  for (const auto &varind : varinds) {
    if (varind < 0) {
      CCTK_VERROR("The varind %i is negative", varind);
    }

    const auto gi{CCTK_GroupIndexFromVarI(varind)};

    if (gi < 0) {
      CCTK_VERROR("The goup index %i for varind %i is negative", gi, varind);
    }

    const auto v0{CCTK_FirstVarIndexI(gi)};

    if (gi < 0) {
      CCTK_VERROR("The first var index %i for varind %i is negative", v0,
                  varind);
    }

    const auto vi{varind - v0};

    if (vi < 0) {
      CCTK_VERROR("The index %i for varind %i is negative", vi, varind);
    }

    cGroup group_data{};
    const auto result{CCTK_GroupData(gi, &group_data)};

    switch (result) {
    case -1:
      CCTK_VERROR("Error while retrieving group data for group index %i in "
                  "varind %i: The "
                  "group index is invalid",
                  gi, varind);
      break;

    case -2:
      CCTK_VERROR("Error while retrieving group data for group index %i in "
                  "varind %i: The group data buffer is null",
                  gi, varind);
      break;

    default:
      break;
    }

    if (group_data.grouptype != CCTK_GF) {
      CCTK_VERROR("The group data member \"grouptype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_GF\"",
                  gi, varind);
    }

    if (group_data.vartype != CCTK_VARIABLE_REAL) {
      CCTK_VERROR("The group data member \"vartype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_VARIABLE_REAL\"",
                  gi, varind);
    }

    if (group_data.disttype != CCTK_DISTRIB_DEFAULT) {
      CCTK_VERROR("The group data member \"disttype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_DISTRIB_DEFAULT\"",
                  gi, varind);
    }

    if (group_data.dim != dim) {
      CCTK_VERROR("The group data member \"dim\" in group index %i in varinds "
                  "%i is not %i",
                  gi, varind, dim);
    }

    // TODO: Check centerings table
  }

  // Step 1: Find coordinates where we need to interpolate

  std::map<Location, PointList> source_mapping{};

#ifdef CCTK_DEBUG
  PointListDebugMap source_mapping_debug_data{};
#endif

  CarpetX::active_levels_t().loop_parallel(
      [&](int patch, int level, int index, int component, const cGH *cctkGH) {
        const Loop::GridDescBase grid(cctkGH);
        const std::array<int, dim> centering{0, 0, 0};
        const Loop::GF3D2layout layout(cctkGH, centering);

        const auto &current_patch{the_patch_system->patches.at(patch)};
        const auto &patch_faces{current_patch.faces};

        const Location location{patch, level, index, component};

        const std::array<Loop::GF3D2<const CCTK_REAL>, dim> vcoords{
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordx"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordy"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordz")))};

        PointList source_points;

#ifdef CCTK_DEBUG
        PointListDebugData source_points_debug{};
#endif

        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          // Skip outer boundaries
          for (int d = 0; d < dim; ++d) {
            if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {
              return;
            }

            if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {
              return;
            }
          }

          for (int d = 0; d < dim; ++d) {
            source_points[d].push_back(vcoords[d](p.I));
#ifdef CCTK_DEBUG
            PointDebugData pd{
                vcoords[d](p.I), {patch_faces[0][d], patch_faces[1][d]}, p};
            source_points_debug[d].push_back(pd);
#endif
          }
        });
#pragma omp critical
        {
          source_mapping[location] = std::move(source_points);
#ifdef CCTK_DEBUG
          source_mapping_debug_data[location] = std::move(source_points_debug);
#endif
        }
      });

#ifdef CCTK_DEBUG
#pragma omp critical
  { dump_point_map(varinds, cctkGH, source_mapping_debug_data); }
#endif

  // Step 2: Interpolate to these coordinates

  // Gather all coordinates
  PointList coords{};
  for (const auto &[location, point_list] : source_mapping) {
    for (int d = 0; d < dim; ++d) {
      coords[d].insert(coords[d].end(), point_list[d].begin(),
                       point_list[d].end());
    }
  }

  assert(coords[0].size() == coords[1].size() &&
         coords[0].size() == coords[2].size() &&
         coords[1].size() == coords[2].size());

  const std::size_t nvars = varinds.size();
  const std::size_t npoints = coords[0].size();

  const std::vector<CCTK_INT> operations(nvars, 0);

  // Allocate memory for values
  std::vector<std::vector<CCTK_REAL> > results(nvars);
  std::vector<CCTK_REAL *> resultptrs(nvars);

  results.reserve(nvars);
  resultptrs.reserve(nvars);

  for (size_t n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }

  // Interpolate
  const bool allow_boundaries = false;
  Interpolate(cctkGH, npoints, coords[0].data(), coords[1].data(),
              coords[2].data(), varinds.size(), varinds.data(),
              operations.data(), allow_boundaries, resultptrs.data());

#ifdef CCTK_DEBUG
  for (size_t n = 0; n < nvars; ++n) {
    for (size_t i = 0; i < results.at(n).size(); ++i) {
      using std::isfinite;
      const auto x = results.at(n).at(i);
      if (!isfinite(x)) {
        CCTK_VERROR("Non finite value encoutered: var=%zu i=%zu "
                    "coord=[%g,%g,%g] value=%g",
                    n, i, coords[0][i], coords[1][i], coords[2][i], x);
      }
    }
  }
#endif

  // Scatter interpolated values
  std::map<Location, std::vector<std::vector<CCTK_REAL> > > result_mapping;

#ifdef CCTK_DEBUG
  InterpResultDebugMap interp_result_dbg_map{};
#endif

  std::size_t pos = 0;
  for (const auto &[location, source_points] : source_mapping) {
    const std::size_t length = source_points[0].size();
    std::vector<std::vector<CCTK_REAL> > result_values(nvars);
    for (std::size_t n = 0; n < nvars; ++n)
      result_values.at(n).insert(result_values.at(n).begin(),
                                 results.at(n).data() + pos,
                                 results.at(n).data() + pos + length);
    pos += length;

#pragma omp critical
    {
#ifdef CCTK_DEBUG
      result_mapping[location] = result_values;
      interp_result_dbg_map[location] = result_values;
#else
      result_mapping[location] = std::move(result_values);
#endif
    }
  }
  assert(pos == results.at(0).size());

#ifdef CCTK_DEBUG
#pragma omp critical
  { dump_interp_result(varinds, cctkGH, interp_result_dbg_map); }
#endif

  // Step 3: Write back results
  CarpetX::active_levels_t().loop_parallel([&](int patch, int level, int index,
                                               int component,
                                               const cGH *cctkGH) {
    const Loop::GridDescBase grid(cctkGH);
    const std::array<int, dim> centering{0, 0, 0};
    const Loop::GF3D2layout layout(cctkGH, centering);

    const auto &current_patch{the_patch_system->patches.at(patch)};
    const auto &patch_faces{current_patch.faces};

    const Location location{patch, level, index, component};
    const std::vector<std::vector<CCTK_REAL> > &result_values =
        result_mapping.at(location);

    for (std::size_t n = 0; n < nvars; ++n) {
      const std::vector<CCTK_REAL> &result_values_n = result_values.at(n);

      const Loop::GF3D2<CCTK_REAL> var(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, varinds.at(n))));

      std::size_t pos = 0;

      // Note: This includes symmetry points
      grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
        // Skip outer boundaries
        for (int d = 0; d < dim; ++d) {
          if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {
            return;
          }

          if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {
            return;
          }
        }

        var(p.I) = result_values_n[pos];

        pos++;
      });
      assert(pos == result_values.at(0).size());
    }
  });
}

static inline void interpolate_single(const CCTK_POINTER_TO_CONST cctkGH_,
                                      const CCTK_INT nvars_,
                                      const CCTK_INT *restrict const varinds_) {
  // Step 0: Check input

  // Function Input checking
  if (cctkGH_ == nullptr) {
    CCTK_VERROR("The cctkGH_ pointer is null. Unable to continue");
  }

  if (nvars_ < 0) {
    CCTK_VERROR("The nvars_ variable is negatie. Unable to continue");
  }

  if (varinds_ == nullptr) {
    CCTK_VERROR("The varinds_ pointer is null. Unable to continue");
  }

  // Cast GH and wrap varinds
  const auto cctkGH{static_cast<const cGH *>(cctkGH_)};
  const std::vector<CCTK_INT> varinds(varinds_, varinds_ + nvars_);

  // Check input varinds validity
  for (const auto &varind : varinds) {
    if (varind < 0) {
      CCTK_VERROR("The varind %i is negative", varind);
    }

    const auto gi{CCTK_GroupIndexFromVarI(varind)};

    if (gi < 0) {
      CCTK_VERROR("The goup index %i for varind %i is negative", gi, varind);
    }

    const auto v0{CCTK_FirstVarIndexI(gi)};

    if (gi < 0) {
      CCTK_VERROR("The first var index %i for varind %i is negative", v0,
                  varind);
    }

    const auto vi{varind - v0};

    if (vi < 0) {
      CCTK_VERROR("The index %i for varind %i is negative", vi, varind);
    }

    cGroup group_data{};
    const auto result{CCTK_GroupData(gi, &group_data)};

    switch (result) {
    case -1:
      CCTK_VERROR("Error while retrieving group data for group index %i in "
                  "varind %i: The "
                  "group index is invalid",
                  gi, varind);
      break;

    case -2:
      CCTK_VERROR("Error while retrieving group data for group index %i in "
                  "varind %i: The group data buffer is null",
                  gi, varind);
      break;

    default:
      break;
    }

    if (group_data.grouptype != CCTK_GF) {
      CCTK_VERROR("The group data member \"grouptype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_GF\"",
                  gi, varind);
    }

    if (group_data.vartype != CCTK_VARIABLE_REAL) {
      CCTK_VERROR("The group data member \"vartype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_VARIABLE_REAL\"",
                  gi, varind);
    }

    if (group_data.disttype != CCTK_DISTRIB_DEFAULT) {
      CCTK_VERROR("The group data member \"disttype\" in group index %i in "
                  "varinds %i is not of type \"CCTK_DISTRIB_DEFAULT\"",
                  gi, varind);
    }

    if (group_data.dim != dim) {
      CCTK_VERROR("The group data member \"dim\" in group index %i in varinds "
                  "%i is not %i",
                  gi, varind, dim);
    }

    // TODO: Check centerings table
  }

  // 1: Gather source coordinates and their locations
  source_points src_pts{};

  CarpetX::active_levels_t().loop_parallel(
      [&](int patch, int level, int index, int component, const cGH *cctkGH) {
        const Loop::GridDescBase grid(cctkGH);
        const std::array<int, dim> centering{0, 0, 0};
        const Loop::GF3D2layout layout(cctkGH, centering);

        const auto &current_patch{the_patch_system->patches.at(patch)};
        const auto &patch_faces{current_patch.faces};

        const std::array<Loop::GF3D2<const CCTK_REAL>, dim> vcoords{
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordx"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordy"))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<const CCTK_REAL *>(CCTK_VarDataPtr(
                            cctkGH, 0, "CoordinatesX::vcoordz")))};

        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          // Skip outer boundaries
          for (int d = 0; d < dim; ++d) {
            if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {
              return;
            }

            if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {
              return;
            }
          }

          point_location loc{patch, level, component, p.i, p.j, p.k};

#pragma omp single
          {
            src_pts.global_xs.push_back(vcoords[0](p.I));
            src_pts.global_ys.push_back(vcoords[1](p.I));
            src_pts.global_zs.push_back(vcoords[2](p.I));

            src_pts.local_as.push_back(p.x);
            src_pts.local_bs.push_back(p.y);
            src_pts.local_cs.push_back(p.z);

            src_pts.locations.push_back(loc);
          }
        });
      });

  assert(src_pts.global_xs.size() == src_pts.global_ys.size());
  assert(src_pts.global_xs.size() == src_pts.global_zs.size());
  assert(src_pts.global_ys.size() == src_pts.global_zs.size());

  assert(src_pts.global_xs.size() == src_pts.local_as.size());
  assert(src_pts.global_xs.size() == src_pts.local_bs.size());
  assert(src_pts.global_xs.size() == src_pts.local_cs.size());

  assert(src_pts.global_ys.size() == src_pts.local_as.size());
  assert(src_pts.global_ys.size() == src_pts.local_bs.size());
  assert(src_pts.global_ys.size() == src_pts.local_cs.size());

  assert(src_pts.global_zs.size() == src_pts.local_as.size());
  assert(src_pts.global_zs.size() == src_pts.local_bs.size());
  assert(src_pts.global_zs.size() == src_pts.local_cs.size());

  assert(src_pts.local_as.size() == src_pts.local_bs.size());
  assert(src_pts.local_as.size() == src_pts.local_cs.size());
  assert(src_pts.local_bs.size() == src_pts.local_cs.size());

  const auto num_src_pts{src_pts.global_xs.size()};

#ifdef CCTK_DEBUG
#pragma omp critical
  { dump_src_pts(varinds, cctkGH, src_pts); }
#endif

  // 2: Allocate The results of the interpolation
  std::vector<CCTK_REAL> interp_results(num_src_pts);
  std::vector<CCTK_REAL *> interp_results_ptrs(num_src_pts);

  interp_results.reserve(src_pts.local_as.size());
  interp_results_ptrs.reserve(src_pts.local_as.size());

  std::fill(interp_results.begin(), interp_results.end(), 0.0);
  std::fill(interp_results_ptrs.begin(), interp_results_ptrs.end(), nullptr);

  for (std::size_t i = 0; i < num_src_pts; i++) {
    interp_results_ptrs[i] = &(interp_results[i]);
  }

  // 3: For each function in the group, do
  for (const auto &varind : varinds) {
    // 3.1: Create a sing element array of var indices and operations
    const int varind_array[1] = {varind};
    const int operation_array[1] = {0};

    // 3.2: Interpolate
    Interpolate(cctkGH, num_src_pts, src_pts.global_xs.data(),
                src_pts.global_ys.data(), src_pts.global_zs.data(), 1,
                varind_array, operation_array, false,
                interp_results_ptrs.data());

    // 3.4 Loop over all locations
    CarpetX::active_levels_t().loop_parallel([&](int patch, int level,
                                                 int index, int component,
                                                 const cGH *cctkGH) {
      const Loop::GridDescBase grid(cctkGH);
      const std::array<int, dim> centering{0, 0, 0};
      const Loop::GF3D2layout layout(cctkGH, centering);

      const Loop::GF3D2<CCTK_REAL> var(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, 0, varind)));

      grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
        // 3.5: If this location was stored in the source pointss list, write
        // the interpolated result to the variable at this point
        const auto srx_pt_iter{
            std::find_if(src_pts.locations.begin(), src_pts.locations.end(),
                         [&](const point_location &loc) {
                           return loc.patch == patch && loc.level == level &&
                                  loc.component == component && loc.i == p.i &&
                                  loc.j == p.k && loc.k == p.k;
                         })};

        if (srx_pt_iter == src_pts.locations.end()) {
          return;
        }

        const auto src_pt_idx{
            std::distance(src_pts.locations.begin(), srx_pt_iter)};
        var(p.I) = interp_results[src_pt_idx];
      });
    });

#ifdef CCTK_DEBUG
#pragma omp critical
    { dump_interp_results(varind, cctkGH, src_pts, interp_results); }
#endif
  }
}

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
  interpolate_single(cctkGH_, nvars_, varinds_);
  // interpolate_batch(cctkGH_, nvars_, varinds_);
}

} // namespace MultiPatch
