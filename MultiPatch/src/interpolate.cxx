#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/schedule.hxx"

#include <cstdio>
#include <string>

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

struct PointData {
  CCTK_REAL coord{};
  Loop::PointDesc p;
};

using PointList = std::array<std::vector<PointData>, dim>;

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

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
  using CarpetX::dim;
  using MultiPatch::the_patch_system;

  CCTK_VINFO("MultiPatch1_Interpolate called");

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

  // Step 1: gather source points
  std::map<Location, PointList> source_mapping{};

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

        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          for (int d = 0; d < dim; ++d) {
            // Skip outer boundaries
            const auto is_outer_boundary{
                (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) ||
                (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary)};

            if (is_outer_boundary) {
              return;
            } else {
              PointData pd{vcoords[d](p.I), p};
              source_points[d].push_back(pd);
            }
          }
        });
#pragma omp critical
        { source_mapping[location] = std::move(source_points); }
      });

// Debug step: Dump source points
#pragma omp critical
  {
    using std::fopen, std::fclose, std::fprintf;

    CCTK_VINFO("Dumping source points data");

    for (const auto &varind : varinds) {
      const auto var_name{std::string(CCTK_FullVarName(varind))};
      const auto iter_idx{std::to_string(cctkGH->cctk_iteration)};
      const auto out_file_name{var_name + "_src_pts_iter_" + iter_idx +
                               std::string(".txt")};

      auto src_pts_file{fopen(out_file_name.c_str(), "w")};

      fprintf(src_pts_file, "#1:patch\t"
                            "2:level\t"
                            "3:index\t"
                            "4:component\t"
                            "5:dim\t"
                            "6:i\t"
                            "7:j\t"
                            "8:k\t"
                            "9:x\t"
                            "10:y\t"
                            "11:z\t"
                            "12:coord\n");

      for (const auto &[location, point_list] : source_mapping) {
        for (int d = 0; d < dim; ++d) {
          for (const auto &poit_data : point_list[d]) {
            fprintf(src_pts_file,
                    "%i\t"
                    "%i\t"
                    "%i\t"
                    "%i\t"
                    "%i\t"
                    "%i\t"
                    "%i\t"
                    "%i\t"
                    "%.16f\t"
                    "%.16f\t"
                    "%.16f\t"
                    "%.16f\n",
                    location.patch, location.level, location.index,
                    location.component, d, poit_data.p.i, poit_data.p.j,
                    poit_data.p.k, poit_data.p.x, poit_data.p.y, poit_data.p.z,
                    poit_data.coord);
          }
        }
      }

      fclose(src_pts_file);
    }
  }

  // Step 2:
}

} // namespace MultiPatch