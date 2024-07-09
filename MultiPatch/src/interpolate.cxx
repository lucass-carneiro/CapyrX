#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/schedule.hxx"

#include <cstdio>
#include <string>

struct coord_data {
  std::vector<CCTK_REAL> coords{};
  std::vector<CCTK_INT> dims{};

  std::vector<CCTK_INT> i{};
  std::vector<CCTK_INT> j{};
  std::vector<CCTK_INT> k{};

  std::vector<CCTK_REAL> x{};
  std::vector<CCTK_REAL> y{};
  std::vector<CCTK_REAL> z{};

  std::vector<CCTK_INT> patches{};
  std::vector<CCTK_INT> levels{};
  std::vector<CCTK_INT> idxs{};
  std::vector<CCTK_INT> components{};
};

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

  // Step 1: Gather target coodinates
  coord_data src_points{};

  CarpetX::active_levels_t().loop_parallel([&](int patch, int level, int index,
                                               int component,
                                               const cGH *cctkGH) {
    const Loop::GridDescBase grid(cctkGH);
    const std::array<int, dim> centering{0, 0, 0};
    const Loop::GF3D2layout layout(cctkGH, centering);

    const auto &current_patch{the_patch_system->patches.at(patch)};
    const auto &patch_faces{current_patch.faces};

    const std::array<const char *, dim> coord_var_names{
        "CoordinatesX::vcoordx", "CoordinatesX::vcoordy",
        "CoordinatesX::vcoordz"};

    for (auto d = 0; d < dim; d++) {
      auto var_data_ptr{static_cast<const CCTK_REAL *>(
          CCTK_VarDataPtr(cctkGH, 0, coord_var_names[d]))};
      Loop::GF3D2<const CCTK_REAL> var{layout, var_data_ptr};

      // Note: This includes symmetry points
      grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
        const auto is_outer_boundary{
            (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) ||
            (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary)};

        // Skip outer boundaries
        if (is_outer_boundary) {
          return;
        } else {
          src_points.coords.push_back(var(p.I));
          src_points.dims.push_back(d);

          src_points.i.push_back(p.i);
          src_points.j.push_back(p.j);
          src_points.k.push_back(p.k);

          src_points.x.push_back(p.x);
          src_points.y.push_back(p.y);
          src_points.z.push_back(p.z);

          src_points.patches.push_back(patch);
          src_points.levels.push_back(level);
          src_points.idxs.push_back(index);
          src_points.components.push_back(component);
        }
      });
    }
  });

#pragma omp critical
  {
    using std::fopen, std::fclose, std::fprintf;

    CCTK_VINFO("Dumping source points data");

    std::string out_file_name{"src_pts_iter_" +
                              std::to_string(cctkGH->cctk_iteration) +
                              std::string(".txt")};

    auto src_pts_file{fopen(out_file_name.c_str(), "w")};

    fprintf(src_pts_file, "# 1:direction 2:patch 3:level 4:index 5:component "
                          "6:i 7:j 8:k 9:x 10:y 11:z 12:coord\n");

    for (std::size_t i = 0; i < src_points.coords.size(); i++) {
      fprintf(src_pts_file, "%d %d %d %d %d %d %d %d %f %f %f %f\n",
              src_points.dims[i], src_points.patches[i], src_points.levels[i],
              src_points.idxs[i], src_points.components[i], src_points.i[i],
              src_points.j[i], src_points.k[i], src_points.x[i],
              src_points.y[i], src_points.z[i], src_points.coords[i]);
    }

    fclose(src_pts_file);
  }
}