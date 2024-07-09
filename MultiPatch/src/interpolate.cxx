#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"

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
using SourcePoints = std::array<std::vector<CCTK_REAL>, dim>;

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
  assert(cctkGH_);
  assert(nvars_ >= 0);
  assert(varinds_);
  const cGH *const cctkGH = static_cast<const cGH *>(cctkGH_);
  const std::vector<int> varinds(varinds_, varinds_ + nvars_);

  // Step 0: Check input

  for (const int varind : varinds) {
    assert(varind >= 0);
    const int gi = CCTK_GroupIndexFromVarI(varind);
    assert(gi >= 0);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = varind - v0;
    assert(vi >= 0);
    cGroup group_data;
    const int ierr = CCTK_GroupData(gi, &group_data);
    assert(!ierr);
    assert(group_data.grouptype == CCTK_GF);
    assert(group_data.vartype == CCTK_VARIABLE_REAL);
    assert(group_data.disttype == CCTK_DISTRIB_DEFAULT);
    assert(group_data.dim == dim);
    // TODO: Check centering
  }

  // Step 1: Find coordinates where we need to interpolate

  std::map<Location, SourcePoints> source_mapping;
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

        SourcePoints source_points;
        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          // Skip outer boundaries
          for (int d = 0; d < dim; ++d) {
            if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary)
              return;
            if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary)
              return;
          }

          for (int d = 0; d < dim; ++d)
            source_points[d].push_back(vcoords[d](p.I));
        });
#pragma omp critical
        source_mapping[location] = std::move(source_points);
      });

  // Step 2: Interpolate to these coordinates

  // Gather all coordinates
  std::array<std::vector<CCTK_REAL>, dim> coords;
  for (const auto &[location, source_points] : source_mapping) {
    for (int d = 0; d < dim; ++d) {
      coords[d].insert(coords[d].end(), source_points[d].begin(),
                       source_points[d].end());
    }
  }

  assert(coords[0].size() == coords[1].size() &&
         coords[0].size() == coords[2].size());

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
  Interpolate(cctkGH, coords[0].size(), coords[0].data(), coords[1].data(),
              coords[2].data(), varinds.size(), varinds.data(),
              operations.data(), allow_boundaries, resultptrs.data());

  for (size_t n = 0; n < nvars; ++n) {
    for (size_t i = 0; i < results.at(n).size(); ++i) {
      using std::isfinite;
      const auto x = results.at(n).at(i);
      if (!isfinite(x)) {
        CCTK_VINFO("var=%zu i=%zu coord=[%g,%g,%g] value=%g", n, i,
                   coords[0][i], coords[1][i], coords[2][i], x);
      }
      assert(isfinite(x));
    }
  }

  // Scatter interpolated values
  std::map<Location, std::vector<std::vector<CCTK_REAL> > > result_mapping;
  std::size_t pos = 0;
  for (const auto &[location, source_points] : source_mapping) {
    const std::size_t length = source_points[0].size();
    std::vector<std::vector<CCTK_REAL> > result_values(nvars);
    for (std::size_t n = 0; n < nvars; ++n)
      result_values.at(n).insert(result_values.at(n).begin(),
                                 results.at(n).data() + pos,
                                 results.at(n).data() + pos + length);
    pos += length;
    result_mapping[location] = std::move(result_values);
  }
  assert(pos == results.at(0).size());

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

} // namespace MultiPatch
