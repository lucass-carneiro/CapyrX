#include "multipatch.hxx"

#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/interp.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/timer.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#ifdef __CUDACC__
#include <nvtx3/nvToolsExt.h>
#endif

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// System-agnostic owner lookup, provided by this same thorn
// (CapyrX_MultiPatch/src/multipatch.cxx). For a batch of global coordinates it
// returns, per point, the index of the single patch that owns that coordinate
// (the same classifier global2local uses to pick a donor). Declared directly
// rather than through the aliased MultiPatch_GetGlobalToLocal2 registration
// because it lives in the same thorn.
extern "C" void
MultiPatch1_GlobalToLocal2(CCTK_INT npoints, const CCTK_REAL *globalsx,
                           const CCTK_REAL *globalsy, const CCTK_REAL *globalsz,
                           CCTK_INT *patches, CCTK_REAL *localsx,
                           CCTK_REAL *localsy, CCTK_REAL *localsz);

namespace CapyrX::MultiPatch {

struct Location {
  int patch{0};
  int level{0};
  int index{0};
  int component{0};
};

using PointList = std::array<std::vector<CCTK_REAL>, dim>;

} // namespace CapyrX::MultiPatch

// See https://stackoverflow.com/a/2595226
static constexpr inline auto hash_combine(std::size_t h1,
                                          std::size_t h2) -> std::size_t {
  return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}

namespace std {
using namespace CapyrX;

template <> struct equal_to<MultiPatch::Location> {
  bool operator()(const MultiPatch::Location &x,
                  const MultiPatch::Location &y) const {
    return std::equal_to<std::array<int, 4> >()(
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

namespace CapyrX::MultiPatch {

struct ComponentSlice {
  std::size_t offset;
  std::size_t length;
};

struct InterpolationCache {
  // Epoch at which this cache was built. -1 = never built.
  CCTK_INT epoch{-1};

  // AMR-B2: the caller's active (level, patch) range this cache was built for,
  // and part of its key. THE EPOCH ALONE IS NOT A KEY once the range is
  // honoured: `regrid_interpatch_repair` and the syncs that follow it inside
  // one regrid share an epoch and ask for DIFFERENT ranges
  // (`CarpetX/src/schedule.cxx:1503` and `:1916` against `:1975`), and reusing
  // one cache across the two is wrong in both directions -- a narrow cache
  // reused by a wide caller has no slice for the levels the write-back then
  // loops over, so `slices.at()` throws, and a wide cache reused by a narrow
  // caller writes levels the caller excluded, which is the defect this step
  // removes. -1 = never built.
  int min_level{-1}, max_level{-1}, min_patch{-1}, max_patch{-1};

  // AMR-B3: the (patch, level) pairs this cache did NOT collect, because the
  // driver does not call their `CoordinatesX::vertex_coords` interior valid.
  // IT IS PART OF THE KEY, and that is a correctness requirement and not
  // bookkeeping: a cache built while a level was skipped carries no slice for
  // that level, so the write-back's `slices.at(location)` would throw the
  // moment a later caller with the SAME (epoch, range) found the level valid.
  // That is not hypothetical -- it is the regrid path, where
  // `CCTK_Traverse("CCTK_BASEGRID")` runs between the repair and the fill that
  // follows it and turns exactly these pairs valid. Empty = nothing skipped.
  std::vector<std::pair<int, int> > skipped_pl;

  // Ghost-zone source points in deterministic insertion order.
  std::vector<std::pair<Location, PointList> > ordered_components;
  // O(1) slot lookup used during the parallel fill phase.
  std::unordered_map<Location, std::size_t> component_index;

  // Slaved interior cells per component (same slot indexing as
  // ordered_components), in the order their coordinates were appended to the
  // source points. Computed once per epoch; reused by the write-back pass so
  // the owner classification is not recomputed every sync. Empty per slot
  // unless slave_overlap is enabled.
  std::vector<std::vector<Arith::vect<int, dim> > > slaved_indices;

  // Flat coordinate arrays fed to InterpolationSetup.
  PointList coords;

  // The CarpetX interpolation setup (particle container + distribution).
  std::optional<CarpetX::InterpolationSetup> setup;

  // Per-patch outer-boundary policy for InterpolateFromSetup.
  std::vector<Arith::vect<Arith::vect<bool, 3>, 2> > policy;

  // Per-Location (offset, length) into the flat coords/results arrays.
  // Built once per epoch alongside ordered_components; avoids scatter copy.
  std::unordered_map<Location, ComponentSlice> slices;

  // Result buffer reused across calls; resized only on epoch change.
  std::vector<std::vector<CCTK_REAL> > results;
};

static InterpolationCache g_interp_cache;

////////////////////////////////////////////////////////////////////////////////
//
// AMR-B2 INSTRUMENT.  WHICH LEVELS DID THE CALLER ASK FOR, AND WHICH ONES DOES
// THIS FILL TOUCH?
//
// WHY THIS EXISTS.  `MultiPatch_Interpolate` has two callers with two different
// meanings (`CarpetX/src/schedule.cxx`):
//
//   * `SyncGroupsByDirI`, once per sync, under whatever level range the
//     traverse it is inside set -- during evolution `(0, num_levels)`,
//     because nothing subcycles (`schedule.cxx:1975`);
//   * `regrid_interpatch_repair`, once per regrid that modified a level, under
//     `(first_modified_level, last_modified_level + 1)` (`schedule.cxx:1503`
//     in `Initialise`, `:1916` in `Evolve`) -- i.e. under the levels the
//     regrid actually changed, which on a fresh level 1 is the single-element
//     range `{1}`.
//
// The second range is the whole subject of AMR-B2, and NOTHING in this file
// read it before this instrument was written.  The three loops that make up
// the repair disagreed about it: `regrid_interpatch_repair`'s validity filter
// (`schedule.cxx:191`) and its corners-only boundary pass (`:232`) both honour
// `active_levels`, while the interpolation between them -- this function --
// looped every level of every patch.  So the instrument prints BOTH numbers on
// one line: what the driver asked for (`range=`) and what this call will write
// (`touched=`).  On a regrid repair, before AMR-B2, they do not match, and
// that mismatch is the step's failing-before test.
//
// WHAT IT PRINTS, one line per call, all fields on one line so a stream can be
// read with `grep MPLEVELS`:
//
//   MPLEVELS call=N epoch=E src=caller|default range=L[a,b)P[c,d)
//            nlevels=X npatches=Y rebuilt=0|1 reason=first|epoch|range|none
//            nvars=V npoints=P nslaved=S origin_pts=O touched=p:l:n,...
//
//   src         `caller` if the driver's `active_levels` was set, `default`
//               if this function fell back to "every level of every patch".
//               Both callers set it, so `default` is unreachable today; it is
//               printed rather than assumed so that a future caller reaching
//               here outside a traverse is visible instead of silent.
//   range       the caller's range.  `L[1,2)P[0,7)` is "level 1 only, all
//               seven patches".
//   rebuilt     whether THIS call rebuilt the cache, and `reason` why.
//               `first` = the cache had never been built; `epoch` = the AMR
//               epoch moved; `range` = the caller's range moved (AMR-B2 only:
//               before AMR-B2 the range is not part of the cache key and this
//               value cannot appear); `none` = the cache was reused.
//   npoints     the number of cells this call's write-back will write, per
//               variable, on this process.  It is read out of the cache's own
//               slices, so it is what the call WILL do, not an estimate.
//   nslaved     how many of those are slaved interior cells rather than
//               interpatch ghosts.
//   origin_pts  how many query coordinates are EXACTLY (0, 0, 0).  This is
//               `[P246]`'s signature: an interpatch query at the origin means
//               the coordinate it was collected from had not been written when
//               the cache was built.  One legitimate hit is conceivable (the
//               cube's centre vertex is a real coordinate, though never an
//               interpatch ghost), so the number to read is the difference
//               between two columns, not the value itself.
//   coords_invalid
//               how many (patch, level) pairs in `range` hold a
//               `CoordinatesX::vertex_coords` interior the driver does NOT
//               call valid -- i.e. how many of them this call must not read a
//               coordinate from.  AMR-B3's predicate, and AMR-D5's subject.
//               BEFORE AMR-B3's FIX THIS IS WHAT THE CALL READ ANYWAY; AFTER
//               IT, IT IS WHAT THE CALL SKIPPED.  The column does not move
//               across the fix -- `touched` does, and that is the whole test.
//   invalid_pl  those pairs as `patch:level`, or `-`.  Every call outside the
//               regrid path reports `-`: `CCTK_BASEGRID` has run on every
//               level a sync covers.
//   touched     per (patch, level), the point count -- the level census this
//               step is about.  A `:0` entry is printed rather than dropped,
//               so that "the level was in range and had nothing" is
//               distinguishable from "the level was not in range".
//
// COST.  One `std::getenv` per process when off, and nothing else: no line, no
// map, no sweep.  When ON it is one stderr line per call plus one O(npoints)
// pass over the cached coordinates for `origin_pts` and one O(ncomponents)
// pass for `touched` -- 22 lines on `a2_v1_Pno_Sno_R0.par`'s whole run, against
// the 12472-line floods the `CAPYRX_LOG_DONORS` blocks in this file emit for
// one par.  It gets its OWN environment variable (`[P27]`: `CAPYRX_LOG_DONORS`
// already carries eight blocks and cannot be enabled alone) and it is ALWAYS
// COMPILED rather than `CCTK_DEBUG`-gated (AMR-A4's deviation, and for the same
// reason: the numbers have to be readable in the optimized build, which is the
// only build in which the production geometry runs in minutes).
//
// IT IS ON STDERR, DELIBERATELY.  `CCTK_VINFO` is root-rank-only by default --
// the flesh reopens non-root stdout to the null device
// (`Cactus/src/main/CommandLine.c:782-785`) -- and every count here is
// RANK-LOCAL, because which rank holds a component is the decomposition's
// business.  On stdout a four-rank run would silently report rank 0 alone
// (`[P223]`).  The line is assembled in an `ostringstream` and emitted with a
// single `<<`, which the flood instruments below deliberately do NOT do
// (`[P26]`, `[P32]`: their streams are already measured and must stay
// comparable); this one has no history to preserve and is emitted from a
// serial point anyway, since `MultiPatch_Interpolate` is called in global mode.
//
////////////////////////////////////////////////////////////////////////////////

// The active (level, patch) range the CALLER set before invoking us.
//
// `CarpetX::active_levels` is the driver's own record of which levels the
// current traverse applies to (`CarpetX/src/schedule.hxx:93`, defined at
// `schedule.cxx:82`).  It is an `optional` and is empty outside a traverse.
// Both callers of `MultiPatch_Interpolate` set it -- `regrid_interpatch_repair`
// asserts it and `SyncGroupsByDirI` dereferences it several times before
// reaching us -- so the fallback here is unreachable today.  It is NOT a silent
// default: the instrument records which of the two was used.
static CarpetX::active_levels_t caller_active_levels() {
  if (CarpetX::active_levels)
    return *CarpetX::active_levels;
  return CarpetX::active_levels_t();
}

////////////////////////////////////////////////////////////////////////////////
//
// AMR-B3 PREDICATE.  WHICH (patch, level) PAIRS IN THE CALLER'S RANGE HOLD
// VERTEX COORDINATES THIS FUNCTION IS ALLOWED TO READ?
//
// WHY THIS EXISTS.  Both collection passes below read
// `CoordinatesX::vcoordx/y/z`: the ghost pass to get the global coordinate of
// every interpatch ghost point (`grid.loop_bnd`), and
// `collect_slaved_interior` to classify the owner of every interior cell.
// Those arrays are written ONCE per level, `(everywhere)`, by
// `CoordinatesX_Setup` at `CCTK_BASEGRID` (`CoordinatesX/schedule.ccl:3-8`).
// On the regrid path `regrid_interpatch_repair` is called BEFORE that traverse
// -- `CarpetX/src/schedule.cxx:1518` against `:1520` in `Initialise`, `:1931`
// against `:1933` in `Evolve` -- and its active range is exactly the levels
// the regrid modified.  On such a level the coordinate array holds whatever
// the allocator left: poison with `CarpetX::poison_undefined_values = yes`,
// residue without.  This function then reads it as a physical coordinate.
//
// THAT IS AMR-D5, AND IT HAS BEEN MEASURED IN THREE SHAPES:
//
//   * `[P212]` (A3): poison x `slave_overlap` -- `MultiPatch1_GlobalToLocal2`
//     refuses the classified NaN with `Unable to compute global2local:
//     Unknown patch piece`, `cubed_sphere.cxx:177`, exit 1, 4 of 4, ~19 s.
//   * `[P246]` (A5): residue of zeros -- 9126 interpatch queries carry the
//     coordinate (0,0,0), are answered from level 1 at the cube's centre, and
//     the leg exits 0 with no message at all.
//   * `[P257]` (A6): arbitrary residue -- the donor anchor lands outside the
//     box and `CarpetX/src/interpolate.cxx:293`'s
//     `assert(all(i >= 0 && i + order < grid.lsh))` aborts with exit 134,
//     non-deterministically, 13 of 24 attempts on the production geometry, in
//     the OPTIMIZED build with poisoning OFF.
//
// So this is not a poisoning artefact and not a debug-build nuisance: the
// residue decides which of the three you get.
//
// THE PREDICATE IS VALIDITY, AND IT IS D3's PREDICATE VERBATIM.  CarpetX's own
// regrid repair already asks the driver's flags the same question about the
// groups it is about to FILL (`CarpetX/src/schedule.cxx:189-205`, with the
// argument at `:167-175`: "it held poison" and "its interior is not valid" are
// the same statement, because `poison_invalid_gf` poisons exactly the regions
// the flags call invalid).  Before AMR-B3 this file disagreed with its own
// caller: the repair excluded `CoordinatesX` from the variable list BECAUSE
// its interior was not valid, and then this function read that same array for
// every coordinate it interpolated at.  Asking the flags asks the driver's own
// record rather than compiling another thorn's schedule into this one, and it
// fails safe if a fourth caller ever reaches here.
//
// WHAT IT RETURNS.  The (patch, level) pairs in `levels` that EXIST and whose
// `CoordinatesX::vertex_coords` interior the driver does not call valid, in
// ascending (patch, level) order so that two calls' lists can be compared with
// `==` -- which is what makes the list usable as part of the interpolation
// cache's key.  A (patch, level) that does not exist for that patch is NOT in
// the list and is not "skipped": `active_levels_t`'s loops already leave it
// out (`CarpetX/src/schedule.cxx:951-963`), so listing it would make the cache
// key depend on how many levels a patch happens to carry.
//
// COST.  `npatches * nlevels` iterations of a loop over a handful of `bool`s,
// once per `MultiPatch_Interpolate` call: 14 iterations on the seven-patch
// two-level production geometry.  No grid traversal, no allocation beyond the
// (normally empty) result.  It is ALWAYS COMPILED and never env-gated,
// because it is the predicate the fix acts on rather than an instrument.
//
////////////////////////////////////////////////////////////////////////////////
static std::vector<std::pair<int, int> >
invalid_coord_levels(const CarpetX::active_levels_t &levels) {
  std::vector<std::pair<int, int> > invalid;

  // CoordinatesX must be active: both collection passes below fetch
  // `CoordinatesX::vcoordx` BY NAME through `CCTK_VarDataPtr`, so a missing
  // group index is a broken configuration and not an inactive feature.
  // Refusing here rather than returning an empty list follows the same rule as
  // the `MultiPatch_GetBoundarySpecification2` check further down: silently
  // defaulting to "every level's coordinates are valid" would put this defect
  // straight back with nobody able to see it.
  const int coords_gi = CCTK_GroupIndex("CoordinatesX::vertex_coords");
  if (coords_gi < 0)
    CCTK_VERROR("MultiPatch1_Interpolate: CCTK_GroupIndex("
                "\"CoordinatesX::vertex_coords\") returned %d, but this "
                "function reads CoordinatesX::vcoordx/y/z as the coordinates "
                "of every point it interpolates at. Refusing to assume that "
                "every active level's coordinates have been written "
                "(AMR-B3 / AMR-D5)",
                coords_gi);

  for (int patch = levels.min_patch; patch < levels.max_patch; ++patch) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    for (int level = levels.min_level; level < levels.max_level; ++level) {
      if (level >= int(patchdata.leveldata.size()))
        continue; // this patch has no such level, and the loops skip it too
      const auto &leveldata = patchdata.leveldata.at(level);
      const auto &groupdata_ptr = leveldata.groupdata.at(coords_gi);
      if (!groupdata_ptr) {
        // A CCTK_GF group with no level data. Not reachable today; counted as
        // invalid because that is the direction that fails safe.
        invalid.emplace_back(patch, level);
        continue;
      }
      const auto &groupdata = *groupdata_ptr;
      // `sync_tl` exactly as CarpetX's own repair computes it
      // (`schedule.cxx:196-197`): where there is more than one time level the
      // oldest is not synced and is not required to hold a value.
      const int ntls = int(groupdata.mfab.size());
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;
      bool interior_is_valid = true;
      for (int tl = 0; tl < sync_tl; ++tl)
        for (int vi = 0; vi < groupdata.numvars; ++vi)
          if (!groupdata.valid.at(tl).at(vi).get().valid_int)
            interior_is_valid = false;
      if (!interior_is_valid)
        invalid.emplace_back(patch, level);
    }
  }
  return invalid;
}

// AMR-B3: render a (patch, level) list for one log field -- `0:1,0:2`, or `-`
// for the empty list, which is what every call outside the regrid path reports.
static std::string
format_patch_levels(const std::vector<std::pair<int, int> > &pl) {
  if (pl.empty())
    return "-";
  std::ostringstream out;
  bool first = true;
  for (const auto &[patch, level] : pl) {
    if (!first)
      out << ",";
    first = false;
    out << patch << ":" << level;
  }
  return out.str();
}

// AMR-B5: render a (patch, level) -> count census for one log field --
// `0:0:11548,1:0:42602`, or `-` for the empty census.  Same shape as
// `touched=` above so the two fields can be read side by side; `slaved_pl` is
// a SUBSET of `touched`, because `collect_slaved_interior`'s coordinates are
// appended to the same per-component point list the ghost points are in.
static std::string format_census(
    const std::map<std::pair<int, int>, std::size_t> &census) {
  if (census.empty())
    return "-";
  std::ostringstream out;
  bool first = true;
  for (const auto &[patch_level, count] : census) {
    if (!first)
      out << ",";
    first = false;
    out << patch_level.first << ":" << patch_level.second << ":" << count;
  }
  return out.str();
}

static void log_active_levels(
    const CarpetX::active_levels_t &levels, const bool from_caller,
    const CCTK_INT epoch, const bool rebuilt, const char *const reason,
    const std::size_t nvars,
    const std::vector<std::pair<int, int> > &coords_invalid) {
  DECLARE_CCTK_PARAMETERS;

  static const bool log_levels = std::getenv("CAPYRX_LOG_LEVELS") != nullptr;
  if (!log_levels)
    return;

  // BUGFIX_TODO.md B10's rule: an unsynchronised mutable static is written
  // only on the path that also reads it.  This function is called from global
  // mode (`assert(in_global_mode(cctkGH))` in both callers), so the counter is
  // serial; it counts LOGGED calls, which is every call when the instrument is
  // on.
  static long call_counter = 0;
  ++call_counter;

  // AMR-B5: `slaved_pl` is the per-(patch, level) SLAVED census and
  // `slaved_maxlevel` is the deepest level carrying one.  They cost one
  // `std::map` insert per component that has a slaved cell, inside a loop that
  // already walks both vectors, and they exist because the B9 guard below
  // hardcodes `leveldata.at(0)` on the premise that no slaved cell lives above
  // level 0.  `[P217]` measured that premise once, on one geometry, out of an
  // 8.5 GB `CAPYRX_LOG_DONORS` flood; these two fields make it one line of any
  // run.  `slaved_maxlevel = -1` means no slaved cell on this rank at all.
  std::map<std::pair<int, int>, std::size_t> per_patch_level;
  std::map<std::pair<int, int>, std::size_t> slaved_per_patch_level;
  int slaved_maxlevel = -1;
  std::size_t npoints = 0, nslaved = 0;
  for (std::size_t slot = 0;
       slot < g_interp_cache.ordered_components.size(); ++slot) {
    const Location &location = g_interp_cache.ordered_components[slot].first;
    const std::size_t n =
        g_interp_cache.ordered_components[slot].second[0].size();
    per_patch_level[std::make_pair(location.patch, location.level)] += n;
    npoints += n;
    const std::size_t ns = slot < g_interp_cache.slaved_indices.size()
                               ? g_interp_cache.slaved_indices[slot].size()
                               : 0;
    nslaved += ns;
    if (ns > 0) {
      slaved_per_patch_level[std::make_pair(location.patch, location.level)] +=
          ns;
      slaved_maxlevel = std::max(slaved_maxlevel, location.level);
    }
  }

  std::size_t origin_pts = 0;
  for (std::size_t i = 0; i < g_interp_cache.coords[0].size(); ++i)
    if (g_interp_cache.coords[0][i] == 0 && g_interp_cache.coords[1][i] == 0 &&
        g_interp_cache.coords[2][i] == 0)
      ++origin_pts;

  std::ostringstream line;
  line << "MPLEVELS call=" << call_counter << " epoch=" << epoch
       << " src=" << (from_caller ? "caller" : "default") << " range=L["
       << levels.min_level << "," << levels.max_level << ")P["
       << levels.min_patch << "," << levels.max_patch << ")"
       << " nlevels=" << CarpetX::ghext->num_levels()
       << " npatches=" << CarpetX::ghext->num_patches()
       << " rebuilt=" << (rebuilt ? 1 : 0) << " reason=" << reason
       << " nvars=" << nvars << " npoints=" << npoints
       << " nslaved=" << nslaved << " origin_pts=" << origin_pts
       // `[N13]`: an empty census has two causes and they are not the same
       // measurement.  `slave_overlap = no` means the classifier never ran;
       // `slave_overlap = yes` with an empty census means it ran and found
       // nothing.  Say which.
       << " slaving=" << (slave_overlap ? "on" : "off")
       << " slaved_maxlevel=" << slaved_maxlevel
       << " slaved_pl=" << format_census(slaved_per_patch_level)
       << " coords_invalid=" << coords_invalid.size()
       << " invalid_pl=" << format_patch_levels(coords_invalid)
       << " touched=";
  if (per_patch_level.empty()) {
    line << "-";
  } else {
    bool first = true;
    for (const auto &[patch_level, count] : per_patch_level) {
      if (!first)
        line << ",";
      first = false;
      line << patch_level.first << ":" << patch_level.second << ":" << count;
    }
  }
  line << "\n";
  std::cerr << line.str();
}

// One interior overlap-band cell that this patch does not own: its grid index
// (for write-back) and its global coordinate (for the interpolation query).
struct SlavePoint {
  Arith::vect<int, dim> I;
  std::array<CCTK_REAL, dim> x;
};

// Enumerate the interior cells of one component whose true global owner is a
// *different* patch (the "slave" overlap band). These are the cells that today
// are dual-evolved: both this patch and the owner independently RK4-integrate
// them. When slave_overlap is enabled they are treated exactly like ghost
// cells -- overwritten with an interpolated read of the owner's value -- so
// that exactly one numerical solution survives per physical point.
//
// Ownership is resolved through MultiPatch1_GlobalToLocal2, the same
// coordinate-based classifier the donor lookup uses, so every slaved cell
// mirrors straight from its true owner in a single hop (no corner chaining --
// see mp_noise_1.md §4). The returned list has a deterministic order (the
// loop_int traversal order), so the collection and write-back passes, which
// call this independently, agree cell-for-cell.
static std::vector<SlavePoint> collect_slaved_interior(
    const Loop::GridDescBase &grid, const int patch, const int level,
    const int component,
    const std::array<Loop::GF3D2<const CCTK_REAL>, dim> &vcoords) {
  // Gather every interior cell as a candidate, then batch-classify.
  std::vector<Arith::vect<int, dim> > cand_I;
  std::array<std::vector<CCTK_REAL>, dim> cand_x;
  grid.loop_int<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
    cand_I.push_back(p.I);
    for (int d = 0; d < dim; ++d)
      cand_x[d].push_back(vcoords[d](p.I));
  });

  const CCTK_INT ncand = static_cast<CCTK_INT>(cand_I.size());
  std::vector<SlavePoint> slaved;
  if (ncand == 0)
    return slaved;

  std::vector<CCTK_INT> owner(ncand);
  std::array<std::vector<CCTK_REAL>, dim> local; // discarded
  for (int d = 0; d < dim; ++d)
    local[d].resize(ncand);

  MultiPatch1_GlobalToLocal2(ncand, cand_x[0].data(), cand_x[1].data(),
                             cand_x[2].data(), owner.data(), local[0].data(),
                             local[1].data(), local[2].data());

#ifdef CCTK_DEBUG
  // INSTRUMENT (BUGFIX_TODO.md R2 / B10), debug builds only and OFF unless
  // `CAPYRX_LOG_DONORS` is set.  KEPT.
  //
  // mp_slave_3.md §4 instrumentation: log the exact coordinate and owner
  // decision this classifier saw at collection time (which happens once, at
  // the interpolation cache's single serial rebuild pass), for every
  // candidate interior cell. This was compared offline against the same
  // cell's coordinate in the final output TSV (written much later, at I/O
  // time) to test the coordinate-staleness hypothesis: that
  // CoordinatesX::vcoordx/y/z for a wedge's overlap-band cells is not yet in
  // its final settled state when collect_slaved_interior classifies it,
  // so a coordinate that ends up just outside cartesian's cube by output
  // time was seen just inside (or closer to) it here.
  //
  // THAT HYPOTHESIS IS RETRACTED, AND THIS NOTE IS HERE SO NOBODY RE-DERIVES
  // IT.  mp_slave_5.md, confirmed exhaustively over the full 696-cell set in
  // mp_slave_6.md: `CoordinatesX` is never resynced by `MultiPatch_Interpolate`
  // at all and its ghost coordinate is byte-identical from basegrid through
  // cache rebuild to final output, 52504/52504.  The "staleness" that
  // mp_slave_3/4 measured was an unapplied `+nghostzones` index offset in their
  // own analysis scripts.  The instrument is kept because it is cheap when off
  // and because it is the only per-candidate record of the owner decision --
  // not because the result it was built for stands.
  //
  // COST WHEN ON: one stderr line per CANDIDATE INTERIOR CELL of every
  // component, i.e. the whole interior of the grid, and it fires only when
  // `slave_overlap` is on (the collection pass is inside that test).  Measured
  // 8991 lines on `color_slave.par`, 0 on every `slave_overlap = no` leg
  // (evidence/fix/b10/before/h_report.txt).  It shares one environment variable
  // with seven other blocks, so it cannot be enabled alone (`[P27]`).
  {
    static const bool log_donors = std::getenv("CAPYRX_LOG_DONORS") != nullptr;
    if (log_donors) {
      for (CCTK_INT i = 0; i < ncand; ++i) {
        std::cerr << "COLLECT patch=" << patch << " level=" << level
                  << " component=" << component << " I=(" << cand_I[i][0]
                  << "," << cand_I[i][1] << "," << cand_I[i][2] << ")"
                  << " x=" << std::setprecision(17) << cand_x[0][i]
                  << " y=" << cand_x[1][i] << " z=" << cand_x[2][i]
                  << " owner=" << owner[i] << "\n";
      }
    }
  }
#endif

  for (CCTK_INT i = 0; i < ncand; ++i) {
    if (owner[i] != patch) {
      SlavePoint sp;
      sp.I = cand_I[i];
      for (int d = 0; d < dim; ++d)
        sp.x[d] = cand_x[d][i];
      slaved.push_back(sp);
    }
  }
  return slaved;
}

// AMR-B4 (`[P123]`).  DOES THIS FUNCTION WRITE INTERIOR CELLS?
//
// WHY THE DRIVER HAS TO ASK.  With `slave_overlap = yes` the write-back in
// Step 3 does two things, not one: it fills the interpatch ghost points, and it
// OVERWRITES interior cells this patch holds but does not own (the slave band).
// The second write happens after AMReX has already copied those same interior
// cells into the neighbouring boxes' ghost regions, and nothing re-runs that
// copy, so for the rest of the sync every inter-box ghost copy of a slaved cell
// is one slave-write stale -- 400 / 1624 / 5040 cells on `color` /
// `color_ghost` / `color_ghost_overlap`, and exactly zero with slaving off
// (`[P123]`, C1 gate 5b; `[P157]` measured the same signature on the BBH).
// The repair is a `FillBoundary` on the driver's side, because the driver is
// the side that holds the `MultiFab`s.  It is not free, so the driver has to be
// told whether this call writes an interior cell at all.
//
// WHY THIS IS NOT `slave_overlap` READ BY NAME FROM CarpetX.  CarpetX reaches
// this thorn only through the `MultiPatch_*` alias set, and what it needs is a
// statement about the CONTRACT of `MultiPatch_Interpolate` -- "not every cell I
// write is outside your valid region" -- not the name of one implementation's
// boolean.  A different multipatch thorn that writes interior cells for some
// other reason answers this correctly without knowing the word "slave".
//
// IT IS A PARAMETER READ AND NOTHING ELSE, which is load-bearing rather than
// incidental: the driver uses the answer to decide whether to enter a
// COLLECTIVE `FillBoundary`, so the answer must be the same on every rank.
// `slave_overlap` is not steerable, so it is also constant for the run.
extern "C" CCTK_INT MultiPatch1_InterpolateWritesInterior() {
  DECLARE_CCTK_PARAMETERS;
  return static_cast<CCTK_INT>(slave_overlap);
}

extern "C" void
MultiPatch1_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                        const CCTK_INT nvars_,
                        const CCTK_INT *restrict const varinds_) {
#ifdef __CUDACC__
  const nvtxRangeId_t range =
      nvtxRangeStartA("CapyrX::MultiPatch1_Interpolate");
#endif

  static CarpetX::Timer timer("CapyrX::MultiPatch1_Interpolate");
  CarpetX::Interval interval(timer);

  DECLARE_CCTK_PARAMETERS;

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

  // Check input varinds validity.
  //
  // B6: the first non-vertex-centred group found here, if any. The refusal is
  // deferred to the top of Step 2, where `npoints` is known -- see the
  // centering block below, and the comment at the refusal itself.
  int nonvertex_gi{-1};
  CCTK_INT nonvertex_varind{-1};
  std::array<CCTK_INT, dim> nonvertex_centering{0, 0, 0};

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
                  "%i is not %lu",
                  gi, varind, dim);
    }

    // BUGFIX_TODO.md step B6: the check this TODO asked for.
    //
    // Every GF3D2 layout in this function is built with a hard-coded VERTEX
    // centering -- the ghost-point collection pass, the slave-overlap pass and
    // the write-back all say `centering{0, 0, 0}`. CarpetX allocates a group
    // whose indextype[d] is 1 one point SMALLER in direction d
    // (`gash[d] = ash[d] - groupdata.indextype.at(d)` in CarpetX's
    // schedule.cxx; `ash[d] = cctk_ash[d] - indextype[d]` in Loop's loop.hxx),
    // so a vertex
    // layout over such a group has larger `dj`, `dk` and `np` than the
    // allocation: the write-back at the end of this function stores PAST THE
    // END of the array. That is an out-of-bounds write, not a wrong value, and
    // GF3D2's own debug asserts cannot catch it -- they test the vertex
    // `imin`/`imax`, which the loop respects. Measured with valgrind on the
    // `{ccc}` rig named below, before this check existed: 51 x "Invalid write
    // of size 8 ... 8 bytes after a block of size 55,296", all of them at the
    // write-back store, and an uninstrumented optimized run then dies inside
    // CarpetX's own `why_valid_t::set_ghosts`. On a group with several
    // variables most of the overrun is worse than out of bounds and yet
    // invisible: the variables are consecutive inside one fab, so all but the
    // last overrun into the NEXT VARIABLE'S data. The interpolation points are
    // vertex coordinates (`CoordinatesX::vcoord*`) too, so the values would be
    // wrong for a non-vertex group even if the store were in bounds.
    //
    // `SyncGroupsByDirI` flattens every synced group of a sync into ONE call to
    // this function, and the two regrid call sites pass EVERY `CCTK_GF` group
    // on the grid (CarpetX/src/schedule.cxx), so neither a sync nor a regrid
    // can route a non-vertex group around it.
    //
    // The centering is read from the group's centering table, which is the same
    // source CarpetX's own `get_group_indextype` reads
    // (CarpetX/src/driver.cxx), so the two cannot drift apart. No key means
    // vertex-centred, as it does there.
    // Failing-before test: CapyrX_TestMultiPatch/par/centering_ccc.par.
    //
    // The refusal itself is NOT here: it is at the top of Step 2, where
    // `npoints` is known, because it must not fire on a configuration that
    // works today. See the comment at that site.
    std::array<CCTK_INT, dim> centering{0, 0, 0};
    const auto centering_table{CCTK_GroupCenteringTableI(gi)};

    if (centering_table < 0) {
      CCTK_VERROR("Group index %i in varind %i has no centering table (error "
                  "%i); refusing to assume it is vertex-centred",
                  gi, varind, centering_table);
    }

    const auto centering_ret{Util_TableGetIntArray(
        centering_table, static_cast<int>(dim), centering.data(), "centering")};

    if (centering_ret != UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      if (centering_ret != static_cast<int>(dim)) {
        CCTK_VERROR("Could not read the \"centering\" key of group %s (group "
                    "index %i, varind %i): Util_TableGetIntArray returned %i, "
                    "expected %lu",
                    CCTK_FullGroupName(gi), gi, varind, centering_ret, dim);
      }

      for (int d = 0; d < dim; ++d) {
        if (centering[d] != 0) {
          if (nonvertex_gi < 0) {
            nonvertex_gi = gi;
            nonvertex_varind = varind;
            nonvertex_centering = centering;
          }
          break;
        }
      }
    }
  }

  // Step 1: Rebuild the interpolation cache if the AMR epoch or the caller's
  // active level range has changed (AMR-B2).

  const CCTK_INT current_epoch = GetEpoch();

  // AMR-B2: honour the caller's active (level, patch) range, and make it part
  // of the cache key. See `InterpolationCache`'s key comment for why the epoch
  // alone stops being a key the moment the range is honoured.
  const CarpetX::active_levels_t caller_levels = caller_active_levels();
  const bool levels_from_caller = CarpetX::active_levels.has_value();
  const bool cache_never_built = g_interp_cache.epoch < 0;
  const bool epoch_changed = current_epoch != g_interp_cache.epoch;
  const bool range_changed =
      caller_levels.min_level != g_interp_cache.min_level ||
      caller_levels.max_level != g_interp_cache.max_level ||
      caller_levels.min_patch != g_interp_cache.min_patch ||
      caller_levels.max_patch != g_interp_cache.max_patch;

  // AMR-B3: the (patch, level) pairs whose `CoordinatesX::vertex_coords`
  // interior the driver does not call valid. The four loops below skip them --
  // a grid function whose interior is not valid must not be a SOURCE, and this
  // function's source is a coordinate -- and the cache records them, because a
  // cache built while a level was skipped is not the cache a later caller of
  // the same range wants. See `InterpolationCache::skipped_pl`.
  const std::vector<std::pair<int, int> > coords_invalid =
      invalid_coord_levels(caller_levels);
  const auto coords_are_invalid = [&coords_invalid](const int patch,
                                                    const int level) {
    return std::find(coords_invalid.begin(), coords_invalid.end(),
                     std::make_pair(patch, level)) != coords_invalid.end();
  };
  const bool skip_changed = coords_invalid != g_interp_cache.skipped_pl;

  if (epoch_changed || range_changed || skip_changed) {
    CCTK_VINFO("Interpolation cache out of date (cache epoch = %d, active "
               "levels [%d,%d) patches [%d,%d); current epoch = %d, active "
               "levels [%d,%d) patches [%d,%d)). Rebuilding",
               g_interp_cache.epoch, g_interp_cache.min_level,
               g_interp_cache.max_level, g_interp_cache.min_patch,
               g_interp_cache.max_patch, current_epoch,
               caller_levels.min_level, caller_levels.max_level,
               caller_levels.min_patch, caller_levels.max_patch);

    // AMR-B3: say it out loud. A skip that nobody can see is the same defect
    // one indirection further away (`[P135]`, `[P184]`: a zero is only a zero
    // if the instrument spoke). This is bounded by the REBUILD count, not the
    // call count -- one or two lines per regrid.
    if (!coords_invalid.empty())
      CCTK_VINFO("MultiPatch_Interpolate: skipping %zu (patch, level) pair(s) "
                 "%s whose CoordinatesX::vertex_coords interior the driver "
                 "does not call valid. This call interpolates AT vertex "
                 "coordinates, so a level whose coordinates have not been "
                 "written cannot be collected from (AMR-B3 / AMR-D5). Nothing "
                 "on such a level is filled by this call",
                 coords_invalid.size(),
                 format_patch_levels(coords_invalid).c_str());

    // Collect ghost-zone coordinates.
    // Serial pass: assign each component a fixed slot so the parallel fill
    // phase can write without any synchronisation.
    g_interp_cache.ordered_components.clear();
    g_interp_cache.component_index.clear();
    caller_levels.loop_serially(
        [&](int patch, int level, int index, int component, const cGH *) {
          // AMR-B3: no slot, so no later pass visits this (patch, level) --
          // the ghost fill, the slave collection and the write-back all index
          // through `component_index`/`slices`, which are built here.
          if (coords_are_invalid(patch, level))
            return;
          const Location location{patch, level, index, component};
          g_interp_cache.component_index[location] =
              g_interp_cache.ordered_components.size();
          g_interp_cache.ordered_components.emplace_back(location, PointList{});
        });

    // Per-slot slaved-index storage, sized to match ordered_components so the
    // parallel pass can write each slot without synchronisation.
    g_interp_cache.slaved_indices.assign(
        g_interp_cache.ordered_components.size(), {});

    // Parallel pass: each component writes to its pre-assigned slot.
    {
      caller_levels.loop_parallel([&](int patch, int level, int index,
                                     int component, const cGH *cctkGH) {
        if (coords_are_invalid(patch, level))
          return; // AMR-B3: `vcoords` below is exactly the unwritten array
        const Loop::GridDescBase grid(cctkGH);
        const std::array<int, dim> centering{0, 0, 0};
        const Loop::GF3D2layout layout(cctkGH, centering);

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

        const auto &current_patch{g_patch_system->patches.at(patch)};
        const auto &patch_faces{current_patch.faces};

        const Location location{patch, level, index, component};
        const std::size_t slot = g_interp_cache.component_index.at(location);

        PointList source_points;

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
          }

#ifdef CCTK_DEBUG
          // INSTRUMENT (BUGFIX_TODO.md R2 / B10), debug builds only and OFF
          // unless `CAPYRX_LOG_DONORS` is set.  KEPT.
          //
          // mp_slave_3.md §4 instrumentation, part 2: bucket (a)'s victims
          // turned out to be ordinary ghost points (this loop), not members
          // of the `slaved` list (collect_slaved_interior, logged
          // separately as COLLECT) -- so the coordinate-staleness question
          // had to be tested against *this* collection point, the one that
          // actually feeds InterpolationSetup's own donor-routing call for
          // ghost points, not collect_slaved_interior's separate,
          // redundant classification pass.  THAT QUESTION IS CLOSED AND THE
          // ANSWER WAS A RETRACTION: see the note on COLLECT above -- the
          // coordinate is byte-identical from basegrid to output, 52504/52504
          // (mp_slave_5.md, mp_slave_6.md).  Kept for the
          // BASEGRID_COORD x GHOSTCOORD join, not for the retracted result.
          //
          // COST WHEN ON: one line per interpatch ghost point, per component,
          // on a cache rebuild.  Measured 12472 on `color.par`, 31552 on
          // `color_ghost.par` (evidence/fix/b10/before/h_report.txt).
          {
            static const bool log_donors =
                std::getenv("CAPYRX_LOG_DONORS") != nullptr;
            if (log_donors) {
              std::cerr << "GHOSTCOORD patch=" << patch << " level=" << level
                        << " component=" << component << " I=(" << p.I[0]
                        << "," << p.I[1] << "," << p.I[2] << ")"
                        << " x=" << std::setprecision(17) << vcoords[0](p.I)
                        << " y=" << vcoords[1](p.I) << " z=" << vcoords[2](p.I)
                        << "\n";
            }
          }
#endif
        });

        g_interp_cache.ordered_components[slot].second =
            std::move(source_points);
      });
    }

    // Slave overlap band: append the interior cells this patch does not own,
    // in loop_int order, right after the ghost cells collected above. Done in
    // a *serial* pass because collect_slaved_interior calls
    // MultiPatch1_GlobalToLocal2, whose CarpetX timer is only safe to enter
    // single-threaded (its handle is shared across threads). This runs only on
    // a cache rebuild (epoch change), not every sync, so the serial cost is
    // paid rarely. The write-back pass reuses the cached slaved_indices, so the
    // owner classification here is never repeated per-step.
    if (slave_overlap) {
      caller_levels.loop_serially([&](int patch, int level, int index,
                                     int component, const cGH *cctkGH) {
        if (coords_are_invalid(patch, level))
          return; // AMR-B3: this is `[P212]`'s reader -- the classified
                  // coordinate that reads `-nan(0x80000deadbeef)`
        const Loop::GridDescBase grid(cctkGH);
        const std::array<int, dim> centering{0, 0, 0};
        const Loop::GF3D2layout layout(cctkGH, centering);

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

        const Location location{patch, level, index, component};
        const std::size_t slot = g_interp_cache.component_index.at(location);

        const auto slaved =
            collect_slaved_interior(grid, patch, level, component, vcoords);

        PointList &source_points =
            g_interp_cache.ordered_components[slot].second;
        auto &slot_indices = g_interp_cache.slaved_indices[slot];
        slot_indices.reserve(slaved.size());
        for (const auto &sp : slaved) {
          for (int d = 0; d < dim; ++d)
            source_points[d].push_back(sp.x[d]);
          slot_indices.push_back(sp.I);
        }
      });
    }

    // Flatten into coordinate arrays and record per-Location slices.
    g_interp_cache.coords = {};
    g_interp_cache.slices.clear();
    std::size_t flat_offset = 0;
    for (const auto &[location, point_list] :
         g_interp_cache.ordered_components) {
      const std::size_t length = point_list[0].size();
      g_interp_cache.slices[location] = {flat_offset, length};
      flat_offset += length;
      for (int d = 0; d < dim; ++d) {
        g_interp_cache.coords[d].insert(g_interp_cache.coords[d].end(),
                                        point_list[d].begin(),
                                        point_list[d].end());
      }
    }

    assert(g_interp_cache.coords[0].size() == g_interp_cache.coords[1].size() &&
           g_interp_cache.coords[0].size() == g_interp_cache.coords[2].size() &&
           g_interp_cache.coords[1].size() == g_interp_cache.coords[2].size());

    // Build InterpolationSetup (expensive — skipped when epoch is unchanged)
    const std::size_t npoints_cache = g_interp_cache.coords[0].size();
    for (auto &r : g_interp_cache.results)
      r.resize(npoints_cache);
    // C-AMR2 (CarpetX `src/interpolate.cxx`): these query points are the
    // interpatch ghost points of a patch boundary -- and, when
    // `slave_overlap` is on, the non-owned interior overlap cells that are
    // filled the same way. Every one of them must be answered from a level-0
    // box. If a refined level ever covers the region a seam point is drawn
    // from, `Redistribute` will silently start answering it from prolongated
    // fine data instead of evolved coarse data, at a moment set by wherever
    // the refinement boxes have travelled to. CarpetX cannot tell that apart
    // from a legitimate level-1 answer to somebody else's query, so the
    // caller declares it; the flag defaults to false there.
    g_interp_cache.setup.emplace(cctkGH, static_cast<CCTK_INT>(npoints_cache),
                                 g_interp_cache.coords[0].data(),
                                 g_interp_cache.coords[1].data(),
                                 g_interp_cache.coords[2].data(),
                                 /*require_level0_donors=*/true);

    // Build per-patch outer-boundary policy
    const int npatches = cctkGH->cctk_npatches;
    g_interp_cache.policy.resize(npatches);

    static const bool have_boundary_spec =
        CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification2");

    // This routine is itself provided by CapyrX_MultiPatch, which also
    // provides MultiPatch_GetBoundarySpecification2 (interface.ccl). The two
    // are always wired together, so reaching here without the alias means
    // the aliasing is broken (e.g. a malformed interface.ccl), not that
    // multipatch is inactive. Silently defaulting to policy=true would open
    // every outer face to donor stencils reading outer-BC ghosts into
    // interpatch values without anyone noticing.
    if (!have_boundary_spec) {
      CCTK_VERROR(
          "MultiPatch1_Interpolate: MultiPatch_GetBoundarySpecification2 is "
          "not aliased even though CapyrX_MultiPatch (which provides it) is "
          "active. Refusing to silently default every outer face's donor "
          "policy to \"anchor allowed\" (this would let outer-boundary data "
          "leak into interpatch values on every face). Check interface.ccl "
          "wiring.");
    }

    std::array<CCTK_INT, 2 * dim> spec;
    for (int p = 0; p < npatches; ++p) {
      MultiPatch_GetBoundarySpecification2(p, 2 * dim, spec.data());
      for (int f = 0; f < 2; ++f) {
        for (int d = 0; d < dim; ++d) {
          g_interp_cache.policy[p][f][d] = !spec[2 * d + f];
        }
      }
    }

    g_interp_cache.epoch = current_epoch;
    g_interp_cache.min_level = caller_levels.min_level;
    g_interp_cache.max_level = caller_levels.max_level;
    g_interp_cache.min_patch = caller_levels.min_patch;
    g_interp_cache.max_patch = caller_levels.max_patch;
    g_interp_cache.skipped_pl = coords_invalid;
  }

  // AMR-B3 adds a fourth term. THE VOCABULARY OF THE FIRST THREE IS B2's AND
  // IS UNCHANGED -- `first`, `epoch`, `range`, `epoch+range`, `none` -- so a
  // b2 stream and a b3 stream can be read side by side; `skip` joins them with
  // `+` in the order the tests are written.
  std::string rebuild_reason;
  if (cache_never_built) {
    rebuild_reason = "first";
  } else {
    if (epoch_changed)
      rebuild_reason = "epoch";
    if (range_changed)
      rebuild_reason += rebuild_reason.empty() ? "range" : "+range";
    if (skip_changed)
      rebuild_reason += rebuild_reason.empty() ? "skip" : "+skip";
    if (rebuild_reason.empty())
      rebuild_reason = "none";
  }
  log_active_levels(caller_levels, levels_from_caller, current_epoch,
                    epoch_changed || range_changed || skip_changed,
                    rebuild_reason.c_str(), varinds.size(), coords_invalid);

  // Step 2: Interpolate using the cached setup.

  const std::size_t nvars = varinds.size();
  const std::size_t npoints = g_interp_cache.coords[0].size();

  // BUGFIX_TODO.md step B6: refuse a non-vertex-centred group HERE, and not in
  // the input-validation loop above, because `npoints == 0` is a configuration
  // that works today and an additive commit may not break it.
  //
  // `npoints` is the number of cells this call will write, and it is the same
  // set for every group in the call (the target coordinates come from
  // CoordinatesX, not from the group). It is ZERO whenever no cell needs
  // filling -- in particular for the single-patch `Cartesian` patch system,
  // where all six faces are outer boundaries and the collection loop skips
  // every ghost point. That configuration DOES reach this function with
  // non-vertex groups, and not by anyone's choice: the two regrid call sites
  // build their variable list from every `CCTK_GF` group on the grid, so
  // `CoordinatesX::cell_coords`, `CoordinatesX::cell_volume` and
  // `CarpetXRegrid::regrid_error` -- all `{ccc}` -- arrive here on any run with
  // `max_num_levels > 1`. With `npoints == 0` nothing is stored through the
  // wrong layout and there is nothing to refuse; measured on
  // evidence/fix/a8/pars/a8_cart_l2.par, which an earlier revision of this
  // check turned from exit 0 into exit 1.
  //
  // Multi-rank note: `npoints` is this rank's count, so on a decomposition
  // where one rank owns no interpatch ghost cell the refusal is raised by the
  // other ranks. `CCTK_VERROR` aborts the job, so the run still stops; what is
  // not guaranteed is which rank prints it.
  if (nonvertex_gi >= 0 && npoints > 0) {
    CCTK_VERROR(
        "Group %s (group index %i, varind %i) has centering [%i,%i,%i], and "
        "MultiPatch_Interpolate only supports vertex-centred [0,0,0] grid "
        "functions: its interpolation points are vertex coordinates and every "
        "layout it builds assumes a vertex-centred allocation, so filling this "
        "group's %zu interpatch cells would write past the end of its array. "
        "Do not sync this group on a patch system with interpatch faces "
        "(CarpetX::SyncGroupsByDirI passes every synced group of a sync to "
        "this function, and the regrid path passes every CCTK_GF group, so it "
        "cannot be excluded one group at a time)",
        CCTK_FullGroupName(nonvertex_gi), nonvertex_gi, int(nonvertex_varind),
        int(nonvertex_centering[0]), int(nonvertex_centering[1]),
        int(nonvertex_centering[2]), npoints);
  }

  // BUGFIX_TODO.md step B9: REFUSE `slave_overlap` WHEN A VARIABLE BEING
  // FILLED HAS A PHYSICAL OUTER FACE THAT NOTHING WRITES.
  //
  // THE DEFECT.  A slaved write puts an interpolated value into an INTERIOR
  // cell -- a cell CarpetX marks `valid_int` -- and the interpolation's donor
  // stencil can reach the donor patch's physical outer ghost zone.  If nothing
  // writes that ghost zone, the interior takes whatever was there.  With
  // `CarpetX::poison_undefined_values = yes` that is poison and the arithmetic
  // carries it (`0 x NaN = NaN`, so a near-zero interpolation weight does not
  // save it); with poisoning off it is stale memory and the corruption is
  // quiet at the moment it happens.  Measured at `evidence/it7/i7_e5.par`: 32
  // NaN interior cells on the wedges' outermost interior radial vertex plane,
  // per patch 1/1/5/5/10/10 -- a direct fingerprint of `get_owner_patch`'s
  // `x > y > z` tie-break -- and 0 NaN, no message, on the same leg with
  // poisoning off.
  //
  // THE PREDICATE HAS TWO TERMS AND BOTH ARE NECESSARY.
  //
  //   (1) The face carries neither a symmetry nor a boundary condition:
  //       `all_faces_have_symmetries_or_boundaries()`, which is CarpetX's own
  //       function, CALLED and not re-derived.  `boundary_x` and its eleven
  //       siblings do not answer this question, because
  //       `get_group_boundaries` applies the 24
  //       `{dirichlet,linear_extrapolation,neumann,robin}_{,upper_}{x,y,z}_vars`
  //       per-group overrides on top of them; a parameter scan therefore says
  //       "no boundary condition" for a group that has one.  There is no
  //       `none_*_vars`, so that error only ever runs in the direction of
  //       refusing a configuration that works.
  //
  //   (2) The driver has not certified this variable's outer region:
  //       `valid.at(0).at(vi).get().valid_outer` is false.  TERM (1) ALONE IS
  //       NOT ENOUGH, and that is measured, not reasoned: with term (1) alone
  //       this guard refused `CAPYRX_TESTMULTIPATCH::COLOR` on
  //       `smooth_z_none.par` (`evidence/fix/b9/probe/`), a group whose own
  //       thorn declares `WRITES: color(everywhere)` and fills its outer ghost
  //       zone itself.  Nothing is unwritten there and nothing is poisoned;
  //       refusing it would have broken a rig that works, which is the one
  //       thing an additive commit may not do.  The outer validity bit is
  //       precisely the record of such a write, and it is the bit
  //       `poison_invalid_gf` acts on -- so term (2) is not a proxy for the
  //       defect, it is the same bit that decides whether there is anything
  //       poisonous in the donor zone to import.
  //
  // WHAT TERM (2) DOES NOT PROMISE.  `valid_outer` says that something wrote
  // that region and the driver is willing to certify it.  It does not say the
  // contents are current: a group written everywhere at initial and evolved in
  // the interior afterwards keeps the bit while its outer data goes stale.
  // That is CarpetX's own validity semantics, and this guard defers to it
  // rather than inventing a second, stricter notion of validity downstream of
  // the driver that maintains the first.
  //
  // WHY THE PREDICATE HAS A `slave_overlap` TERM AT ALL, which is the part of
  // this guard most likely to be questioned.  WITHOUT slaving the same import
  // lands only in interpatch GHOST cells, and CarpetX already declines to
  // certify those: `SyncGroupsByDirI`'s postcondition calls `set_outer(true)`
  // only if `all_faces_have_symmetries_or_boundaries()` -- term (1) again --
  // and on a patch system the outer bit is what covers the interpatch ghosts.
  // So the non-slaved case writes wrong values into cells that are MARKED
  // INVALID, which is CarpetX working as designed.  Slaving is what moves
  // those values into cells marked valid.  The control for this sentence is
  // B7's `b7_wit_none.par` with slaving off: it must NOT be refused here, and
  // it is not -- it still aborts exactly where it did before, at the
  // `CCTK_DEBUG` `contains_nan()` sweep in `CarpetX/src/schedule.cxx`.
  //
  // WHY HERE AND NOT AT PARAMCHECK.  Term (1) is a `GroupData` member reading
  // `ghext->patchdata[p].symmetries` and `groupdata.boundaries`; term (2) is a
  // runtime bit that a thorn's own write can set.  Neither exists at
  // PARAMCHECK -- there is no grid yet.  This is the first point at which both
  // are known AND slaved cells are about to be written.
  //
  // WHEN IT FIRES, MEASURED RATHER THAN ASSUMED, AND IT IS NOT ALWAYS THE
  // FIRST FILL OF THE RUN.  Term (2) is a runtime bit, so a group whose
  // initial-data routine declares `WRITES: state(everywhere)` -- CapyrX_WaveToy
  // does, `schedule.ccl` -- has a CERTIFIED outer zone at the first fill.  That
  // fill imports real initial data, there is nothing wrong with it, and it is
  // allowed.  The evolution then writes `rhs(interior)` only, the bit goes
  // false, and the NEXT fill is the first one that would import an uncertified
  // zone.  That is the one refused.  Measured on A9's pair (evidence/fix/b9):
  // `ov1_none` (`cctk_itlast = 0`, one fill) runs to completion with its donor
  // census unmoved in every field; `ov1_none_it1` is the same rig one iteration
  // longer and is refused at its SECOND call.
  //
  // The alignment with the defect is exact, and that is the argument for
  // firing here rather than earlier: the 32 NaNs of `evidence/it7/i7_e5.par`
  // appear at `SyncGroupsByDirI call #2`, iteration 1's `ODESolvers_PostStep`
  // sync -- precisely the fill this guard stops.  So do not describe this as a
  // configuration check.  It refuses the first fill that would actually do the
  // damage and lets a legitimate one through, which is why a rig whose only
  // fill is legitimate keeps working.
  //
  // WHY `leveldata.at(0)`, AND IT IS NOT BECAUSE MESH REFINEMENT IS REFUSED.
  // This block used to read: "a slaved cell requires a patch overlap, and B8's
  // PARAMCHECK guard in `CarpetX/src/driver.cxx` refuses mesh refinement on a
  // multi-patch grid, so `n_slaved > 0` implies a single-level hierarchy".
  // AMR-B1 DELETED that refusal and replaced it with two contracts, each
  // evaluated at the site that would violate it, so the sentence named
  // something that no longer exists.  The reason has two legs; they are
  // independent, and THE FIRST IS THE ONE THAT CARRIES THE WEIGHT.
  //
  //   (i) BOTH TERMS DESCRIBE THE DONOR, NOT THE CELL BEING WRITTEN, AND
  //       C-AMR2 PINS THE DONOR TO LEVEL 0.  Term (1) is level-free by
  //       construction: `all_faces_have_symmetries_or_boundaries()` reads
  //       `ghext->patchdata[p].symmetries` and `groupdata.boundaries`
  //       (`CarpetX/src/driver.cxx:1060`), and `get_group_boundaries(gi,
  //       patch)` takes no level (`:368`), so every level's `GroupData`
  //       answers it identically.  Term (2) is the certification of the outer
  //       zone THE DONOR STENCIL WILL READ -- the defect at the top of this
  //       block is an import FROM an unwritten physical outer ghost zone, so
  //       the bit that decides it is the donor patch's, at the level the donor
  //       is on.  C-AMR2 is "no interpatch query point may be answered from a
  //       `level > 0` box", and it is enforced here rather than assumed:
  //       `InterpolationSetup::RefuseAboveLevel0Donors` is a PRE-PASS inside
  //       the `Interpolate` call below, this thorn opts into it
  //       (`/*require_level0_donors=*/true` at the `setup.emplace` above, and
  //       the comment there says why the slaved points are part of its
  //       subject), and it sweeps the whole particle container -- ONE array
  //       carrying this patch's interpatch ghost points and its slaved points
  //       together.  So the level of the cell being WRITTEN never enters this
  //       predicate, and `leveldata.at(0)` is right because the DONOR is
  //       there.
  //
  //  (ii) AND AT EVERY GEOMETRY OF RECORD THERE IS NO SLAVED CELL ABOVE LEVEL
  //       0 TO ASK ABOUT.  `[P217]`: on A2's V1 at `max_num_levels = 2`,
  //       `collect_slaved_interior` classified 35,937 level-1 candidates and
  //       slaved 0 of them, while level 0 carried 11,548 on the cube and
  //       42,602 / 47,458 / 52,133 on the wedge pairs -- a real zero, because
  //       the classifier demonstrably ran there.  `[P313]` re-measured the
  //       same seven numbers from a single `CAPYRX_LOG_LEVELS` line instead of
  //       an 8.5 GB donor flood: `slaved_pl=` is that per-(patch, level)
  //       census and `slaved_maxlevel=` its deepest level.  `[P276]`: the
  //       count is a PER-RANGE quantity, so once this function honours the
  //       caller's active levels the regrid repair's own `L[1,2)` call reports
  //       `nslaved = 0` where it used to report 295,934.
  //
  // LEG (ii) IS A MEASUREMENT AND NOT A THEOREM.  Do not restate it as "C-AMR
  // implies no slaved cell above level 0"; `[P316]` is the counterexample.
  // `evidence/amr/b5/pars/b5_c_p37_r012.par` widens `patch_overlap` to 8 and
  // puts a small OFF-CENTRE refined box inside the resulting overlap band.  It
  // runs to completion; C-AMR reports "holds" with a clearance of 1 coarse
  // cell at both sites it enters; C-AMR2 reports "holds"; and the census reads
  // `slaved_maxlevel=1 slaved_pl=0:1:567` -- all 567 interior vertices of the
  // level-1 box are slaved.  They are the SAME 567 that C-AMR2's line reports
  // as answered from level 0.  That is leg (i) working exactly where leg (ii)
  // does not, and `[P319]` is the same geometry with `boundary = none`, where
  // this guard then REFUSES and names a WEDGE's unwritten face while every
  // slaved cell sits on patch 0's level 1.
  //   Why the band is reachable there and not on the rigs of record: the
  //   overlap band lies OUTSIDE `r0` while the wedges' interpatch ghost zone
  //   lies INSIDE it, so an off-centre box can be in the band without meeting
  //   the zone C-AMR2 protects -- and a wide enough overlap leaves room for
  //   such a box to also stay clear of the patch face, which is C-AMR's
  //   condition.  The arithmetic, and it is a DERIVATION rather than a
  //   measurement: with vertex centring, `ghost_size = g`, ratio 2 and a
  //   prolongation stencil of `s = (prolongation_order + 1) / 2` coarse cells
  //   (`CarpetX/src/prolongate_3d_rf2_impl.hxx:398`, `:1243-1253`), C-AMR
  //   needs a fine nodal top `F <= 2N - 2s - g`, and a refined region is a
  //   union of COARSE cells so `F` is even; hence C-AMR alone excludes a band
  //   cell from level 1 only while `patch_overlap <= s + ceil(g/2)`, which at
  //   the defaults `prolongation_order = 1` and `ghost_size = 3` is
  //   `patch_overlap <= 3`.  `ho_slave.par` has exactly 3 and A2's V1 has 2 --
  //   both inside the bound, production AT it.  Three points corroborate the
  //   derivation (V1's and `[P316]`'s clearances came out at the predicted +5
  //   and +1, and `[P317]` measured the window it predicts, one coarse cell
  //   wide, bounded by C-AMR outward and C-AMR2 inward), and none of them is
  //   a proof.
  //
  // AND THE ONE PARAMETER THAT SUSPENDS BOTH LEGS AT ONCE:
  // `CarpetX::multipatch_amr_contract = "warn"`.  Under it a run reaches here
  // with slaved cells above level 0 AND donors above level 0, and term (2) is
  // then read off a level the donor is not on.  `[P318]` is that state,
  // measured: 92,384 level-1 slaved cells with both contracts reporting a
  // violation and the run exiting 0.  The keyword's own description says it is
  // for diagnosing a violation and not for running with one; this guard is one
  // of the places that depends on that being true.
  //
  // WHERE THE GATE FIRST EVALUATES AT TWO LEVELS, AND THIS PART IS A CODE
  // READING RATHER THAN A MEASUREMENT.  The gate below is `n_slaved > 0`, so
  // on a geometry whose level 1 carries no slaved cell this guard is NOT
  // evaluated at the regrid repair's `L[1,2)` call and is first evaluated at
  // the sync that follows it (`[P276]`).  It still fires before any slaved
  // cell is written, because after AMR-B2 the repair writes none -- that
  // clause is read off the code, not measured.  `[P316]`'s geometry is the
  // other case and it IS measured: there the repair's own `L[1,2)` call
  // reports `nslaved = 567` and the guard is evaluated at it.
  //
  // ONE PROPERTY OF THE CENSUS THAT ITS ABSOLUTE NUMBERS DEPEND ON (`[P321]`,
  // and it is `[P304]` again): `slaved_pl` sums per-COMPONENT list lengths,
  // and a vertex-centred group's boxes OVERLAP on their shared planes, so a
  // slaved vertex on a plane two boxes of one patch share is counted in both.
  // `nslaved` has always had that property and so does `[P276]`'s 295,934.
  // The zero-versus-non-zero reading this block rests on is unaffected; a
  // count compared across two decompositions is not.
  //
  // THE GATE IS `n_slaved > 0`, for B6's reason one block up: a single-patch
  // `Cartesian` patch system reaches this function with no interpatch cell and
  // no slaved cell, and refusing there would break configurations that work.
  // On a multi-rank decomposition the count is this rank's, so a rank that
  // owns no slaved cell stays silent and another rank raises the refusal;
  // `CCTK_VERROR` stops the job either way.
  //
  // WHAT THIS DOES AND DOES NOT BUY.  On this branch the refused configuration
  // does not run silently to completion today either -- with poisoning on it
  // dies in `check_valid_gf` AFTER the interior has already been corrupted,
  // and with poisoning off it dies later at `valid.cxx`'s `error_if_invalid`
  // when a reader asks for the outer boundary that `set_outer` declined to
  // certify.  Both measured (`evidence/fix/b9/before/`).  What changes here is
  // WHEN the run stops and WHAT IT SAYS.  Do not describe this as closing a
  // hole through which wrong answers were leaving; describe it as refusing,
  // before it writes anything, a configuration that cannot be made right.
  if (slave_overlap) {
    std::size_t n_slaved = 0;
    for (const auto &slot : g_interp_cache.slaved_indices) {
      n_slaved += slot.size();
    }

    if (n_slaved > 0) {
      const int npatches_gh = CarpetX::ghext->num_patches();
      for (const auto &varind : varinds) {
        const int gi = CCTK_GroupIndexFromVarI(varind);
        for (int p = 0; p < npatches_gh; ++p) {
          const auto &patchdata = CarpetX::ghext->patchdata.at(p);
          if (patchdata.leveldata.empty()) {
            continue;
          }
          const auto &leveldata = patchdata.leveldata.at(0);
          if (gi < 0 || std::size_t(gi) >= leveldata.groupdata.size()) {
            continue;
          }
          const auto &groupdata = leveldata.groupdata.at(gi);
          if (!groupdata) {
            continue; // not a CCTK_GF group on this grid
          }
          if (groupdata->all_faces_have_symmetries_or_boundaries()) {
            continue; // term (1) is false
          }
          const int vi = int(varind) - groupdata->firstvarindex;
          if (vi < 0 || vi >= groupdata->numvars) {
            continue; // not this group's variable after all
          }
          // BUGFIX_TODO.md B10.  Unreachable: `GroupData`'s constructor
          // sizes `valid` to `group.numtimelevels`
          // (`CarpetX/src/driver.cxx:1013`, ctor at :935), which is at least 1
          // for any declared group.  It used to `continue`, which is the one
          // thing a guard like this may not do -- silently allow a state it
          // does not understand, in the middle of refusing states it does.  If
          // it ever fires, this is not a `GroupData` the driver built and
          // nothing below can be trusted, term (2) least of all.
          if (groupdata->valid.empty()) {
            CCTK_VERROR(
                "MultiPatch1_Interpolate: group %s has an empty `valid` "
                "vector. GroupData's constructor sizes it to the group's time "
                "level count, so this cannot happen; reaching it means this is "
                "not a GroupData the driver built. Refusing rather than "
                "skipping the slave_overlap outer-face check for it.",
                CCTK_FullGroupName(gi));
          }
          if (groupdata->valid.at(0).at(vi).get().valid_outer) {
            continue; // term (2) is false: something wrote the outer zone
          }

          // Name every offending face, by the parameter that would fix it,
          // rather than saying "some face": the user's fix is a specific
          // parameter and there is no reason to make them find out which.
          static_assert(dim == 3, "the face names below assume three "
                                  "directions");
          const char *const dirname[dim] = {"x", "y", "z"};
          std::string faces;
          for (int f = 0; f < 2; ++f) {
            for (int d = 0; d < dim; ++d) {
              if (patchdata.symmetries.at(f).at(d) ==
                      CarpetX::symmetry_t::none &&
                  groupdata->boundaries.at(f).at(d) ==
                      CarpetX::boundary_t::none) {
                if (!faces.empty()) {
                  faces += ", ";
                }
                faces += (f == 0 ? "CarpetX::boundary_"
                                 : "CarpetX::boundary_upper_");
                faces += dirname[d];
              }
            }
          }

          CCTK_VERROR(
              "MultiPatch1_Interpolate: CapyrX_MultiPatch::slave_overlap is "
              "on, and variable %s has a physical outer face that nothing "
              "writes on patch %d: %s is \"none\" there, no per-group override "
              "names %s, and the driver has not marked this variable's outer "
              "region valid, so nothing else wrote it either. Slaved cells are "
              "INTERIOR cells: this call would write %zu of them from an "
              "interpolation whose donor stencil can reach that unwritten "
              "outer ghost zone, and the driver marks the result valid. "
              "Refusing before this fill writes anything. Give that face a "
              "boundary condition, or name %s in one of CarpetX's "
              "{dirichlet,linear_extrapolation,neumann,robin}_{,upper_}"
              "{x,y,z}_vars, or write the group everywhere, or set "
              "CapyrX_MultiPatch::slave_overlap = no.",
              CCTK_FullVarName(int(varind)), p, faces.c_str(),
              CCTK_FullGroupName(gi), n_slaved, CCTK_FullGroupName(gi));
        }
      }
    }
  }

  const std::vector<CCTK_INT> operations(nvars, 0);

  auto &results = g_interp_cache.results;
  results.resize(nvars);
  std::vector<CCTK_REAL *> resultptrs(nvars);

  for (size_t n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }

  g_interp_cache.setup.value().Interpolate(
      cctkGH, nvars, varinds.data(), operations.data(), g_interp_cache.policy,
      resultptrs.data());

// Diagnostic: count NaN values in the interpolated results. Non-zero means
// source data in neighboring patches contains NaN (e.g. their interior cells
// were not yet filled before MultiPatch_Interpolate ran, or
// poison_undefined_values=yes is leaving interior ghosts uninitialized).
#ifdef CCTK_DEBUG
  {
    static bool g_nan_result_reported = false;
    if (!g_nan_result_reported && nvars > 0 && npoints > 0) {
      std::size_t n_nan = 0;
      for (std::size_t n = 0; n < nvars; ++n)
        for (std::size_t p = 0; p < npoints; ++p)
          if (std::isnan(results[n][p])) {
            ++n_nan;
          }
      if (n_nan > 0) {
        g_nan_result_reported = true;
        CCTK_VINFO("MultiPatch1_Interpolate: %zu NaN values in interpolated "
                   "results (out of %zu points x %zu vars = %zu total). Source "
                   "patches may have poisoned interior/ghost cells. ",
                   n_nan, npoints, nvars, npoints * nvars);
      }
    }
  }
#endif // CCTK_DEBUG

  // Step 3: Write back results
  {
    caller_levels.loop_parallel([&](int patch, int level, int index,
                                   int component, const cGH *cctkGH) {
      // AMR-B3: the collection passes assigned this (patch, level) no slot, so
      // `slices.at(location)` below would throw. The skip is the same
      // predicate they used, and the cache key records it so that this call
      // and the build that produced the cache can never disagree.
      if (coords_are_invalid(patch, level))
        return;
      const Loop::GridDescBase grid(cctkGH);
      const std::array<int, dim> centering{0, 0, 0};
      const Loop::GF3D2layout layout(cctkGH, centering);

      std::vector<Loop::GF3D2<CCTK_REAL> > vars;
      vars.reserve(nvars);
      for (const auto &varind : varinds) {
        vars.emplace_back(layout, static_cast<CCTK_REAL *>(
                                      CCTK_VarDataPtrI(cctkGH, 0, varind)));
      }

      const auto &current_patch{g_patch_system->patches.at(patch)};
      const auto &patch_faces{current_patch.faces};

      const Location location{patch, level, index, component};
      const ComponentSlice &slice = g_interp_cache.slices.at(location);

      // Slaved interior cells for this component, computed once at cache
      // rebuild (see the collection pass). Empty unless slave_overlap is on.
      const std::vector<Arith::vect<int, dim> > &slaved =
          g_interp_cache.slaved_indices.at(
              g_interp_cache.component_index.at(location));

// Count ghost cells skipped due to outer boundary (includes pure outer
// cells and corner cells at the outer+interpatch intersection).
// Corner cells are NOT filled here and require a 2nd BC pass in
// SyncGroupsByDirI after this function returns.
#ifdef CCTK_DEBUG
      int n_outer_skipped = 0;
#endif // CCTK_DEBUG

// mp_slave_2.md §4 instrumentation: dump every victim cell's identity
// (Location + grid index + the flat query index `n` fed to
// InterpolationSetup, i.e. `idata(1)` on the CarpetX side) so it can be
// cross-referenced against the CarpetX "DONOR" log lines for the same `n`
// to find which donor cell (and whether it is a genuine ghost zone or an
// interior overlap-band cell) actually fed a given failing point. Logged
// once (var loop index 0), not once per variable. Opt-in via env var.
// INSTRUMENTS (BUGFIX_TODO.md R2 / B10), debug builds only and OFF unless
// `CAPYRX_LOG_DONORS` is set: `GHOSTSKIP` (lo and hi twins) and `VICTIM`
// (buckets `ghost` and `slaved`) below.  ALL KEPT, and `VICTIM bucket=slaved`
// IS LOAD-BEARING: it is half of A9's donor-census join and all of C6's
// slaved-cell identification.  Deleting it in a future cleanup silently removes
// the only record of WHICH interior cells a slaved write replaced.
//
// COST WHEN ON, measured (evidence/fix/b10/before/h_report.txt): per sync and
// per component, one `GHOSTSKIP` line per outer-boundary ghost point and one
// `VICTIM` line per filled cell -- 6084 + 37416 on `color.par`, 11664 + 94656
// on `color_ghost.par`, and `bucket=slaved` adds 6258 on `color_slave.par`.
//
// READING THE STREAM: `OMP_NUM_THREADS=1` AND `MPIEXEC=none`, and they are two
// different preconditions.  These are bare `std::cerr` chains, one `<<` per
// field, emitted from inside an `omp parallel` region: at 16 threads they
// INTERLEAVE, and the damage is invisible to `wc -l` because each thread still
// writes its own newline (`[P26]`, `[P32]`: 19.8M of 21.2M lines malformed at
// 16 threads, 0 at 1).  Separately, at ONE thread, mpiexec's stderr forwarder
// DROPS BYTES from the head of a line under an instrument flood (`[P31]`), which
// no amount of atomicity here would fix.  B10 considered building each line in
// an `ostringstream` and emitting it with a single `<<`; it was rejected,
// because it removes only the first hazard, leaves both operational gates in
// place unchanged, and would make every stream A8 and A9 already measured
// incomparable with the next one.  `tools/run_split.sh` pins both.
#ifdef CCTK_DEBUG
      static const bool log_donors = std::getenv("CAPYRX_LOG_DONORS") != nullptr;
#endif // CCTK_DEBUG

      for (std::size_t n = 0; n < nvars; n++) {
        std::size_t pos = 0;

        // Note: This includes symmetry points
        grid.loop_bnd<0, 0, 0>(grid.nghostzones, [&](const Loop::PointDesc &p) {
          // Skip outer boundaries (pure outer ghost cells and corner cells at
          // the outer+interpatch face intersection). Corner cells must be
          // filled by the 2nd BC pass in SyncGroupsByDirI.
          for (int d = 0; d < dim; ++d) {
            if (p.NI[d] < 0 && patch_faces[0][d].is_outer_boundary) {

#ifdef CCTK_DEBUG
              if (n == 0) {
                ++n_outer_skipped;
                // mp_slave_7.md §5 / mp_slave_8.md instrumentation: log which
                // axis/face triggered the skip for this ghost point, so it
                // can be cross-referenced against bucket (b)'s known cell
                // list to confirm this is the mechanism that makes those
                // cells invisible to this loop.
                if (log_donors) {
                  std::cerr << "GHOSTSKIP patch=" << patch
                            << " level=" << level << " component=" << component
                            << " I=(" << p.I[0] << "," << p.I[1] << ","
                            << p.I[2] << ")"
                            << " axis=" << d << " face=lo\n";
                }
              }
#endif // CCTK_DEBUG

              return;
            }

            if (p.NI[d] > 0 && patch_faces[1][d].is_outer_boundary) {

#ifdef CCTK_DEBUG
              if (n == 0) {
                ++n_outer_skipped;
                if (log_donors) {
                  std::cerr << "GHOSTSKIP patch=" << patch
                            << " level=" << level << " component=" << component
                            << " I=(" << p.I[0] << "," << p.I[1] << ","
                            << p.I[2] << ")"
                            << " axis=" << d << " face=hi\n";
                }
              }
#endif // CCTK_DEBUG

              return;
            }
          }

#ifdef CCTK_DEBUG
          if (n == 0 && log_donors) {
            std::cerr << "VICTIM bucket=ghost patch=" << patch
                      << " level=" << level << " index=" << index
                      << " component=" << component << " I=(" << p.I[0] << ","
                      << p.I[1] << "," << p.I[2] << ")"
                      << " n=" << (slice.offset + pos) << "\n";
          }
#endif // CCTK_DEBUG

          vars[n](p.I) = results[n][slice.offset + pos];

          pos++;
        });

        // Slave overlap band: overwrite non-owned interior cells with the
        // owner's interpolated value, immediately after the ghost cells and
        // in the same order used when their coordinates were collected.
        //
        // These writes used to be gated on an apply_slave_writes argument
        // (mp_slave_2.md §7 Fix #1), because SyncGroupsByDirI ran this
        // function twice per sync and slaving on the first call mutated
        // interior cells that were the second call's own donor source data --
        // the two-pass non-idempotency behind the 696 + 312 corruption.
        // BUGFIX_TODO.md step B3 deleted that first call, so there is exactly
        // one interpolate per sync, no call can read another call's output,
        // and the gate has nothing left to protect against.
        for (const auto &I : slaved) {
#ifdef CCTK_DEBUG
          if (n == 0 && log_donors) {
            std::cerr << "VICTIM bucket=slaved patch=" << patch
                      << " level=" << level << " index=" << index
                      << " component=" << component << " I=(" << I[0] << ","
                      << I[1] << "," << I[2] << ")"
                      << " n=" << (slice.offset + pos) << "\n";
          }
#endif // CCTK_DEBUG

          vars[n](I) = results[n][slice.offset + pos];

          pos++;
        }

        assert(pos == slice.length);
      }

// Report once globally: a non-zero count confirms corner cells exist
// and were left unfilled (outer+interpatch intersection).
#ifdef CCTK_DEBUG
      {
        if (n_outer_skipped > 0 && component == 0) {
          CCTK_VINFO("MultiPatch1_Interpolate [patch %d level %d]: "
                     "skipped %d outer-boundary ghost cells including "
                     "corner cells at outer+interpatch face intersections; "
                     "these are filled by the 2nd BC pass in SyncGroupsByDirI",
                     patch, level, n_outer_skipped);
        }
      }
#endif // CCTK_DEBUG
    });
  }

#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace CapyrX::MultiPatch
