#!/usr/bin/env python3
"""
Interpret output from the CapyrX_TestMultiPatch coloring test.

The coloring test stamps every grid point with

    color = 10 + patch   (interior points, grid.loop_int_device)
    color = 20 + patch   (boundary points, grid.loop_bnd_device)
    color = 30 + patch   (ghost points, before SYNC overwrites them)
    color = 40 + patch   (overlap band: the patch_overlap-wide slice of a
                          patch's own interior, at faces bordering another
                          patch, that exists purely to give that neighbor
                          enough source data to interpolate its ghost zone)

then runs `SYNC: color`, which should fill every ghost point via CarpetX's
interpatch interpolation from the *valid* (interior/boundary) data of
whichever patch geometrically owns that point.

This script does not trust CarpetX's own patch/region bookkeeping. For every
ghost/exterior point it independently recomputes, from the point's *global*
physical coordinates, which patch should own it (a port of
`CubedSphere::get_owner_patch`, see
repos/CapyrX/CapyrX_MultiPatch/src/cubed_sphere/cubed_sphere.cxx), and
compares that ground truth against the decoded color value.

IMPORTANT: the (x, y, z) columns CarpetX writes into the color TSV itself are
per-patch *local logical* coordinates (computed generically from the AMReX
box geometry in CarpetX/src/io_tsv.cxx, ProbLo + I*CellSize) -- NOT the
multipatch-global physical coordinates. They are not comparable across
patches and must not be fed to get_owner_patch. The real global coordinates
live in the companion `coordinatesx-vertex_coords` TSV (the `vcoordx`,
`vcoordy`, `vcoordz` columns), which this script joins in by
(patch, level, component, i, j, k).

Two patch systems are supported, both a cartesian cube + 6 wedge shells:
cubed_sphere (the cube reuses angular_cells) and Llama (the cube has its own
cartesian_ncells_{i,j,k}). They differ in two ways: how the cube patch's index
range is sized, and the cube/wedge ownership boundary -- cubed_sphere uses a box
(max|coord| <= R), Llama a sphere (r < R, matching real Llama). get_owner_patch
takes a `spherical` flag selected from --patch-system.

Two fields can be checked. The default is the two-digit `color` field (tens =
region, ones = patch), which only works at interpolation_order 0: at order > 0
the interpatch SYNC blends a donor's interior/overlap markers to a non-integer.
With --owner the input is the `owner` field (value = 1 + patch, constant over a
patch), which an order > 0 SYNC returns intact -- use it for the Llama system,
whose Check_Parameters guard forbids order 0.

If a patch is split across multiple AMReX boxes (CarpetX::max_grid_size_*
smaller than the patch's own extent), grid.loop_ghosts_device also fires on
the intra-patch, box-to-box faces (see coloring_1.md) and a point that's
geometrically interior to the *patch* can show up as more than one row here
-- once from the box that owns it as true interior, and once from each
neighboring box's ghost copy of it (distinguished by the `component` column).
Before SYNC, that ghost copy legitimately still holds the not-yet-refilled
ghost marker; pass --pre-sync when checking such a file so this is treated as
expected rather than flagged as a `write_color` bug. Without --pre-sync
(checking post-SYNC output, the default), the same state is a real bug --
the intra-patch box exchange failed to refill it -- and is flagged as such.

Usage:
    check_color.py exe/color/capyrx_testmultipatch-color.it000000.p0000.tsv

    check_color.py --pre-sync exe/color/capyrx_testmultipatch-color_pre.it000000.p0000.tsv

    check_color.py --inner-boundary 0.5 --outer-boundary 2.5 \\
        --angular-cells 8 --radial-cells 16 COLOR.tsv

If --inner-boundary/--outer-boundary/--angular-cells/--radial-cells are not
given, the script recovers them from the parameter file referenced in the
TSV header comment ("# parameter file: ..."). The companion vertex_coords
file is located automatically next to the color file unless --coords is
given explicitly.
"""

import argparse
import math
import re
import sys
from pathlib import Path

CARTESIAN, PLUS_X, MINUS_X, PLUS_Y, MINUS_Y, PLUS_Z, MINUS_Z = range(7)
PATCH_NAMES = {
    CARTESIAN: "cartesian",
    PLUS_X: "plus_x",
    MINUS_X: "minus_x",
    PLUS_Y: "plus_y",
    MINUS_Y: "minus_y",
    PLUS_Z: "plus_z",
    MINUS_Z: "minus_z",
}

INTERIOR, BOUNDARY, GHOST, OVERLAP = 10, 20, 30, 40

# vertex_coords are written with ~17 significant digits; the sqrt()/pow() in
# the cubed-sphere mapping accumulates a little roundoff, so allow a tiny
# relative tolerance when comparing against r0/r1.
REL_TOL = 1e-9

# If more than this fraction of the rows read fail to join against the coords
# TSV, the join itself is broken (coords file mismatch, changed `component`
# numbering, truncated output) and any "no evidence" verdict would be vacuous.
# main() treats crossing this threshold -- or checking zero rows at all -- as a
# hard error (exit 2, distinct from the exit 1 used for a genuine content
# mismatch). See coloring_2 finding 2.
MAX_MISSING_COORDS_FRACTION = 0.5


def get_owner_patch(x, y, z, r0, spherical=False):
    """Port of the C++ get_owner_patch. The wedge selection is identical for both
    patch systems; only the cube region differs:

      * cubed_sphere (spherical=False): box, max(|x|,|y|,|z|) <= R -> cube
        (CubedSphere::get_owner_patch, cubed_sphere.cxx).
      * Llama (spherical=True): sphere, x^2+y^2+z^2 < R^2 -> cube, matching real
        Llama's global_to_local_Thornburg04 (thornburg04.cc) and CapyrX
        Llama::get_owner_patch (llama/llama.cxx).
    """
    ax, ay, az = abs(x), abs(y), abs(z)
    if spherical:
        if ax * ax + ay * ay + az * az < r0 * r0:
            return CARTESIAN
    elif ax <= r0 and ay <= r0 and az <= r0:
        return CARTESIAN
    coords = [ax, ay, az]
    max_idx = coords.index(max(coords))  # first-max wins ties, matches std::max_element
    if max_idx == 0:
        return PLUS_X if x > 0.0 else MINUS_X
    if max_idx == 1:
        return PLUS_Y if y > 0.0 else MINUS_Y
    return PLUS_Z if z > 0.0 else MINUS_Z


def in_overlap_band(patch, i, j, k, angular_cells, radial_cells, patch_overlap):
    """True if (i, j, k) falls in the patch_overlap-wide band a patch's own
    (possibly already-widened) index range grows by at each face that
    borders another patch -- write_color's overlap_marker region -- as
    opposed to the patch's "true" (non-overlap) interior.

    Mirrors cubed_sphere.cxx's make_patch: wedges widen both ends of both
    angular axes and only the cartesian-facing (inner) end of the radial
    axis (its outer end is the true Dirichlet boundary, never widened); the
    cartesian patch widens both ends of all three axes.
    """
    if patch_overlap == 0:
        return False
    angular_lo, angular_hi = patch_overlap, angular_cells + patch_overlap
    in_angular_core = angular_lo <= i <= angular_hi and angular_lo <= j <= angular_hi
    if patch == CARTESIAN:
        return not (in_angular_core and angular_lo <= k <= angular_hi)
    return not (in_angular_core and patch_overlap <= k)


def find_param_file(tsv_path):
    with open(tsv_path) as f:
        for _ in range(5):
            line = f.readline()
            m = re.search(r'parameter file:\s*"([^"]+)"', line)
            if m:
                return m.group(1)
    return None


def parse_par_params(par_path):
    """Extract CapyrX_MultiPatch::{inner_boundary_radius,outer_boundary_radius,
    angular_cells,radial_cells,patch_overlap} from a par file.

    Values may be numeric literals or a cross-parameter reference like
    `CapyrX_MultiPatch::patch_overlap = CarpetX::ghost_size` (precedented in
    repos/CapyrX/CapyrX_MultiPatch/par/canudax_bench.par) -- resolved against
    every other `Thorn::param = value` assignment in the same file.
    """
    try:
        text = Path(par_path).read_text()
    except OSError:
        return {}

    raw = {}
    pattern = re.compile(
        r"([A-Za-z_][A-Za-z0-9_]*::[A-Za-z_][A-Za-z0-9_]*)\s*=\s*"
        r"([A-Za-z_][A-Za-z0-9_:]*|[0-9.eE+-]+)"
    )
    for m in pattern.finditer(text):
        raw[m.group(1)] = m.group(2)

    def resolve(value, seen=()):
        try:
            return float(value)
        except ValueError:
            pass
        if value in raw and value not in seen:
            return resolve(raw[value], seen + (value,))
        return None

    params = {}
    for name in ("inner_boundary_radius", "outer_boundary_radius",
                 "angular_cells", "radial_cells", "patch_overlap",
                 "cartesian_ncells_i", "cartesian_ncells_j", "cartesian_ncells_k"):
        val = resolve(raw.get("CapyrX_MultiPatch::" + name, ""))
        if val is not None:
            params[name] = val
    return params


def parse_patch_system(par_path):
    """Read CapyrX_MultiPatch::patch_system (a quoted keyword, possibly with a
    space, e.g. "Cubed sphere" or "Llama"). Returns "llama" for the Llama system
    and "cubed_sphere" otherwise (the two the owner/color tests distinguish)."""
    try:
        text = Path(par_path).read_text()
    except OSError:
        return None
    m = re.search(r'CapyrX_MultiPatch::patch_system\s*=\s*"([^"]*)"', text)
    if m is None:
        return None
    return "llama" if m.group(1).strip().lower() == "llama" else "cubed_sphere"


def guess_coords_path(color_path):
    p = Path(color_path)
    # "...-color_pre..." must lose the "_pre" too -- a plain substring
    # replace() of "capyrx_testmultipatch-color" would leave it dangling
    # (producing "coordinatesx-vertex_coords_pre...", which never exists),
    # since color_pre's vertex_coords companion is the same file color's is.
    # The owner field shares the same vertex_coords companion.
    name = re.sub(r"^capyrx_testmultipatch-(?:color(?:_pre)?|owner)",
                  "coordinatesx-vertex_coords", p.name)
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def read_tsv_rows(path):
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            yield line.rstrip("\n").split("\t")


def load_global_coords(coords_path):
    """key (patch, level, component, i, j, k) -> (vcoordx, vcoordy, vcoordz)."""
    table = {}
    for fields in read_tsv_rows(coords_path):
        if len(fields) < 14:
            continue
        _it, _t, patch, level, comp, i, j, k, _x, _y, _z, vx, vy, vz = fields[:14]
        key = (int(patch), int(level), int(comp), int(i), int(j), int(k))
        table[key] = (float(vx), float(vy), float(vz))
    return table


def decode_color(color):
    if not math.isfinite(color):
        return None, None  # NaN/inf poison (poison_undefined_values) -> NON-INTEGER-COLOR
    c = round(color)
    if abs(c - color) > 1e-6:
        return None, None  # not an integer marker at all
    if c == 0:
        return 0, None
    return c - (c % 10), c % 10


class Report:
    def __init__(self):
        self.total = 0
        self.checked = 0
        self.passed = 0
        self.missing_coords = 0
        self.by_kind = {}

    def flag(self, kind, row):
        self.by_kind.setdefault(kind, []).append(row)

    def summary(self):
        lines = [
            f"Total rows read:              {self.total}",
            f"Rows missing a coords match:  {self.missing_coords}",
            f"Rows checked against ground truth: {self.checked}",
            f"Passed:                       {self.passed}",
        ]
        if self.total and (self.checked == 0 or
                           self.missing_coords > MAX_MISSING_COORDS_FRACTION * self.total):
            lines.append("  ^^ coords join looks broken (see RESULT below): with these "
                         "checked / missing-coords counts the verdict is not trustworthy.")
        return "\n".join(lines)


def process_file(color_path, coords_path, r0, r1, angular_cells, radial_cells, report,
                 pre_sync=False, patch_overlap=0, spherical=False):
    coords = load_global_coords(coords_path)

    # With patch_overlap > 0 (repos/CapyrX/CapyrX_MultiPatch/src/cubed_sphere/
    # cubed_sphere.cxx:1074-1127, make_patch), each patch's own stored index
    # range grows: wedge patches gain `overlap` extra cells on *both* angular
    # axes' ends and on the inner (cartesian-facing) end of the radial axis
    # only (the outer end is the true Dirichlet boundary, never overlapped);
    # the cartesian patch gains `overlap` extra cells on all three axes (it
    # borders wedges on all six faces). The overlap band is functionally
    # interior to whichever patch(es) now store it -- each patch's own rows
    # there must still show that patch's own marker, so own_interior below
    # must use these widened extents, not the bare non-overlap cell counts.
    angular_extent = angular_cells + 2 * patch_overlap
    radial_extent = radial_cells + patch_overlap

    for fields in read_tsv_rows(color_path):
        if len(fields) < 12:
            continue
        _it, _t, patch_s, level_s, comp_s, i_s, j_s, k_s, _x, _y, _z, color_s = fields[:12]
        patch, level, comp = int(patch_s), int(level_s), int(comp_s)
        i, j, k = int(i_s), int(j_s), int(k_s)
        color = float(color_s)
        report.total += 1

        key = (patch, level, comp, i, j, k)
        if key not in coords:
            report.missing_coords += 1
            continue
        vx, vy, vz = coords[key]

        category, stored_patch = decode_color(color)
        if category is None:
            report.checked += 1
            report.flag("NON-INTEGER-COLOR (poison value never overwritten?)",
                         (patch, i, j, k, vx, vy, vz, color))
            continue

        if patch == CARTESIAN:
            own_interior = 0 <= i <= angular_extent and 0 <= j <= angular_extent and 0 <= k <= angular_extent
        else:
            own_interior = 0 <= i <= angular_extent and 0 <= j <= angular_extent and 0 <= k <= radial_extent

        report.checked += 1

        if own_interior:
            # Geometrically interior to the *patch*. If the patch is a
            # single AMReX box this row can only be that box's own true
            # interior, written directly by loop_int_device -- the only
            # possible correct value is this patch's own marker (INTERIOR,
            # or OVERLAP if this point also falls in the patch_overlap band
            # write_color carves out near faces bordering another patch),
            # independent of any geometric tie-break at shared-interface
            # vertices.
            #
            # If the patch is split across multiple boxes, this same
            # (patch, i, j, k) can also show up as a *different* box's ghost
            # copy of that point (component differs) -- ground truth for
            # that row's own component is still "this patch's own marker",
            # but only once the intra-patch box exchange has run. Pre-SYNC,
            # a ghost marker here is the expected, not-yet-refilled state,
            # not a write_color bug; post-SYNC it means the intra-patch
            # exchange silently failed to refill it.
            overlap = in_overlap_band(patch, i, j, k, angular_cells, radial_cells, patch_overlap)
            expected_category = OVERLAP if overlap else INTERIOR
            mismatch_kind = "OWN-OVERLAP-MISMATCH" if overlap else "OWN-INTERIOR-MISMATCH"
            if category == expected_category and stored_patch == patch:
                report.passed += 1
            elif category == GHOST and stored_patch == patch and pre_sync:
                report.flag("informational: intra-patch ghost not yet filled (expected pre-SYNC)",
                             (patch, i, j, k, vx, vy, vz, color))
            elif category == GHOST and stored_patch == patch:
                report.flag("INTRA-PATCH-GHOST-NOT-SYNCED (box exchange never refilled this point)",
                             (patch, i, j, k, vx, vy, vz, color))
            else:
                report.flag(f"{mismatch_kind} (write_color itself is wrong here)",
                             (patch, i, j, k, vx, vy, vz, color, "expected category", expected_category))
            continue

        # Ghost point (of this patch). Ground truth comes from where its
        # *global* coordinate actually sits, independent of CarpetX's own
        # bookkeeping.
        if pre_sync and category == BOUNDARY and stored_patch == patch:
            # loop_bnd_device fires on every bbox=true face of this patch --
            # an interpatch border exactly as much as a true physical/
            # Dirichlet edge (bbox carries no such distinction, see the
            # module docstring) -- so pre-SYNC, *every* not-yet-filled ghost
            # or exterior point here still holds write_color's own
            # never-overwritten boundary marker, whether SYNC is about to
            # interpatch-interpolate it or zero it as a physical BC point.
            # That's the expected pre-SYNC state, not a bug.
            report.flag("informational: not yet filled by SYNC (expected pre-SYNC)",
                         (patch, i, j, k, vx, vy, vz, color))
            continue

        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        inside_cube = max(abs(vx), abs(vy), abs(vz)) <= r0 * (1 + REL_TOL)
        inside_domain = inside_cube or r <= r1 * (1 + REL_TOL)

        if not inside_domain:
            # Genuinely outside the whole multipatch domain: filled by
            # CarpetX's physical (Dirichlet) BC, not interpatch
            # interpolation. `color` has no dirichlet_values tag, so 0 is
            # the only correct fill.
            if color != 0.0:
                report.flag("EXTERIOR-NOT-ZEROED (physical BC did not overwrite this point)",
                             (patch, i, j, k, vx, vy, vz, color))
            else:
                report.passed += 1
            continue

        expected_patch = get_owner_patch(vx, vy, vz, r0, spherical=spherical)
        out_of_range = sum([
            not (0 <= i <= angular_extent),
            not (0 <= j <= angular_extent),
            not (0 <= k <= (angular_extent if patch == CARTESIAN else radial_extent)),
        ])
        corner_tag = " [corner/edge: %d axes out of range]" % out_of_range if out_of_range > 1 else ""

        if category == GHOST:
            report.flag("GHOST-MARKER-LEAK (sourced from a not-yet-filled ghost cell)" + corner_tag,
                         (patch, i, j, k, vx, vy, vz, color, "expected patch", expected_patch, PATCH_NAMES[expected_patch]))
            continue

        if category == 0:
            report.flag("DEFAULT-LEAK (Dirichlet-zero leaked inside the valid domain)" + corner_tag,
                         (patch, i, j, k, vx, vy, vz, color, "expected patch", expected_patch, PATCH_NAMES[expected_patch]))
            continue

        if stored_patch != expected_patch:
            report.flag("WRONG-NEIGHBOR (interpolated from the wrong patch)" + corner_tag,
                         (patch, i, j, k, vx, vy, vz, color, "expected patch", expected_patch, PATCH_NAMES[expected_patch]))
            continue

        if category == BOUNDARY:
            report.flag("boundary-marker-survived (informational -- not observed to occur; investigate if it starts)",
                         (patch, i, j, k, vx, vy, vz, color))
            continue

        # category is INTERIOR or OVERLAP here, either is a pass: whether an
        # interpolation stencil for this ghost point actually reaches into
        # the donor's overlap band (vs. its plain interior) depends on
        # CarpetX::interpolation_order (see required_overlap in
        # CapyrX_MultiPatch_Check_Parameters, multipatch.cxx) -- e.g. with
        # order 0 (all current par files) it never does, so OVERLAP never
        # shows up here even when patch_overlap > 0. Only stored_patch
        # (checked above) is meaningful ground truth for this row.
        report.passed += 1


def process_file_owner(owner_path, coords_path, r0, r1, angular_cells,
                       radial_cells, cube_ncells, report, patch_overlap=0,
                       spherical=False):
    """Check the ownership-marker field (interface.ccl `owner`).

    Each valid cell holds `1 + patch` (0 reserved for unfilled/exterior), a
    constant across a whole patch, so an order>0 interpatch SYNC returns a ghost
    its single donor patch's exact integer index -- unlike the two-digit color
    markers, whose interior/overlap digits blend to a non-integer at order>0.
    This is the Llama-capable ownership check (the color test only works at
    interpolation_order 0, which the Llama Check_Parameters guard forbids).

    Ground truth for a ghost is get_owner_patch of its *global* coordinate; for
    an own-interior cell it is the cell's own patch. The cube patch's index range
    is sized from cube_ncells (Llama's independent h_cartesian), not angular_cells
    -- the one geometric difference from cubed_sphere here.
    """
    coords = load_global_coords(coords_path)

    angular_extent = angular_cells + 2 * patch_overlap
    radial_extent = radial_cells + patch_overlap
    cube_extent = tuple(n + 2 * patch_overlap for n in cube_ncells)

    for fields in read_tsv_rows(owner_path):
        if len(fields) < 12:
            continue
        _it, _t, patch_s, level_s, comp_s, i_s, j_s, k_s, _x, _y, _z, val_s = fields[:12]
        patch, level, comp = int(patch_s), int(level_s), int(comp_s)
        i, j, k = int(i_s), int(j_s), int(k_s)
        value = float(val_s)
        report.total += 1

        key = (patch, level, comp, i, j, k)
        if key not in coords:
            report.missing_coords += 1
            continue
        vx, vy, vz = coords[key]

        report.checked += 1

        hi_i = cube_extent[0] if patch == CARTESIAN else angular_extent
        hi_j = cube_extent[1] if patch == CARTESIAN else angular_extent
        hi_k = cube_extent[2] if patch == CARTESIAN else radial_extent
        own_interior = 0 <= i <= hi_i and 0 <= j <= hi_j and 0 <= k <= hi_k
        out_of_range = sum([not (0 <= i <= hi_i), not (0 <= j <= hi_j),
                            not (0 <= k <= hi_k)])
        # A ghost with >=2 axes out of range sits at a patch edge/corner, where
        # the point is covered by two or three donor patches at once. An order>N
        # centered stencil there cannot stay inside one donor, so its owner value
        # legitimately blends across patches (or lands on the "wrong" one). The
        # test's real assertion is the single-donor face ghosts (<=1 axis out) and
        # the valid cells; overset corner/edge ambiguity is tolerated as
        # informational. See Step 8 / R2 (face-centre donor availability).
        is_corner = (not own_interior) and out_of_range >= 2
        corner_tag = " [corner/edge: %d axes out of range]" % out_of_range if out_of_range > 1 else ""

        if not math.isfinite(value):
            if is_corner:
                report.flag("informational: corner/edge ghost non-finite (overset, expected)"
                            + corner_tag, (patch, i, j, k, vx, vy, vz, value))
            else:
                report.flag("NON-FINITE-OWNER (poison value never overwritten?)",
                             (patch, i, j, k, vx, vy, vz, value))
            continue
        pv = round(value)
        if abs(pv - value) > 1e-6:
            # A single donor patch is a constant field, so a face ghost must come
            # back an exact integer; a blend there means the stencil spanned more
            # than one source value -- a real defect. At a corner/edge it is the
            # expected overset ambiguity.
            if is_corner:
                report.flag("informational: corner/edge ghost blends across patches (overset, expected)"
                            + corner_tag, (patch, i, j, k, vx, vy, vz, value))
            else:
                report.flag("NON-INTEGER-OWNER (interpolation blended across patches)",
                             (patch, i, j, k, vx, vy, vz, value))
            continue

        if own_interior:
            if pv == 1 + patch:
                report.passed += 1
            else:
                report.flag("OWN-OWNER-MISMATCH (write_owner wrong or intra-patch "
                            "box exchange failed here)",
                             (patch, i, j, k, vx, vy, vz, value, "expected", 1 + patch))
            continue

        # Ghost/exterior point: ground truth is where its global coordinate sits.
        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        inside_cube = max(abs(vx), abs(vy), abs(vz)) <= r0 * (1 + REL_TOL)
        inside_domain = inside_cube or r <= r1 * (1 + REL_TOL)

        if not inside_domain:
            # Outside the whole domain: filled by the physical (Dirichlet) BC.
            # `owner` has no dirichlet_values tag, so 0 is the only correct fill.
            if pv == 0:
                report.passed += 1
            else:
                report.flag("EXTERIOR-NOT-ZEROED (physical BC did not overwrite this point)",
                             (patch, i, j, k, vx, vy, vz, value))
            continue

        expected_patch = get_owner_patch(vx, vy, vz, r0, spherical=spherical)
        if pv == 1 + expected_patch:
            report.passed += 1
        elif is_corner:
            report.flag("informational: corner/edge ghost from an adjacent donor (overset, expected)"
                        + corner_tag, (patch, i, j, k, vx, vy, vz, value,
                         "expected patch", expected_patch, PATCH_NAMES[expected_patch],
                         "got patch", pv - 1))
        elif pv == 0:
            report.flag("DEFAULT-LEAK (Dirichlet-zero leaked inside the valid domain)",
                         (patch, i, j, k, vx, vy, vz, value,
                          "expected patch", expected_patch, PATCH_NAMES[expected_patch]))
        else:
            report.flag("WRONG-OWNER (interpolated from the wrong patch)",
                         (patch, i, j, k, vx, vy, vz, value,
                          "expected patch", expected_patch, PATCH_NAMES[expected_patch],
                          "got patch", pv - 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("color_tsv", nargs="+",
                    help="capyrx_testmultipatch-color (or -owner, with --owner) TSV file(s)")
    ap.add_argument("--coords", default=None,
                    help="matching coordinatesx-vertex_coords TSV (auto-detected by default; "
                         "only valid with a single color_tsv argument)")
    ap.add_argument("--owner", action="store_true",
                    help="the input is the ownership-marker field (interface.ccl `owner`, "
                         "value = 1 + patch), not the two-digit color field. Use this for any "
                         "run with interpolation_order > 0 (e.g. the Llama system, whose "
                         "Check_Parameters guard forbids order 0): the color markers blend to a "
                         "non-integer at order > 0, the owner marker does not.")
    ap.add_argument("--patch-system", choices=["cubed_sphere", "llama"], default=None,
                    help="which patch system produced the TSV (defaults to the value found in "
                         "the parameter file, or cubed_sphere). Only affects how the cube "
                         "patch's index range is sized: cubed_sphere reuses angular_cells, "
                         "Llama uses its independent cartesian_ncells_{i,j,k}.")
    ap.add_argument("--inner-boundary", type=float, default=None)
    ap.add_argument("--outer-boundary", type=float, default=None)
    ap.add_argument("--angular-cells", type=int, default=None)
    ap.add_argument("--radial-cells", type=int, default=None)
    ap.add_argument("--cube-ncells-i", type=int, default=None,
                    help="Llama cube resolution (CapyrX_MultiPatch::cartesian_ncells_i); "
                         "defaults to the par-file value, or angular_cells for cubed_sphere")
    ap.add_argument("--cube-ncells-j", type=int, default=None)
    ap.add_argument("--cube-ncells-k", type=int, default=None)
    ap.add_argument("--patch-overlap", type=int, default=None,
                    help="CapyrX_MultiPatch::patch_overlap (defaults to the value "
                         "found in the parameter file, or 0 if absent -- matches "
                         "the thorn's own default)")
    ap.add_argument("--max-examples", type=int, default=10, help="examples to print per issue kind")
    ap.add_argument("--pre-sync", action="store_true",
                    help="checking a pre-SYNC snapshot (e.g. color_pre): a not-yet-refilled "
                         "intra-patch ghost inside a patch's own interior range is expected, "
                         "not a bug. Omit when checking post-SYNC output. (color field only.)")
    args = ap.parse_args()

    if args.owner and args.pre_sync:
        ap.error("--pre-sync is a color-field notion; the owner check is post-SYNC only")

    r0, r1 = args.inner_boundary, args.outer_boundary
    angular_cells, radial_cells = args.angular_cells, args.radial_cells
    patch_overlap = args.patch_overlap
    patch_system = args.patch_system
    cube_ni, cube_nj, cube_nk = args.cube_ncells_i, args.cube_ncells_j, args.cube_ncells_k

    needs_par = (None in (r0, r1, angular_cells, radial_cells, patch_overlap)
                 or patch_system is None
                 or None in (cube_ni, cube_nj, cube_nk))
    if needs_par:
        par_path = find_param_file(args.color_tsv[0])
        if par_path:
            params = parse_par_params(par_path)
            r0 = r0 if r0 is not None else params.get("inner_boundary_radius")
            r1 = r1 if r1 is not None else params.get("outer_boundary_radius")
            angular_cells = angular_cells if angular_cells is not None else params.get("angular_cells")
            radial_cells = radial_cells if radial_cells is not None else params.get("radial_cells")
            patch_overlap = patch_overlap if patch_overlap is not None else params.get("patch_overlap")
            patch_system = patch_system if patch_system is not None else parse_patch_system(par_path)
            cube_ni = cube_ni if cube_ni is not None else params.get("cartesian_ncells_i")
            cube_nj = cube_nj if cube_nj is not None else params.get("cartesian_ncells_j")
            cube_nk = cube_nk if cube_nk is not None else params.get("cartesian_ncells_k")

    missing = [name for name, v in [
        ("--inner-boundary", r0), ("--outer-boundary", r1),
        ("--angular-cells", angular_cells), ("--radial-cells", radial_cells),
    ] if v is None]
    if missing:
        ap.error(f"could not determine {', '.join(missing)}; pass explicitly")

    angular_cells, radial_cells = int(angular_cells), int(radial_cells)
    patch_overlap = int(patch_overlap) if patch_overlap is not None else 0
    if patch_system is None:
        patch_system = "cubed_sphere"

    # The cube's index range: cubed_sphere reuses angular_cells; Llama uses its
    # own cartesian_ncells. Fall back to angular_cells for any axis not found.
    if patch_system == "llama":
        cube_ncells = (int(cube_ni) if cube_ni is not None else angular_cells,
                       int(cube_nj) if cube_nj is not None else angular_cells,
                       int(cube_nk) if cube_nk is not None else angular_cells)
    else:
        cube_ncells = (angular_cells, angular_cells, angular_cells)

    if args.coords and len(args.color_tsv) > 1:
        ap.error("--coords can only be used with a single color_tsv argument")

    # Llama splits cube/wedge ownership at the sphere r=R; cubed_sphere uses the
    # box. The ground-truth get_owner_patch must match the system that wrote the
    # data, or every corner-shell ghost is misjudged.
    spherical = (patch_system == "llama")

    report = Report()
    for color_path in args.color_tsv:
        coords_path = Path(args.coords) if args.coords else guess_coords_path(color_path)
        if coords_path is None or not Path(coords_path).exists():
            ap.error(f"could not find a matching coordinatesx-vertex_coords TSV for {color_path}; "
                     f"pass --coords explicitly")
        if args.owner:
            process_file_owner(color_path, coords_path, r0, r1, angular_cells,
                               radial_cells, cube_ncells, report, patch_overlap=patch_overlap,
                               spherical=spherical)
        else:
            process_file(color_path, coords_path, r0, r1, angular_cells, radial_cells, report,
                         pre_sync=args.pre_sync, patch_overlap=patch_overlap, spherical=spherical)

    print(report.summary())
    print()

    # Vacuous-pass guard (coloring_2 finding 2). The any_critical verdict below
    # inspects only report.by_kind, so a run that validated *nothing* -- because
    # the coords join produced no matches -- would otherwise print "no evidence"
    # and exit 0 while having checked zero rows. A broken/empty join is a
    # "checker couldn't run" condition, not a clean result: fail with exit 2
    # (distinct from the exit 1 used for a genuine content mismatch, so CI can
    # tell the two apart), with an explicit message separate from any content
    # mismatch so the failure mode is unambiguous.
    if report.checked == 0:
        print(f"RESULT: nothing was validated -- 0 of {report.total} rows read were "
              "checked against ground truth. The coords join produced no matches "
              "(wrong/missing coords TSV, mismatched component numbering, or "
              "truncated output); this is not a clean pass.", file=sys.stderr)
        sys.exit(2)
    if report.total and report.missing_coords > MAX_MISSING_COORDS_FRACTION * report.total:
        pct = 100.0 * report.missing_coords / report.total
        print(f"RESULT: coords join largely broke -- {report.missing_coords} of "
              f"{report.total} rows ({pct:.1f}%) had no coords match, over the "
              f"{MAX_MISSING_COORDS_FRACTION:.0%} threshold. Only {report.checked} rows "
              "were actually checked; the verdict is not trustworthy. Treating as a "
              "hard failure.", file=sys.stderr)
        sys.exit(2)

    any_critical = False
    for kind, items in sorted(report.by_kind.items(), key=lambda kv: -len(kv[1])):
        critical = "informational" not in kind
        any_critical = any_critical or (critical and items)
        print(f"=== {kind} ({len(items)}) ===")
        for row in items[: args.max_examples]:
            print("  ", row)
        if len(items) > args.max_examples:
            print(f"   ... and {len(items) - args.max_examples} more")
        print()

    if any_critical:
        print("RESULT: interpatch interpolation is producing incorrect ghost data.")
        sys.exit(1)
    else:
        print("RESULT: no evidence of incorrect interpatch interpolation found.")
        sys.exit(0)


if __name__ == "__main__":
    main()
