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

Only the cubed-sphere patch system (cartesian cube + 6 wedge shells) is
supported; the owner-patch and outer-boundary logic below is specific to it.

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


def get_owner_patch(x, y, z, r0):
    """Faithful port of CubedSphere::get_owner_patch (cubed_sphere.cxx:23-83)."""
    ax, ay, az = abs(x), abs(y), abs(z)
    if ax <= r0 and ay <= r0 and az <= r0:
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
                 "angular_cells", "radial_cells", "patch_overlap"):
        val = resolve(raw.get("CapyrX_MultiPatch::" + name, ""))
        if val is not None:
            params[name] = val
    return params


def guess_coords_path(color_path):
    p = Path(color_path)
    # "...-color_pre..." must lose the "_pre" too -- a plain substring
    # replace() of "capyrx_testmultipatch-color" would leave it dangling
    # (producing "coordinatesx-vertex_coords_pre...", which never exists),
    # since color_pre's vertex_coords companion is the same file color's is.
    name = re.sub(r"^capyrx_testmultipatch-color(?:_pre)?",
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
        return "\n".join(lines)


def process_file(color_path, coords_path, r0, r1, angular_cells, radial_cells, report,
                 pre_sync=False, patch_overlap=0):
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

        expected_patch = get_owner_patch(vx, vy, vz, r0)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("color_tsv", nargs="+", help="capyrx_testmultipatch-color TSV file(s) to check")
    ap.add_argument("--coords", default=None,
                    help="matching coordinatesx-vertex_coords TSV (auto-detected by default; "
                         "only valid with a single color_tsv argument)")
    ap.add_argument("--inner-boundary", type=float, default=None)
    ap.add_argument("--outer-boundary", type=float, default=None)
    ap.add_argument("--angular-cells", type=int, default=None)
    ap.add_argument("--radial-cells", type=int, default=None)
    ap.add_argument("--patch-overlap", type=int, default=None,
                    help="CapyrX_MultiPatch::patch_overlap (defaults to the value "
                         "found in the parameter file, or 0 if absent -- matches "
                         "the thorn's own default)")
    ap.add_argument("--max-examples", type=int, default=10, help="examples to print per issue kind")
    ap.add_argument("--pre-sync", action="store_true",
                    help="checking a pre-SYNC snapshot (e.g. color_pre): a not-yet-refilled "
                         "intra-patch ghost inside a patch's own interior range is expected, "
                         "not a bug. Omit when checking post-SYNC output.")
    args = ap.parse_args()

    r0, r1 = args.inner_boundary, args.outer_boundary
    angular_cells, radial_cells = args.angular_cells, args.radial_cells
    patch_overlap = args.patch_overlap

    if None in (r0, r1, angular_cells, radial_cells, patch_overlap):
        par_path = find_param_file(args.color_tsv[0])
        if par_path:
            params = parse_par_params(par_path)
            r0 = r0 if r0 is not None else params.get("inner_boundary_radius")
            r1 = r1 if r1 is not None else params.get("outer_boundary_radius")
            angular_cells = angular_cells if angular_cells is not None else params.get("angular_cells")
            radial_cells = radial_cells if radial_cells is not None else params.get("radial_cells")
            patch_overlap = patch_overlap if patch_overlap is not None else params.get("patch_overlap")

    missing = [name for name, v in [
        ("--inner-boundary", r0), ("--outer-boundary", r1),
        ("--angular-cells", angular_cells), ("--radial-cells", radial_cells),
    ] if v is None]
    if missing:
        ap.error(f"could not determine {', '.join(missing)}; pass explicitly")

    angular_cells, radial_cells = int(angular_cells), int(radial_cells)
    patch_overlap = int(patch_overlap) if patch_overlap is not None else 0

    if args.coords and len(args.color_tsv) > 1:
        ap.error("--coords can only be used with a single color_tsv argument")

    report = Report()
    for color_path in args.color_tsv:
        coords_path = Path(args.coords) if args.coords else guess_coords_path(color_path)
        if coords_path is None or not Path(coords_path).exists():
            ap.error(f"could not find a matching coordinatesx-vertex_coords TSV for {color_path}; "
                     f"pass --coords explicitly")
        process_file(color_path, coords_path, r0, r1, angular_cells, radial_cells, report,
                     pre_sync=args.pre_sync, patch_overlap=patch_overlap)

    print(report.summary())
    print()

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
