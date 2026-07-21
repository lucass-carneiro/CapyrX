#!/usr/bin/env python3
"""
Interpret output from the CapyrX_TestMultiPatch physical-outer-boundary
NaN-injection test (see mp_nan_test_design.md).

`CapyrX_TestMultiPatch_write_nan_test` fills `nan_test = 1.0` everywhere,
except at cells selected by `nan_injection_mode` ("ghost" / "ghost_and_interior"
/ "interior") relative to each patch's *genuine physical-outer* face (as
opposed to an interpatch face -- the same distinction `write_color` already
draws via `MultiPatch_GetBoundarySpecification2`), where it writes NaN
instead. `nan_test_pre` is a snapshot taken immediately afterwards, before
`SYNC: nan_test` runs -- an exact "what did we inject" ground truth. `SYNC`
then runs CarpetX's full boundary-condition + `MultiPatch_Interpolate`
pipeline, using whichever `CarpetX::boundary_*` kind the `.par` file selects.

This script classifies every grid point by comparing `nan_test_pre` (pre-sync)
against `nan_test` (post-sync):

    pre    post          meaning
    1.0    1.0           untouched -- expected for the bulk of the domain
    NaN    NaN           corruption survives at the injection site itself
    NaN    finite        BC/interpolation cleaned the injected corruption
    1.0    NaN           LEAK -- a clean cell came out corrupted
    1.0    finite!=1.0   unexpected-value leak (BC fill diverged from the
                         sentinel, independent of the NaN mechanism)

Every `1.0 -> NaN` leak is further classified by *where* it sits, reusing
`check_color.py`'s `get_owner_patch` port and its own-patch-interior index
arithmetic, plus explicit cubed-sphere face geometry (`cubed_sphere.cxx`'s
`make_patch`: each wedge patch has exactly one genuine physical-outer face,
the high side of its radial (k) axis; every other face -- both angular axes,
lo and hi, and the low/cartesian-facing side of the radial axis -- is
interpatch; the cartesian patch has no physical-outer face at all, since it
is bordered by wedges on all six faces and is therefore never injected):

    LEAK-LOCAL        an untargeted own-patch cell next to an injected one.
                      Lowest severity.
    LEAK-CORNER       a cell simultaneously near the physical-outer face and
                      near an interpatch face -- the "corner cell
                      catastrophe" boundaries_impl.hxx's own comments name --
                      AND whose true owner (get_owner_patch) is itself (not a
                      genuine ownership mismatch). Tracked but not
                      exit-critical unless --strict.
    LEAK-INTERIOR     a deep-interior cell (more than --band cells from every
                      face) turned NaN. Always critical.
    LEAK-CROSS-PATCH  a cell whose geometric owner (get_owner_patch, from its
                      global coordinates) is a *different* patch than the one
                      reporting the leak -- corruption crossed an interpatch
                      seam via MultiPatch_Interpolate's donor lookup. Always
                      critical. Checked *before* LEAK-CORNER's geometric
                      test, not after: a cell can simultaneously be near the
                      physical-outer face and genuinely cross-owned, and an
                      ownership mismatch is always the more serious finding
                      -- see mp_nan_test_10.md for the false negative this
                      ordering fixes (checking the geometric corner test
                      first let a whole order's worth of genuine
                      cross-patch leak hide inside LEAK-CORNER, which
                      `--strict` is required to catch).

`UNEXPECTED-VALUE-LEAK` (the `1.0 -> finite!=1.0` row) is reported flat, not
sub-classified -- mirrors check_color.py's DEFAULT-LEAK: always critical,
regardless of where it sits.

Usage:
    check_nan_test.py exe/nan_dirichlet_ghost/capyrx_testmultipatch-nan_test.it000000.p0000.tsv

    check_nan_test.py --strict COLOR.tsv   # also fail on LEAK-CORNER

    check_nan_test.py --inner-boundary 0.5 --outer-boundary 2.5 \\
        --angular-cells 8 --radial-cells 16 --band 2 NAN_TEST.tsv

If --inner-boundary/--outer-boundary/--angular-cells/--radial-cells/--band are
not given, the script recovers them from the parameter file referenced in the
TSV header comment ("# parameter file: ..."), the same way check_color.py
does (--band from `CarpetX::ghost_size`, defaulting to 2 if neither is found).
The companion `nan_test_pre` and `coordinatesx-vertex_coords` TSVs are located
automatically next to the `nan_test` file unless --pre/--coords are given
(only valid with a single positional argument).
"""

import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_color import (  # noqa: E402
    CARTESIAN,
    PATCH_NAMES,
    REL_TOL,
    find_param_file,
    get_owner_patch,
    load_global_coords,
    parse_par_params,
    read_tsv_rows,
)

SENTINEL = 1.0
# nan_test's clean value is written bit-exact as 1.0; nan_test_pre's snapshot
# is a plain copy, so pre-sync comparisons need no tolerance. Post-sync values
# have gone through interpolation/BC arithmetic, so a small absolute
# tolerance avoids flagging harmless floating-point roundoff as a leak.
SENTINEL_ATOL = 1e-8


def guess_pre_path(nan_test_path):
    p = Path(nan_test_path)
    if "-nan_test_pre." in p.name:
        return None  # caller passed the pre-sync file itself, not post-sync
    name = p.name.replace("-nan_test.", "-nan_test_pre.")
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def guess_coords_path(nan_test_path):
    p = Path(nan_test_path)
    name = re.sub(r"^capyrx_testmultipatch-nan_test(?:_pre)?",
                  "coordinatesx-vertex_coords", p.name)
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def parse_ghost_size(par_path):
    """CarpetX::ghost_size from the .par file, if set explicitly (not the
    -1 "use ghost_size_[xyz]" sentinel, which this script does not resolve
    further)."""
    try:
        text = Path(par_path).read_text()
    except OSError:
        return None
    m = re.search(r"CarpetX::ghost_size\s*=\s*([0-9]+)", text)
    return int(m.group(1)) if m else None


def load_value_table(path):
    """key (patch, level, component, i, j, k) -> nan_test(_pre) value."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) < 12:
            continue
        _it, _t, patch, level, comp, i, j, k, _x, _y, _z, value = fields[:12]
        key = (int(patch), int(level), int(comp), int(i), int(j), int(k))
        table[key] = float(value)
    return table


def near_lo(idx, band):
    return idx <= band - 1


def near_hi(idx, extent, band):
    return idx >= extent - (band - 1)


def is_own_interior(patch, i, j, k, angular_extent, radial_extent):
    if patch == CARTESIAN:
        return 0 <= i <= angular_extent and 0 <= j <= angular_extent and 0 <= k <= angular_extent
    return 0 <= i <= angular_extent and 0 <= j <= angular_extent and 0 <= k <= radial_extent


def classify_leak(patch, i, j, k, vx, vy, vz, r0, r1,
                  angular_extent, radial_extent, band):
    """Sub-classify a 1.0 -> NaN leak row. Returns (subkind, expected_patch)
    where expected_patch is None unless subkind is CROSS-PATCH (or a
    same-patch ground-truth check happened to run anyway).

    Ownership is checked *before* the geometric CORNER test (mp_nan_test_10.md):
    a ghost/exterior point whose true owner (get_owner_patch) differs from the
    reporting patch is always CROSS-PATCH, even when it also sits near the
    patch's own physical-outer face. Checking CORNER first -- as this
    function used to -- let a cell that is simultaneously corner-adjacent
    and genuinely cross-owned get bucketed as CORNER instead, which is only
    exit-critical under --strict; at interpolation_order=1 this hid an
    entire 1104-cell cross-patch leak behind a clean (exit 0) default-mode
    result. CROSS-PATCH must win that race, since it is always critical.
    """
    is_interior = is_own_interior(patch, i, j, k, angular_extent, radial_extent)

    if not is_interior:
        # Ghost/exterior point: is it inside the whole multipatch domain (a
        # genuine interpatch-interpolation target, per get_owner_patch) or
        # beyond it (a physical-BC ghost cell, which is this patch's own
        # concern)? Resolve real ownership mismatches first, ahead of the
        # geometric CORNER test below.
        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        inside_cube = max(abs(vx), abs(vy), abs(vz)) <= r0 * (1 + REL_TOL)
        inside_domain = inside_cube or r <= r1 * (1 + REL_TOL)
        if inside_domain:
            expected_patch = get_owner_patch(vx, vy, vz, r0)
            if expected_patch != patch:
                return "CROSS-PATCH", expected_patch

    if patch == CARTESIAN:
        # No genuine physical-outer face at all (bordered by wedges on all
        # six faces) -- corner cells therefore cannot occur here.
        near_outer = False
        near_other = (near_lo(i, band) or near_hi(i, angular_extent, band) or
                      near_lo(j, band) or near_hi(j, angular_extent, band) or
                      near_lo(k, band) or near_hi(k, angular_extent, band))
    else:
        # Wedge patch (cubed_sphere.cxx make_patch): the high side of the
        # radial (k) axis is the only genuine physical-outer face; both
        # angular axes (i, j) and the low/cartesian-facing side of k are all
        # interpatch.
        near_outer = near_hi(k, radial_extent, band)
        near_other = (near_lo(i, band) or near_hi(i, angular_extent, band) or
                      near_lo(j, band) or near_hi(j, angular_extent, band) or
                      near_lo(k, band))
    if near_outer and near_other:
        return "CORNER", None

    if is_interior:
        if patch == CARTESIAN:
            axes = ((i, angular_extent), (j, angular_extent), (k, angular_extent))
        else:
            axes = ((i, angular_extent), (j, angular_extent), (k, radial_extent))
        deep_interior = all(min(idx, extent - idx) > band for idx, extent in axes)
        return ("INTERIOR", None) if deep_interior else ("LOCAL", None)

    # Ghost/exterior, not cross-owned (checked above), not corner: either a
    # physical-BC ghost cell beyond the whole multipatch domain, or -- rare,
    # but not impossible -- an inside-domain point this same patch already
    # genuinely owns.
    if not inside_domain:
        return "LOCAL", None
    return "LOCAL", expected_patch


class Report:
    def __init__(self):
        self.total = 0
        self.checked = 0
        self.untouched = 0
        self.cleaned_at_source = 0
        self.missing_coords = 0
        self.missing_pre = 0
        self.by_kind = {}

    def flag(self, kind, row):
        self.by_kind.setdefault(kind, []).append(row)

    def summary(self):
        return "\n".join([
            f"Total rows read:                  {self.total}",
            f"Rows missing a coords match:       {self.missing_coords}",
            f"Rows missing a nan_test_pre match:  {self.missing_pre}",
            f"Rows checked against ground truth: {self.checked}",
            f"Untouched (1.0 -> 1.0):             {self.untouched}",
            f"Cleaned at source (NaN -> 1.0):     {self.cleaned_at_source}",
        ])


def process_file(nan_path, pre_path, coords_path, r0, r1, angular_cells,
                 radial_cells, band, report, patch_overlap=0):
    coords = load_global_coords(coords_path)
    pre_values = load_value_table(pre_path)

    angular_extent = angular_cells + 2 * patch_overlap
    radial_extent = radial_cells + patch_overlap

    for fields in read_tsv_rows(nan_path):
        if len(fields) < 12:
            continue
        _it, _t, patch_s, level_s, comp_s, i_s, j_s, k_s, _x, _y, _z, post_s = fields[:12]
        patch, level, comp = int(patch_s), int(level_s), int(comp_s)
        i, j, k = int(i_s), int(j_s), int(k_s)
        post = float(post_s)
        report.total += 1

        key = (patch, level, comp, i, j, k)
        if key not in coords:
            report.missing_coords += 1
            continue
        if key not in pre_values:
            report.missing_pre += 1
            continue
        vx, vy, vz = coords[key]
        pre = pre_values[key]
        report.checked += 1

        pre_is_nan = pre != pre
        post_is_nan = post != post
        pre_is_sentinel = (not pre_is_nan) and abs(pre - SENTINEL) <= SENTINEL_ATOL

        if pre_is_sentinel and not post_is_nan and abs(post - SENTINEL) <= SENTINEL_ATOL:
            report.untouched += 1
            continue

        row = (patch, i, j, k, vx, vy, vz, "pre", pre, "post", post)

        if pre_is_nan and post_is_nan:
            report.flag("informational: SURVIVED-AT-SOURCE (corruption persisted through sync)", row)
            continue

        if pre_is_nan and not post_is_nan:
            if abs(post - SENTINEL) > SENTINEL_ATOL:
                report.flag("informational: CLEANED-AT-SOURCE (cleaned, but to a non-1.0 value)", row)
            else:
                report.cleaned_at_source += 1
            continue

        if pre_is_sentinel and post_is_nan:
            subkind, expected_patch = classify_leak(
                patch, i, j, k, vx, vy, vz, r0, r1, angular_extent, radial_extent, band)
            extra = row if expected_patch is None else row + ("expected patch", expected_patch, PATCH_NAMES.get(expected_patch))
            report.flag(f"LEAK-{subkind}", extra)
            continue

        if pre_is_sentinel:
            # 1.0 -> finite, != 1.0: not sub-classified, mirrors
            # check_color.py's DEFAULT-LEAK -- always a plain bug report.
            report.flag("UNEXPECTED-VALUE-LEAK", row)
            continue

        # pre is neither NaN nor (approximately) the 1.0 sentinel: nan_test's
        # injection routine only ever writes 1.0 or NaN, so this indicates
        # nan_test_pre itself was corrupted before SYNC ever ran.
        report.flag("UNEXPECTED-PRE-VALUE (nan_test_pre itself is not 1.0 or NaN)", row)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nan_test_tsv", nargs="+",
                    help="post-SYNC capyrx_testmultipatch-nan_test TSV file(s) to check")
    ap.add_argument("--pre", default=None,
                    help="matching nan_test_pre TSV (auto-detected by default; "
                         "only valid with a single nan_test_tsv argument)")
    ap.add_argument("--coords", default=None,
                    help="matching coordinatesx-vertex_coords TSV (auto-detected by default; "
                         "only valid with a single nan_test_tsv argument)")
    ap.add_argument("--inner-boundary", type=float, default=None)
    ap.add_argument("--outer-boundary", type=float, default=None)
    ap.add_argument("--angular-cells", type=int, default=None)
    ap.add_argument("--radial-cells", type=int, default=None)
    ap.add_argument("--patch-overlap", type=int, default=None,
                    help="CapyrX_MultiPatch::patch_overlap (defaults to the value "
                         "found in the parameter file, or 0 if absent)")
    ap.add_argument("--band", type=int, default=None,
                    help="cell distance from a face counted as \"near\" it, for "
                         "LEAK-CORNER/LEAK-INTERIOR classification (defaults to "
                         "CarpetX::ghost_size found in the parameter file, or 2 "
                         "if neither is available)")
    ap.add_argument("--strict", action="store_true",
                    help="also exit non-zero if any LEAK-CORNER row is found")
    ap.add_argument("--max-examples", type=int, default=10, help="examples to print per issue kind")
    args = ap.parse_args()

    r0, r1 = args.inner_boundary, args.outer_boundary
    angular_cells, radial_cells = args.angular_cells, args.radial_cells
    patch_overlap = args.patch_overlap
    band = args.band

    par_path = find_param_file(args.nan_test_tsv[0])
    if None in (r0, r1, angular_cells, radial_cells, patch_overlap) and par_path:
        params = parse_par_params(par_path)
        r0 = r0 if r0 is not None else params.get("inner_boundary_radius")
        r1 = r1 if r1 is not None else params.get("outer_boundary_radius")
        angular_cells = angular_cells if angular_cells is not None else params.get("angular_cells")
        radial_cells = radial_cells if radial_cells is not None else params.get("radial_cells")
        patch_overlap = patch_overlap if patch_overlap is not None else params.get("patch_overlap")
    if band is None and par_path:
        band = parse_ghost_size(par_path)

    missing = [name for name, v in [
        ("--inner-boundary", r0), ("--outer-boundary", r1),
        ("--angular-cells", angular_cells), ("--radial-cells", radial_cells),
    ] if v is None]
    if missing:
        ap.error(f"could not determine {', '.join(missing)}; pass explicitly")

    angular_cells, radial_cells = int(angular_cells), int(radial_cells)
    patch_overlap = int(patch_overlap) if patch_overlap is not None else 0
    band = int(band) if band is not None else 2

    if (args.pre or args.coords) and len(args.nan_test_tsv) > 1:
        ap.error("--pre/--coords can only be used with a single nan_test_tsv argument")

    report = Report()
    for nan_path in args.nan_test_tsv:
        pre_path = Path(args.pre) if args.pre else guess_pre_path(nan_path)
        if pre_path is None or not Path(pre_path).exists():
            ap.error(f"could not find a matching nan_test_pre TSV for {nan_path}; pass --pre explicitly")
        coords_path = Path(args.coords) if args.coords else guess_coords_path(nan_path)
        if coords_path is None or not Path(coords_path).exists():
            ap.error(f"could not find a matching coordinatesx-vertex_coords TSV for {nan_path}; "
                     f"pass --coords explicitly")
        process_file(nan_path, pre_path, coords_path, r0, r1, angular_cells, radial_cells,
                     band, report, patch_overlap=patch_overlap)

    print(report.summary())
    print()

    any_critical = False
    always_critical = {"UNEXPECTED-VALUE-LEAK", "UNEXPECTED-PRE-VALUE (nan_test_pre itself is not 1.0 or NaN)",
                       "LEAK-INTERIOR", "LEAK-CROSS-PATCH"}
    for kind, items in sorted(report.by_kind.items(), key=lambda kv: -len(kv[1])):
        critical = kind in always_critical or (args.strict and kind == "LEAK-CORNER")
        any_critical = any_critical or (critical and items)
        print(f"=== {kind} ({len(items)}) ===")
        for row in items[: args.max_examples]:
            print("  ", row)
        if len(items) > args.max_examples:
            print(f"   ... and {len(items) - args.max_examples} more")
        print()

    if any_critical:
        print("RESULT: NaN corruption leaked beyond its injection site.")
        sys.exit(1)
    else:
        print("RESULT: no leak beyond the injection site found (LEAK-CORNER present but "
              "not --strict, if any is listed above).")
        sys.exit(0)


if __name__ == "__main__":
    main()
