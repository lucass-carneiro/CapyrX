#!/usr/bin/env python3
"""
Interpret output from the CapyrX_TestMultiPatch smooth-field outer-boundary
A/B test -- the "Channel 1" decisive experiment of mp_noise_3.md sec 8.

Background (mp_noise_3.md sec 3/4/8)
-----------------------------------
`CapyrX_TestMultiPatch_write_smooth_test` fills `smooth` with a smooth analytic
field (`smooth_field`, default z_global) *everywhere* -- interior AND every
ghost cell, so each patch's outer-BC ghost cells (radial k=35,36 in the
r=32/overlap=2 cubed-sphere layout) start out holding the *exact* analytic
value. `smooth_pre` is snapshotted immediately afterwards and never synced, so
it is an in-place, per-cell EXACT ground truth. "SYNC: smooth" then runs
CarpetX's full boundary-condition + MultiPatch_Interpolate pipeline using
whichever CarpetX::boundary_* kind the .par file selected:

    smooth_z_neumann.par    boundary_* = neumann              -> (a), Neumann leg
    smooth_z_linextrap.par  boundary_* = linear extrapolation -> (a), lin-extrap leg
    smooth_z_none.par       boundary_* = none                 -> (c), exact ghosts

Because the interior (including the r=2.5 boundary vertex k=34) is written
identically in all three runs and is *not* touched by any outer BC, the three
runs differ ONLY in their outer-BC ghost cells (k=35,36). Any post-SYNC
difference at a cell that is NOT itself an outer-BC ghost therefore must have
arrived through the interpatch interpolation stencil reading a donor patch's
outer-BC ghosts -- i.e. through Channel 1.

What this script reports
------------------------
For every cell present in all supplied runs it computes, against the exact
`smooth_pre`:

    d_ab      = smooth(neumann)   - smooth(linextrap)     Channel-1 BC-distinctness
    err_neu   = smooth(neumann)   - exact
    err_lin   = smooth(linextrap) - exact
    err_none  = smooth(none)      - exact                 (c): interp error, exact ghosts

and classifies every cell with |d_ab| > --tol (i.e. every BC-sensitive cell):

    OUTER-GHOST-SELF        cell lies OUTSIDE the multipatch domain (r > r1, or
                            outside the cube) -- a patch's own physical outer-BC
                            ghost. Expected to differ (neumann vs lin-extrap fill
                            it differently); NOT the invariant violation.
    INTERPATCH-CROSS-PATCH  cell is INSIDE the domain but its geometric owner
                            (get_owner_patch) is a DIFFERENT patch than the one
                            reporting it: a neighbouring patch's interpatch-ghost
                            value that changed purely because of the outer-BC
                            choice. THIS is the mp_noise_3 sec-3 invariant
                            violation ("Channel 1"). Delta != 0 here is the sec-8
                            primary positive result.
    INTERPATCH-SAME-PATCH   inside the domain, owner == reporting patch: a
                            same-patch overlap/ghost receiver whose fill read an
                            outer-BC ghost. Reported separately (relevant but not
                            a cross-patch seam crossing).

The sec-8 decision, read off the summary:
  * max |d_ab| over INTERPATCH-CROSS-PATCH  >  --tol  ==> Channel 1 is a LIVE
    BC-distinct injector (literal invariant violated with *used* values).
    == 0 ==> Channel 1 is inert here (document-only).
  * At the cross-patch receivers, compare |err_none| (exact ghosts, (c)) against
    |err_lin| / |err_neu| ((a)): if |err_none| < |err_lin|, outer-ghost fidelity
    materially improves the interpatch value -- evidence for the mp_noise_3 sec-7.3
    "raise BC fidelity" remedy. (Measurement (b), the Step-B one-sided stencil, is
    deliberately NOT produced here; it needs the shared-CarpetX clamp deferred per
    mp_noise_3 sec 11.)

Usage
-----
    check_smooth_test.py \\
        --neumann   exe/smooth_z_neumann/capyrx_testmultipatch-smooth.it000000.p0000.tsv \\
        --linextrap exe/smooth_z_linextrap/capyrx_testmultipatch-smooth.it000000.p0000.tsv \\
        --none      exe/smooth_z_none/capyrx_testmultipatch-smooth.it000000.p0000.tsv

The companion smooth_pre and coordinatesx-vertex_coords TSVs are auto-located
next to the --none file (the exact-ghost run) unless --pre/--coords are given.
Geometry (inner/outer radius, angular/radial cells, patch_overlap) is recovered
from the parameter file named in the --none TSV header, or may be overridden.

At least --neumann and --linextrap are required (the primary Delta signal).
--none is optional but needed for the (a)-vs-(c) fidelity comparison.
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
    in_overlap_band,
    load_global_coords,
    parse_par_params,
    read_tsv_rows,
)


def guess_companion(smooth_path, target_stem):
    """Locate a sibling TSV of `smooth_path` whose leading variable name is
    replaced by `target_stem` (e.g. smooth_test -> smooth_test_pre, or the
    coords stem).

    CarpetX names TSVs by *group*, and the smooth field's group is `smooth_test`
    (real files are `capyrx_testmultipatch-smooth_test...` and `...-smooth_test_pre...`),
    so the match must anchor on `-smooth_test`, not `-smooth`: anchoring on the
    shorter `-smooth` left `_test` dangling and synthesized non-existent names
    (`...-smooth_pre_test...`), which is why --pre/--coords had to be passed by
    hand. See 12_9_impl.md Phase 3(a)."""
    p = Path(smooth_path)
    name = re.sub(r"^capyrx_testmultipatch-smooth_test(?:_pre)?", target_stem, p.name)
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def load_value_table(path):
    """key (patch, level, component, i, j, k) -> smooth(_pre) value."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) < 12:
            continue
        _it, _t, patch, level, comp, i, j, k, _x, _y, _z, value = fields[:12]
        key = (int(patch), int(level), int(comp), int(i), int(j), int(k))
        table[key] = float(value)
    return table


def inside_domain(vx, vy, vz, r0, r1):
    r = math.sqrt(vx * vx + vy * vy + vz * vz)
    inside_cube = max(abs(vx), abs(vy), abs(vz)) <= r0 * (1 + REL_TOL)
    return inside_cube or r <= r1 * (1 + REL_TOL)


def classify(patch, vx, vy, vz, r0, r1):
    """Return (category, owner) for a BC-sensitive cell."""
    if not inside_domain(vx, vy, vz, r0, r1):
        return "OUTER-GHOST-SELF", None
    owner = get_owner_patch(vx, vy, vz, r0)
    if owner != patch:
        return "INTERPATCH-CROSS-PATCH", owner
    return "INTERPATCH-SAME-PATCH", owner


def is_strict_interior(patch, i, j, k, vx, vy, vz, r0, r1,
                       angular_cells, radial_cells, patch_overlap):
    """True iff this cell is unambiguously in the reporting patch's own *true*
    (non-ghost, non-overlap) interior: in-domain, geometrically owned by this
    same patch, its index inside the patch's widened interior range, and NOT in
    the patch_overlap band. Such a cell is valid data written analytically and
    is touched by neither the outer BC nor the interpatch fill, so it MUST be
    bit-identical across the neumann/linextrap/none runs -- that identity is the
    A/B experiment's premise. Index/overlap logic mirrors check_color.py's
    own_interior + in_overlap_band (shared source: CapyrX_MultiPatch make_patch);
    the overlap band is deliberately excluded to keep the predicate conservative
    (only cells that cannot be a donor-source slice or an interface tie-break)."""
    if not inside_domain(vx, vy, vz, r0, r1):
        return False
    if get_owner_patch(vx, vy, vz, r0) != patch:
        return False
    angular_extent = angular_cells + 2 * patch_overlap
    radial_extent = radial_cells + patch_overlap
    if patch == CARTESIAN:
        own_interior = (0 <= i <= angular_extent and 0 <= j <= angular_extent
                        and 0 <= k <= angular_extent)
    else:
        own_interior = (0 <= i <= angular_extent and 0 <= j <= angular_extent
                        and 0 <= k <= radial_extent)
    if not own_interior:
        return False
    return not in_overlap_band(patch, i, j, k, angular_cells, radial_cells, patch_overlap)


class Bucket:
    def __init__(self):
        self.rows = []       # (patch, i, j, k, vx, vy, vz, d_ab, err_neu, err_lin, err_none, owner)
        self.max_dab = 0.0
        self.max_dab_row = None

    def add(self, row):
        self.rows.append(row)
        d = abs(row[7])
        if d > self.max_dab:
            self.max_dab = d
            self.max_dab_row = row

    def __len__(self):
        return len(self.rows)


def fmt_row(row):
    patch, i, j, k, vx, vy, vz, d_ab, e_neu, e_lin, e_none, owner = row
    r = math.sqrt(vx * vx + vy * vy + vz * vz)
    owner_s = "" if owner is None else f" owner={owner}({PATCH_NAMES.get(owner)})"
    en = "n/a" if e_none is None else f"{e_none:+.3e}"
    return (f"  patch={patch} ijk=({i},{j},{k}) r={r:.4f}"
            f" d_ab={d_ab:+.3e} err_neu={e_neu:+.3e} err_lin={e_lin:+.3e}"
            f" err_none={en}{owner_s}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--neumann", required=True,
                    help="post-SYNC smooth TSV from smooth_z_neumann.par")
    ap.add_argument("--linextrap", required=True,
                    help="post-SYNC smooth TSV from smooth_z_linextrap.par")
    ap.add_argument("--none", default=None,
                    help="post-SYNC smooth TSV from smooth_z_none.par (exact ghosts, "
                         "measurement (c); optional but enables the (a)-vs-(c) comparison)")
    ap.add_argument("--pre", default=None,
                    help="exact ground-truth smooth_pre TSV (auto-detected next to --none, "
                         "else --neumann)")
    ap.add_argument("--coords", default=None,
                    help="coordinatesx-vertex_coords TSV (auto-detected next to --none, "
                         "else --neumann)")
    ap.add_argument("--inner-boundary", type=float, default=None)
    ap.add_argument("--outer-boundary", type=float, default=None)
    ap.add_argument("--angular-cells", type=int, default=None)
    ap.add_argument("--radial-cells", type=int, default=None)
    ap.add_argument("--patch-overlap", type=int, default=None)
    ap.add_argument("--tol", type=float, default=1e-14,
                    help="|neumann - linextrap| above which a cell is counted as "
                         "BC-sensitive (default 1e-14, i.e. above bit-noise)")
    ap.add_argument("--max-examples", type=int, default=10,
                    help="example rows to print per category (default 10)")
    args = ap.parse_args()

    # Geometry: prefer the --none par file (the exact/reference run), else --neumann.
    r0, r1 = args.inner_boundary, args.outer_boundary
    angular_cells, radial_cells = args.angular_cells, args.radial_cells
    patch_overlap = args.patch_overlap
    geom_src = args.none or args.neumann
    par_path = find_param_file(geom_src)
    if None in (r0, r1, angular_cells, radial_cells) and par_path:
        params = parse_par_params(par_path)
        r0 = r0 if r0 is not None else params.get("inner_boundary_radius")
        r1 = r1 if r1 is not None else params.get("outer_boundary_radius")
        angular_cells = angular_cells if angular_cells is not None else params.get("angular_cells")
        radial_cells = radial_cells if radial_cells is not None else params.get("radial_cells")
        patch_overlap = patch_overlap if patch_overlap is not None else params.get("patch_overlap")
    missing = [n for n, v in [("--inner-boundary", r0), ("--outer-boundary", r1)] if v is None]
    if missing:
        ap.error(f"could not determine {', '.join(missing)}; pass explicitly")
    r0, r1 = float(r0), float(r1)
    angular_cells = int(angular_cells) if angular_cells is not None else None
    radial_cells = int(radial_cells) if radial_cells is not None else None
    patch_overlap = int(patch_overlap) if patch_overlap is not None else 0
    # The interior-invariant guard (fix b) needs the patch cell counts to know
    # each patch's true interior index range. If they couldn't be recovered, the
    # guard is skipped -- but loudly, never silently (harsh-critic requirement).
    interior_check_enabled = angular_cells is not None and radial_cells is not None

    # Companions: exact (pre) + coords, from --none if present else --neumann.
    comp_src = args.none or args.neumann
    pre_path = Path(args.pre) if args.pre else guess_companion(comp_src, "capyrx_testmultipatch-smooth_test_pre")
    coords_path = Path(args.coords) if args.coords else guess_companion(comp_src, "coordinatesx-vertex_coords")
    if pre_path is None or not Path(pre_path).exists():
        ap.error("could not find the exact smooth_pre TSV; pass --pre explicitly")
    if coords_path is None or not Path(coords_path).exists():
        ap.error("could not find the coordinatesx-vertex_coords TSV; pass --coords explicitly")

    coords = load_global_coords(coords_path)
    exact = load_value_table(pre_path)
    neu = load_value_table(args.neumann)
    lin = load_value_table(args.linextrap)
    non = load_value_table(args.none) if args.none else None

    # Sanity: the interior must be bit-identical across the runs. Any interior
    # (in-domain, owner==patch, not a ghost/overlap region) mismatch would break
    # the whole premise, so surface it loudly (fix b).
    interior_mismatch = 0
    interior_examples = []
    nonfinite = 0
    total = 0
    buckets = {"OUTER-GHOST-SELF": Bucket(),
               "INTERPATCH-CROSS-PATCH": Bucket(),
               "INTERPATCH-SAME-PATCH": Bucket()}

    for key, en in neu.items():
        if key not in lin or key not in coords or key not in exact:
            continue
        total += 1
        el = lin[key]
        d_ab = en - el

        patch, _level, _comp, i, j, k = key
        vx, vy, vz = coords[key]
        ex = exact[key]
        err_neu = en - ex
        err_lin = el - ex
        err_none = (non[key] - ex) if (non is not None and key in non) else None

        # Non-finite skip (fix c): with smooth_field=one_over_r the origin cell
        # is singular and d_ab/err_* are NaN/inf. `abs(NaN) <= tol` is False, so
        # without this skip the cell would fall through and be mis-bucketed with
        # a NaN delta. Count them separately; a nonzero count with a non-singular
        # field is itself suspicious (reported in the summary).
        if not math.isfinite(d_ab):
            nonfinite += 1
            continue

        # Interior-invariant guard (fix b): a strict-interior cell must be
        # identical across neumann/linextrap (and none, when supplied). Check it
        # before the |d_ab|<=tol skip below, since a passing interior cell has
        # d_ab==0 and would otherwise never be examined.
        if interior_check_enabled and is_strict_interior(
                patch, i, j, k, vx, vy, vz, r0, r1,
                angular_cells, radial_cells, patch_overlap):
            bad = abs(d_ab) > args.tol
            if non is not None and key in non and math.isfinite(non[key]):
                bad = bad or abs(en - non[key]) > args.tol
            if bad:
                interior_mismatch += 1
                interior_examples.append((patch, i, j, k, vx, vy, vz, d_ab,
                                          err_neu, err_lin, err_none, patch))

        if abs(d_ab) <= args.tol:
            continue

        category, owner = classify(patch, vx, vy, vz, r0, r1)
        buckets[category].add((patch, i, j, k, vx, vy, vz, d_ab,
                               err_neu, err_lin, err_none, owner))

    # A broken interior premise invalidates the entire Channel-1 result, so abort
    # here -- before any bucket/verdict output -- rather than reporting a
    # meaningless signal (fix b). exit(2) distinguishes "premise broken / checker
    # can't trust its input" from the exit(1) "Channel 1 is LIVE" positive.
    if interior_mismatch > 0:
        print("!" * 78, file=sys.stderr)
        print(f"INTERIOR-INVARIANT VIOLATED: {interior_mismatch} strict-interior "
              "cell(s) differ across the runs.", file=sys.stderr)
        print("The A/B premise (each patch's true interior is bit-identical across "
              "neumann/linextrap/none, being untouched by any outer BC or interpatch "
              "fill) is broken; the Channel-1 result would be meaningless. Aborting.",
              file=sys.stderr)
        for row in interior_examples[: args.max_examples]:
            print(fmt_row(row), file=sys.stderr)
        if interior_mismatch > args.max_examples:
            print(f"   ... and {interior_mismatch - args.max_examples} more",
                  file=sys.stderr)
        print("!" * 78, file=sys.stderr)
        sys.exit(2)

    print("=" * 78)
    print("Smooth-field outer-BC (Channel 1) A/B test -- mp_noise_3.md sec 8")
    print("=" * 78)
    print(f"field/geometry: r0={r0} r1={r1} angular_cells={angular_cells} "
          f"radial_cells={radial_cells} patch_overlap={patch_overlap}")
    print(f"cells compared across neumann/linextrap: {total}")
    print(f"non-finite cells skipped (NaN/inf, e.g. one_over_r origin): {nonfinite}")
    print(f"BC-sensitivity tolerance (|d_ab| >): {args.tol:g}")
    print(f"(c) exact-ghost run supplied: {'yes' if non is not None else 'no'}")
    if interior_check_enabled:
        print("interior-invariant guard: PASSED (0 strict-interior mismatches)")
    else:
        print("interior-invariant guard: SKIPPED (angular/radial cell counts "
              "unavailable; pass --angular-cells/--radial-cells to enable)")
    print()

    for name in ("INTERPATCH-CROSS-PATCH", "INTERPATCH-SAME-PATCH", "OUTER-GHOST-SELF"):
        b = buckets[name]
        print(f"=== {name} ({len(b)} BC-sensitive cells, max|d_ab|={b.max_dab:.3e}) ===")
        for row in b.rows[: args.max_examples]:
            print(fmt_row(row))
        if len(b) > args.max_examples:
            print(f"   ... and {len(b) - args.max_examples} more")
        if b.max_dab_row is not None:
            print("  peak |d_ab| cell:")
            print(fmt_row(b.max_dab_row))
        # (a)-vs-(c) accuracy comparison for the interpatch receivers
        if name.startswith("INTERPATCH") and non is not None and len(b):
            n_none_better = sum(1 for r in b.rows
                                if r[10] is not None and abs(r[10]) < abs(r[9]))
            max_gap = max((abs(r[9]) - abs(r[10]) for r in b.rows if r[10] is not None),
                          default=0.0)
            print(f"  (a)-vs-(c): |err_none| < |err_lin| at {n_none_better}/{len(b)} "
                  f"cells; max (|err_lin| - |err_none|) = {max_gap:+.3e}")
        print()

    cross = buckets["INTERPATCH-CROSS-PATCH"]
    print("-" * 78)
    if cross.max_dab > args.tol:
        print(f"RESULT: Channel 1 is LIVE. The outer-BC choice changes {len(cross)} "
              f"cross-patch interpatch value(s); peak |neumann - linextrap| = "
              f"{cross.max_dab:.3e}.")
        print("        => mp_noise_3 sec-3 invariant violated with *used* values "
              "(sec-8 primary positive).")
        if non is not None:
            print("        See the (a)-vs-(c) line above to judge whether outer-ghost "
                  "fidelity (sec 7.3) materially helps.")
        sys.exit(1)
    else:
        print("RESULT: Channel 1 is INERT here. No cross-patch interpatch value changed "
              f"above tol ({args.tol:g}).")
        print("        => sec-8 negative: document Channel 1 as a bounded, benign "
              "boundary-adjacency effect; do not touch the anchor policy.")
        sys.exit(0)


if __name__ == "__main__":
    main()
