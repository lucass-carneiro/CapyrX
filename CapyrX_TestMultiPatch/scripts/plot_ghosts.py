#!/usr/bin/env python3
"""
3D scatter plot of one patch's ghost/boundary points, colored by the
CapyrX_TestMultiPatch coloring test's `color` field, to eyeball whether
post-SYNC interpatch interpolation filled them with geometrically sensible
data.

This is a visual companion to check_color.py (same directory), whose TSV
format, column layout, and cubed-sphere ground truth (get_owner_patch) it
imports and reuses rather than re-deriving. See that script's docstring for
the full explanation of the color/color_pre encoding.

Point selection ("what counts as ghost/boundary"): a row of the requested
patch is plotted iff its color_pre category is NOT (INTERIOR or OVERLAP) of
that same patch. color_pre is write_color's pre-SYNC snapshot, so this is a
direct read of what CarpetX itself classified the point as (loop_bnd_device
or loop_ghosts_device fired there), not a re-implementation of the
angular_cells/radial_cells/patch_overlap index arithmetic. Rows where
color_pre decodes to something else entirely (a poison value SYNC never
touched) are included too, since that is itself a bug worth seeing.

Coloring: each point's hue is the *source* patch decoded from the (post-SYNC)
`color` value -- i.e. which patch's data this ghost/boundary cell was
actually filled from. A correct interpatch sync should paint a contiguous
ghost/boundary shell in the hue of whichever single neighboring patch
geometrically owns that region; a wrong-neighbor bug shows up as a
patch-colored speckle sitting inside a differently-colored region. Marker
shape encodes the pre-SYNC structural category (ghost vs. boundary), and an
optional red-ring overlay (on by default, --no-mismatch-overlay to disable)
independently recomputes the expected owning patch from the point's global
coordinates via check_color.get_owner_patch and rings any point that
disagrees -- the same ground truth check_color.py runs, just drawn instead
of tabulated.

Patch outlines: for spatial context, every patch's own true (pre-SYNC
INTERIOR/OVERLAP, i.e. not-ghost) grid is used to draw a wireframe of that
patch's 6 bounding faces, built from the *actual* vertex_coords data rather
than the analytic cubed-sphere mapping -- this renders the wedge patches'
real curvature and, per this repo's existing philosophy (see check_color.py),
doesn't trust bookkeeping that could itself be the thing under test. The
requested patch's own outline is drawn bold/black; the rest are thin and
faint, just for orientation.

Usage:
    plot_ghosts.py 0 exe/color_ghost_overlap/capyrx_testmultipatch-color.it000000.p0000.tsv

    plot_ghosts.py 3 COLOR.tsv --save patch3.png

The matching color_pre and coordinatesx-vertex_coords TSVs are located
automatically next to the given color TSV (like check_color.py's --coords
auto-detection) unless --color-pre/--coords are given explicitly. The
mismatch overlay's inner/outer boundary radius are likewise recovered from
the parameter file referenced in the TSV header, unless passed explicitly.
"""

import argparse
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_color as cc

STRUCT_MARKERS = {cc.GHOST: "^", cc.BOUNDARY: "o"}
DEFAULT_MARKER = "x"


def guess_color_pre_path(color_path):
    p = Path(color_path)
    if "-color_pre." in p.name:
        return p if p.exists() else None
    name = p.name.replace("-color.", "-color_pre.")
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def patch_color(key, palette):
    if key == "zero":
        return (0.6, 0.6, 0.6)
    if key == "poison":
        return (0.0, 0.0, 0.0)
    return palette[key % len(palette)]


def scan_color_pre(color_pre_path, coords, query_patch):
    """Single pass over color_pre building:
    - per_patch_interior: patch -> [(i, j, k, vx, vy, vz)] for that patch's
      own true (non-ghost) grid, used to draw patch outlines.
    - pre_lookup: (level, comp, i, j, k) -> (category, stored_patch) for
      query_patch only, used to decide which of its rows are ghost/boundary.
    """
    per_patch_interior = defaultdict(list)
    pre_lookup = {}
    for fields in cc.read_tsv_rows(color_pre_path):
        if len(fields) < 12:
            continue
        _it, _t, patch_s, level_s, comp_s, i_s, j_s, k_s, _x, _y, _z, color_s = fields[:12]
        patch = int(patch_s)
        level, comp, i, j, k = int(level_s), int(comp_s), int(i_s), int(j_s), int(k_s)
        category, stored_patch = cc.decode_color(float(color_s))

        if category in (cc.INTERIOR, cc.OVERLAP) and stored_patch == patch:
            key = (patch, level, comp, i, j, k)
            if key in coords:
                vx, vy, vz = coords[key]
                per_patch_interior[patch].append((i, j, k, vx, vy, vz))

        if patch == query_patch:
            pre_lookup[(level, comp, i, j, k)] = (category, stored_patch)

    return per_patch_interior, pre_lookup


def load_ghost_boundary_points(query_patch, color_path, pre_lookup, coords):
    points = []
    for fields in cc.read_tsv_rows(color_path):
        if len(fields) < 12:
            continue
        _it, _t, patch_s, level_s, comp_s, i_s, j_s, k_s, _x, _y, _z, color_s = fields[:12]
        patch = int(patch_s)
        if patch != query_patch:
            continue
        level, comp, i, j, k = int(level_s), int(comp_s), int(i_s), int(j_s), int(k_s)

        pre_cat, pre_patch = pre_lookup.get((level, comp, i, j, k), (None, None))
        if pre_cat in (cc.INTERIOR, cc.OVERLAP) and pre_patch == patch:
            continue  # this patch's own true interior/overlap pre-SYNC -- not ghost/boundary

        key = (patch, level, comp, i, j, k)
        if key not in coords:
            continue
        vx, vy, vz = coords[key]
        post_color = float(color_s)
        post_cat, post_patch = cc.decode_color(post_color)

        points.append(dict(i=i, j=j, k=k, vx=vx, vy=vy, vz=vz, pre_cat=pre_cat,
                            post_cat=post_cat, post_patch=post_patch, post_color=post_color))
    return points


def classify_mismatch(p, r0, r1, owned_index_set, query_patch):
    """Re-derive expected ground truth for one ghost/boundary point, mirroring
    check_color.py's process_file branches (its pre_sync informational cases
    don't apply here since we always look at post-SYNC `color`).

    A ghost/boundary row can arise two different ways, and they need two
    different ground truths (this is exactly the own_interior split
    check_color.py's process_file makes via angular_cells/radial_cells/
    patch_overlap index arithmetic -- reproduced here from the data instead,
    via owned_index_set, to avoid needing those params):

    - INTRA-patch: (i, j, k) is, at some *other* component (AMReX box) of
      this same patch, that component's own true INTERIOR/OVERLAP (i.e. the
      point is genuinely inside this patch's own domain, just owned by a
      neighboring box) -- owned_index_set captures exactly this. The only
      correct fill is this same patch's own marker; get_owner_patch's global
      coordinates are irrelevant here (they'd wrongly re-litigate an already
      patch-internal point against the whole cubed-sphere geometry -- this
      was a real false positive during development, see coloring_viz_0.md).
    - INTER-patch: no component owns (i, j, k) as true interior -- this
      location is genuinely outside this patch's own domain, so the
      cross-check falls back to get_owner_patch's independently-recomputed
      global ground truth, same as check_color.py's ghost-point branch.

    Corner/edge points where get_owner_patch's tie-break is genuinely
    ambiguous are NOT suppressed in the inter-patch branch (unlike
    check_color.py's corner_tag, which needs the index-arithmetic params
    this script avoids) -- look for a *small cluster* of rings sitting
    exactly on a multi-patch seam before treating it as a bug.
    """
    vx, vy, vz, post_cat, post_patch, post_color = (
        p["vx"], p["vy"], p["vz"], p["post_cat"], p["post_patch"], p["post_color"])

    if post_cat is None:
        return True, "NON-INTEGER-COLOR (poison value never overwritten)"

    if (p["i"], p["j"], p["k"]) in owned_index_set:
        if post_cat in (cc.INTERIOR, cc.OVERLAP) and post_patch == query_patch:
            return False, None
        if post_cat == cc.GHOST and post_patch == query_patch:
            return True, "INTRA-PATCH-GHOST-NOT-SYNCED (box exchange never refilled this point)"
        return True, (f"OWN-INTERIOR-MISMATCH (another box of this same patch owns this point "
                       f"as true interior, but this row shows {cc.PATCH_NAMES.get(post_patch, post_patch)})")

    r = math.sqrt(vx * vx + vy * vy + vz * vz)
    inside_cube = max(abs(vx), abs(vy), abs(vz)) <= r0 * (1 + cc.REL_TOL)
    inside_domain = inside_cube or r <= r1 * (1 + cc.REL_TOL)

    if not inside_domain:
        if post_color != 0.0:
            return True, "EXTERIOR-NOT-ZEROED (physical BC did not overwrite this point)"
        return False, None

    if post_cat == cc.GHOST:
        return True, "GHOST-MARKER-LEAK (sourced from a not-yet-filled ghost cell)"
    if post_cat == 0:
        return True, "DEFAULT-LEAK (Dirichlet-zero leaked inside the valid domain)"

    expected_patch = cc.get_owner_patch(vx, vy, vz, r0)
    if post_patch != expected_patch:
        return True, (f"WRONG-NEIGHBOR (got {cc.PATCH_NAMES.get(post_patch, post_patch)}, "
                       f"expected {cc.PATCH_NAMES.get(expected_patch, expected_patch)})")
    return False, None


def resolve_radii(args, color_path):
    r0, r1 = args.inner_boundary, args.outer_boundary
    if r0 is None or r1 is None:
        par_path = cc.find_param_file(color_path)
        if par_path:
            params = cc.parse_par_params(par_path)
            r0 = r0 if r0 is not None else params.get("inner_boundary_radius")
            r1 = r1 if r1 is not None else params.get("outer_boundary_radius")
    return r0, r1


def sparse_values(sorted_unique, target):
    n = len(sorted_unique)
    if n <= target:
        return list(sorted_unique)
    idx = sorted(set(round(i * (n - 1) / (target - 1)) for i in range(target)))
    return [sorted_unique[i] for i in idx]


def draw_patch_outline(ax, points, color, lw, alpha, grid_lines):
    if not points:
        return
    index = {(i, j, k): (vx, vy, vz) for i, j, k, vx, vy, vz in points}
    all_i = sorted(set(i for i, j, k, *_ in points))
    all_j = sorted(set(j for i, j, k, *_ in points))
    all_k = sorted(set(k for i, j, k, *_ in points))
    if not (all_i and all_j and all_k):
        return

    faces = [
        ("i", all_i[0]), ("i", all_i[-1]),
        ("j", all_j[0]), ("j", all_j[-1]),
        ("k", all_k[0]), ("k", all_k[-1]),
    ]
    axis_values = {"i": all_i, "j": all_j, "k": all_k}

    for fixed_axis, fixed_val in faces:
        other_axes = [a for a in ("i", "j", "k") if a != fixed_axis]
        a1, a2 = other_axes
        vals1 = sparse_values(axis_values[a1], grid_lines)
        vals2 = axis_values[a2]
        for v1 in vals1:
            line = []
            for v2 in vals2:
                ijk = {fixed_axis: fixed_val, a1: v1, a2: v2}
                key = (ijk["i"], ijk["j"], ijk["k"])
                if key in index:
                    line.append(index[key])
            if len(line) > 1:
                xs, ys, zs = zip(*line)
                ax.plot(xs, ys, zs, color=color, lw=lw, alpha=alpha)

        vals1b = axis_values[a1]
        vals2b = sparse_values(axis_values[a2], grid_lines)
        for v2 in vals2b:
            line = []
            for v1 in vals1b:
                ijk = {fixed_axis: fixed_val, a1: v1, a2: v2}
                key = (ijk["i"], ijk["j"], ijk["k"])
                if key in index:
                    line.append(index[key])
            if len(line) > 1:
                xs, ys, zs = zip(*line)
                ax.plot(xs, ys, zs, color=color, lw=lw, alpha=alpha)


def plot_scatter(ax, points, palette, size):
    groups = defaultdict(list)
    for p in points:
        if p["post_cat"] is None:
            color_key = "poison"
        elif p["post_cat"] == 0:
            color_key = "zero"
        else:
            color_key = p["post_patch"]
        groups[(color_key, p["pre_cat"])].append(p)

    for (color_key, pre_cat), pts in groups.items():
        xs = [p["vx"] for p in pts]
        ys = [p["vy"] for p in pts]
        zs = [p["vz"] for p in pts]
        color = patch_color(color_key, palette)
        marker = STRUCT_MARKERS.get(pre_cat, DEFAULT_MARKER)
        struct_name = {cc.GHOST: "ghost", cc.BOUNDARY: "boundary"}.get(pre_cat, "other/poison")
        if color_key in ("zero", "poison"):
            source_name = color_key
        else:
            source_name = cc.PATCH_NAMES.get(color_key, str(color_key))
        ax.scatter(xs, ys, zs, color=[color] * len(pts), marker=marker, s=size,
                   depthshade=True, label=f"{source_name} / {struct_name} ({len(pts)})")

    flagged = [p for p in points if p.get("flagged")]
    if flagged:
        xs = [p["vx"] for p in flagged]
        ys = [p["vy"] for p in flagged]
        zs = [p["vz"] for p in flagged]
        ax.scatter(xs, ys, zs, s=size * 2.5, facecolors="none", edgecolors="red",
                   linewidths=1.3, marker="o", label=f"flagged mismatch ({len(flagged)})")


def set_equal_aspect(ax, points, per_patch_interior):
    xs = [p["vx"] for p in points]
    ys = [p["vy"] for p in points]
    zs = [p["vz"] for p in points]
    for pts in per_patch_interior.values():
        xs += [q[3] for q in pts]
        ys += [q[4] for q in pts]
        zs += [q[5] for q in pts]
    if not xs:
        return
    xr, yr, zr = max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)
    ax.set_box_aspect((xr or 1, yr or 1, zr or 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("patch", type=int, help="patch number whose ghost/boundary points to plot")
    ap.add_argument("color_tsv", help="post-SYNC capyrx_testmultipatch-color TSV")
    ap.add_argument("--color-pre", default=None,
                    help="matching color_pre TSV (auto-detected by default)")
    ap.add_argument("--coords", default=None,
                    help="matching coordinatesx-vertex_coords TSV (auto-detected by default)")
    ap.add_argument("--inner-boundary", type=float, default=None)
    ap.add_argument("--outer-boundary", type=float, default=None)
    ap.add_argument("--no-mismatch-overlay", action="store_true",
                    help="skip the get_owner_patch ground-truth ring overlay")
    ap.add_argument("--no-outlines", action="store_true", help="skip drawing patch outlines")
    ap.add_argument("--outline-grid-lines", type=int, default=6,
                    help="grid lines per face, per direction, in the patch outlines (default 6)")
    ap.add_argument("--point-size", type=float, default=14, help="scatter marker size")
    ap.add_argument("--save", default=None, help="also save a PNG to this path")
    ap.add_argument("--no-show", action="store_true",
                    help="don't open an interactive window (implies --save is useful)")
    args = ap.parse_args()

    color_path = Path(args.color_tsv)
    if not color_path.exists():
        ap.error(f"{color_path} does not exist")

    color_pre_path = Path(args.color_pre) if args.color_pre else guess_color_pre_path(color_path)
    if color_pre_path is None or not color_pre_path.exists():
        ap.error(f"could not find a matching color_pre TSV for {color_path}; pass --color-pre "
                 f"explicitly, or add CapyrX_TestMultiPatch::color_pre to CarpetX::out_tsv_vars "
                 f"and rerun")

    coords_path = Path(args.coords) if args.coords else cc.guess_coords_path(color_path)
    if coords_path is None or not coords_path.exists():
        ap.error(f"could not find a matching coordinatesx-vertex_coords TSV for {color_path}; "
                 f"pass --coords explicitly")

    print(f"loading vertex coordinates from {coords_path}...")
    coords = cc.load_global_coords(coords_path)

    print(f"loading pre-SYNC categories from {color_pre_path}...")
    per_patch_interior, pre_lookup = scan_color_pre(color_pre_path, coords, args.patch)

    if args.patch not in per_patch_interior:
        available = sorted(per_patch_interior)
        ap.error(f"no interior rows found for patch {args.patch} in {color_pre_path}; "
                 f"available patches: {available}")

    print(f"loading post-SYNC colors from {color_path}...")
    points = load_ghost_boundary_points(args.patch, color_path, pre_lookup, coords)

    if not points:
        print(f"patch {args.patch} has no ghost/boundary rows (fully interior in this file?); nothing to plot")
        return

    r0 = r1 = None
    if not args.no_mismatch_overlay:
        r0, r1 = resolve_radii(args, color_path)
        if r0 is None or r1 is None:
            print("warning: could not resolve --inner-boundary/--outer-boundary from the "
                  "parameter file; disabling mismatch overlay", file=sys.stderr)
            r0 = r1 = None

    if r0 is not None:
        owned_index_set = {(i, j, k) for i, j, k, *_ in per_patch_interior[args.patch]}
        for p in points:
            p["flagged"], p["reason"] = classify_mismatch(p, r0, r1, owned_index_set, args.patch)
    else:
        for p in points:
            p["flagged"], p["reason"] = False, None

    patch_name = cc.PATCH_NAMES.get(args.patch, str(args.patch))
    print(f"\npatch {args.patch} ({patch_name}): {len(points)} ghost/boundary points")
    if r0 is not None:
        n_flagged = sum(p["flagged"] for p in points)
        print(f"  {n_flagged} flagged as geometrically wrong by the ground-truth overlay")
        reasons = Counter(p["reason"] for p in points if p["flagged"])
        for reason, n in reasons.most_common():
            print(f"    {n:6d}  {reason}")
    print()

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    palette = cm.get_cmap("tab10").colors

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    if not args.no_outlines:
        for p_id, pts in per_patch_interior.items():
            is_query = p_id == args.patch
            draw_patch_outline(ax, pts, color=patch_color(p_id, palette),
                                lw=1.8 if is_query else 0.5,
                                alpha=0.9 if is_query else 0.22,
                                grid_lines=args.outline_grid_lines)

    plot_scatter(ax, points, palette, args.point_size)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    set_equal_aspect(ax, points, per_patch_interior)
    ax.set_title(f"patch {args.patch} ({patch_name}): ghost/boundary color post-SYNC")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8, borderaxespad=0)
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=150, bbox_inches="tight")
        print(f"saved to {args.save}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
