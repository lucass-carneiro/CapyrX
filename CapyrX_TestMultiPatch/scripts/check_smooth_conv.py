#!/usr/bin/env python3
"""
Interpatch smooth-field CONVERGENCE checker -- Step 10 of
llama_patch_system_impl.md.

This is NOT check_smooth_test.py. That script is the single-resolution Channel-1
outer-BC A/B test. This one measures the ORDER at which the interpatch ghost-fill
error falls as the grid is refined, separately on the two interface families:

    R2  cube<->wedge   (one side is the central cartesian patch)
    R1  wedge<->wedge  (both sides are spherical Thornburg wedges)

Mechanism (no thorn change needed -- the machinery is always on)
---------------------------------------------------------------
CapyrX_TestMultiPatch's InterpError group, scheduled AT initial, writes an
analytic field into test_data(interior) and runs ONE interpatch SYNC, so
test_data's interpatch ghosts hold the interpolated fill. We read `u` (group
`test_data`, column 12) directly and recompute the exact analytic field at each
ghost's vertex coords here: error = |u - exact|.

We deliberately do NOT read the thorn's own `interp` field. compute_interp_error
does set interp = |u - exact| everywhere, but the LATER compute_deriv_error does
`SYNC: error`, which re-interpolates interp's ghost zones from its all-zero
interior (the interior error is zero because the interior holds exact samples) --
clobbering every ghost to 0. test_data is synced exactly once and never again, so
its ghosts keep the true interpolated values. Run the par file at N = 16/32/64
(angular = radial = cube ncells, doubled each step) and this script reads the
three `test_data` TSVs with their vertex_coords companions, buckets every
interpatch ghost by family, and fits the order.

Field requirement: the par files use the STANDING WAVE, never the parabola. An
order-4 interpolation reproduces any quadratic exactly (machine zero), which
would give no convergence signal; only a non-polynomial field has a genuine h^4
truncation error.

What counts as an interpatch ghost here
---------------------------------------
Only SINGLE-DONOR FACE ghosts: a cell whose index is outside its reporting
patch's (overlap-widened) interior range in EXACTLY ONE axis, that is inside the
multipatch domain, and whose geometric owner (get_owner_patch) is a different
patch. A cell outside the interior range in >=2 axes is an overset corner/edge
ghost covered by 2-3 donors at once; at order 4 its value legitimately blends
across patches (the Step-8 "informational" set), so it is NOT a clean
single-interface error and is excluded from the order fit -- exactly as the
owner test excludes it.

    family R2  iff  (reporting == cartesian) XOR (owner == cartesian)
    family R1  iff  neither reporting nor owner is the cartesian patch

Order fit
---------
Per family and resolution the grid L2 (RMS) norm is
    E(N) = sqrt( mean_{family face ghosts} error^2 ).
Between two resolutions the order is  log(E_coarse/E_fine) / log(N_fine/N_coarse).
The PRIMARY gate is the FINEST consecutive pair (most asymptotic); the coarser
pair and a least-squares slope over all resolutions are reported as corroboration.

Expected order. CarpetX interpolates with an (interpolation_order + 1)-point
stencil per dimension (interpolate.cxx:691, `wx[order+1]`), i.e. a degree-`order`
polynomial, whose Lagrange truncation error is O(h^(order+1)). So at the mandated
interpolation_order = 4 the nominal interpatch convergence is order 5, NOT 4 --
the impl plan's "converge at interpolation_order" undershoots by one. The gate is
therefore a one-sided FLOOR at interpolation_order + 0.5 (= 4.5): both seams must
reach at least the interpolation accuracy. R2 empirically super-converges (~5.5);
a floor does not penalise that.

Donor availability at the 6 face-centre axis points (design 9.4)
---------------------------------------------------------------
At each (+-R,0,0)/(0,+-R,0)/(0,0,+-R) the sphere only touches the cube face and
the cube-side bracket is one-sided -- the spot most at risk of a dropped donor.
The run uses poison_undefined_values=yes, so a genuinely unfilled interpatch
ghost aborts the run; here we additionally assert, per axis, that at least one
near-axis R2 face ghost exists near r=R and that every R2 face ghost's error is
finite and below --donor-bound (a dropped donor would leave the dirichlet
sentinel ~1138, not a sub-unity error).

Exit codes (match the run_checks.sh contract)
    0  both families reach finest-pair order >= --order-floor AND donors OK
    1  a family's finest-pair order is below the floor, or a donor assert failed
    2  could not run: missing/short inputs, <2 resolutions, or a family empty
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
    find_param_file,
    get_owner_patch,
    load_global_coords,
    parse_par_params,
    parse_patch_system,
    read_tsv_rows,
)

FIELD_COL = 11  # 0-based column of `u` in the group-`test_data` TSV (12th field)


def guess_coords_path(field_path):
    """The vertex_coords companion of a group-`test_data` TSV written next to it."""
    p = Path(field_path)
    name = re.sub(r"^capyrx_testmultipatch-test_data",
                  "coordinatesx-vertex_coords", p.name)
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def load_field(path):
    """key (patch, level, comp, i, j, k) -> u value (column 12)."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) <= FIELD_COL:
            continue
        patch, level, comp, i, j, k = (int(fields[2]), int(fields[3]),
                                       int(fields[4]), int(fields[5]),
                                       int(fields[6]), int(fields[7]))
        table[(patch, level, comp, i, j, k)] = float(fields[FIELD_COL])
    return table


def parse_field_params(par_path):
    """Recover the analytic test field and its parameters from the par file so the
    checker can recompute the exact solution (what the thorn's standing_wave /
    parabola compute at cctk_time = 0)."""
    text = Path(par_path).read_text()

    def num(name, default):
        m = re.search(r"CapyrX_TestMultiPatch::" + name + r"\s*=\s*([0-9.eE+-]+)",
                      text)
        return float(m.group(1)) if m else default

    m = re.search(r'CapyrX_TestMultiPatch::test_data\s*=\s*"([^"]*)"', text)
    kind = m.group(1).strip().lower() if m else "standing wave"
    return {"kind": kind, "A": num("A", 1.0),
            "kx": num("kx", 1.0), "ky": num("ky", 1.0), "kz": num("kz", 1.0)}


def exact_field(fp, x, y, z):
    """The analytic solution at t = 0, matching testmultipatch.cxx."""
    if fp["kind"] == "parabola":
        return x * x + y * y + z * z
    two_pi = 2.0 * math.pi
    return (fp["A"] * math.cos(two_pi * fp["kx"] * x)
            * math.cos(two_pi * fp["ky"] * y) * math.cos(two_pi * fp["kz"] * z))


def axes_outside(patch, i, j, k, angular_extent, radial_extent, cube_extent):
    """Number of axes whose index lies outside this patch's (overlap-widened)
    interior range [0, extent]. 0 = interior, 1 = face ghost (single donor),
    >=2 = corner/edge ghost (multi-donor overset blend)."""
    if patch == CARTESIAN:
        hi = cube_extent  # (ex_i, ex_j, ex_k)
        rng = ((0, hi[0]), (0, hi[1]), (0, hi[2]))
    else:
        rng = ((0, angular_extent), (0, angular_extent), (0, radial_extent))
    n = 0
    for idx, (lo, hhi) in zip((i, j, k), rng):
        if idx < lo or idx > hhi:
            n += 1
    return n


class Family:
    def __init__(self, name):
        self.name = name
        self.sumsq = 0.0
        self.count = 0
        self.max_err = 0.0

    def add(self, err):
        self.sumsq += err * err
        self.count += 1
        if err > self.max_err:
            self.max_err = err

    def rms(self):
        return math.sqrt(self.sumsq / self.count) if self.count else float("nan")


class Resolution:
    def __init__(self, N):
        self.N = N
        self.r1 = Family("R1 wedge<->wedge")
        self.r2 = Family("R2 cube<->wedge")


def resolve_par_path(raw, field_path):
    """The parameter-file path in a TSV header is written however cactus was
    invoked -- absolute under run_checks.sh, but relative to the run directory for
    a hand-launched run. Resolve it: as given, then relative to the TSV's own
    directory, then by basename next to the TSV, then in this thorn's par/ dir."""
    if raw is None:
        return None
    cands = [Path(raw)]
    tsv_dir = Path(field_path).resolve().parent
    cands.append(tsv_dir / raw)
    cands.append(tsv_dir / Path(raw).name)
    cands.append(Path(__file__).resolve().parent.parent / "par" / Path(raw).name)
    for c in cands:
        if c.exists():
            return str(c)
    return None


def process_resolution(field_path, coords_override=None, par_override=None):
    par_path = par_override or resolve_par_path(find_param_file(field_path), field_path)
    if par_path is None:
        raise SystemExit(f"could not find the par file named in {field_path}")
    params = parse_par_params(par_path)
    field = parse_field_params(par_path)
    spherical = (parse_patch_system(par_path) == "llama")
    r0 = params.get("inner_boundary_radius")
    angular = params.get("angular_cells")
    radial = params.get("radial_cells")
    patch_overlap = int(params.get("patch_overlap", 0) or 0)
    if None in (r0, angular, radial):
        raise SystemExit(f"could not recover geometry from {par_path}")
    r0 = float(r0)
    r1_outer = float(params["outer_boundary_radius"])
    angular, radial = int(angular), int(radial)
    cube_ncells = (int(params.get("cartesian_ncells_i", angular) or angular),
                   int(params.get("cartesian_ncells_j", angular) or angular),
                   int(params.get("cartesian_ncells_k", angular) or angular))

    angular_extent = angular + 2 * patch_overlap
    radial_extent = radial + patch_overlap
    cube_extent = tuple(n + 2 * patch_overlap for n in cube_ncells)

    N = cube_ncells[0]  # the doubled resolution knob (angular==radial==cube here)
    res = Resolution(N)

    coords_path = coords_override or guess_coords_path(field_path)
    if coords_path is None or not Path(coords_path).exists():
        raise SystemExit(f"could not find vertex_coords companion for {field_path}")
    coords = load_global_coords(coords_path)
    values = load_field(field_path)

    rel = 1e-9
    donor_bad = []         # non-finite / over-bound R2 face ghosts
    axis_hits = {d: 0 for d in range(6)}  # +x,-x,+y,-y,+z,-z near-axis R2 ghosts
    joined = 0

    for key, u in values.items():
        if key not in coords:
            continue
        joined += 1
        patch, _lvl, _cmp, i, j, k = key
        vx, vy, vz = coords[key]
        r = math.sqrt(vx * vx + vy * vy + vz * vz)

        # must be inside the multipatch domain to have an interpatch donor
        inside_cube = (max(abs(vx), abs(vy), abs(vz)) <= r0 * (1 + rel)
                       if not spherical else r <= r0 * (1 + rel))
        if not (inside_cube or r <= r1_outer * (1 + rel)):
            continue
        owner = get_owner_patch(vx, vy, vz, r0, spherical=spherical)
        if owner == patch:
            continue  # own interior / own overlap, not an interpatch ghost

        nout = axes_outside(patch, i, j, k, angular_extent, radial_extent,
                            cube_extent)
        if nout == 0:
            # inside this patch's index range yet owned by another patch: an
            # overset interior point (r in the cube-corner shell). Not a ghost
            # fill; skip (it is written analytically, not interpolated).
            continue
        if nout >= 2:
            continue  # corner/edge ghost: multi-donor blend, excluded by design

        err = abs(u - exact_field(field, vx, vy, vz))
        is_r2 = (patch == CARTESIAN) != (owner == CARTESIAN)
        fam = res.r2 if is_r2 else res.r1
        fam.add(err)

        if is_r2:
            if not math.isfinite(err) or err > process_resolution.donor_bound:
                donor_bad.append((patch, i, j, k, vx, vy, vz, err, owner))
            # near one of the 6 face-centre axes and near r=R?
            if abs(r - r0) <= 2.0 * (r1_outer - r0) / radial:
                ax = [abs(vx), abs(vy), abs(vz)]
                mx = max(ax)
                others = sum(a for a in ax) - mx
                if mx > 0 and others < 0.25 * mx:  # tight tube around an axis
                    d = ax.index(mx) * 2 + (0 if (vx, vy, vz)[ax.index(mx)] > 0 else 1)
                    axis_hits[d] += 1

    res.joined = joined
    res.donor_bad = donor_bad
    res.axis_hits = axis_hits
    return res


process_resolution.donor_bound = 0.5


def order(E_coarse, N_coarse, E_fine, N_fine):
    if not (E_coarse > 0 and E_fine > 0):
        return float("nan")
    return math.log(E_coarse / E_fine) / math.log(N_fine / N_coarse)


def lsq_slope(points):
    """least-squares slope of log(E) vs log(N); returns -slope (the order)."""
    xs = [math.log(n) for n, e in points if e > 0]
    ys = [math.log(e) for n, e in points if e > 0]
    if len(xs) < 2:
        return float("nan")
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return -num / den if den else float("nan")


AXIS_NAMES = ["+x", "-x", "+y", "-y", "+z", "-z"]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("field", nargs="+",
                    help="group-`test_data` TSVs, one per resolution "
                         "(capyrx_testmultipatch-test_data.it000000.p0000.tsv)")
    ap.add_argument("--order-floor", type=float, default=4.5,
                    help="minimum acceptable finest-pair order per family "
                         "(interpolation_order + 0.5; the interpatch fill error is "
                         "O(h^(order+1)), so order 4 -> nominal 5)")
    ap.add_argument("--donor-bound", type=float, default=0.5,
                    help="max |error| allowed at an R2 face ghost; above this a "
                         "donor is presumed dropped (dirichlet sentinel leaks in)")
    ap.add_argument("--max-examples", type=int, default=10)
    args = ap.parse_args()
    process_resolution.donor_bound = args.donor_bound

    resolutions = [process_resolution(p) for p in args.field]
    resolutions.sort(key=lambda r: r.N)

    print("=" * 78)
    print("Llama smooth-field interpatch CONVERGENCE -- Step 10 (R1/R2)")
    print("=" * 78)
    # run_checks.sh reads this exact label for its row count.
    checked = sum(r.r1.count + r.r2.count for r in resolutions)
    print(f"Rows checked against ground truth: {checked}")
    print(f"resolutions (N): {', '.join(str(r.N) for r in resolutions)}")
    print(f"order gate: finest-pair order >= {args.order_floor}   "
          f"donor-bound: {args.donor_bound}")
    print()

    cannot_run = []
    if len(resolutions) < 2:
        cannot_run.append("need >= 2 resolutions to measure an order")

    for r in resolutions:
        print(f"--- N={r.N} (joined {r.joined} cells) ---")
        for fam in (r.r2, r.r1):
            print(f"  {fam.name:18s}  face-ghosts={fam.count:6d}  "
                  f"L2={fam.rms():.4e}  Linf={fam.max_err:.4e}")
            if fam.count == 0:
                cannot_run.append(f"N={r.N} family '{fam.name}' has zero face ghosts")
        missing_axes = [AXIS_NAMES[d] for d, h in r.axis_hits.items() if h == 0]
        if missing_axes:
            cannot_run.append(f"N={r.N} no near-axis R2 face ghost on axes "
                              f"{','.join(missing_axes)} (face-centre donor coverage)")
        else:
            print(f"  face-centre donor coverage: all 6 axes have a near-axis "
                  f"R2 face ghost near r=R")
        if r.donor_bad:
            print(f"  !! {len(r.donor_bad)} R2 face ghost(s) non-finite or > "
                  f"{args.donor_bound} (dropped donor?):")
            for row in r.donor_bad[: args.max_examples]:
                p, i, j, k, vx, vy, vz, err, ow = row
                print(f"     patch={p} ijk=({i},{j},{k}) r="
                      f"{math.sqrt(vx*vx+vy*vy+vz*vz):.4f} err={err:.3e} "
                      f"owner={ow}({PATCH_NAMES.get(ow)})")
        print()

    if cannot_run:
        print("-" * 78)
        print("RESULT: CANNOT-RUN -- " + "; ".join(cannot_run))
        sys.exit(2)

    # donor failures are a correctness failure, not a could-not-run
    donor_fail = sum(len(r.donor_bad) for r in resolutions)

    print("=" * 78)
    print("convergence order (log ratio across consecutive resolutions)")
    print("=" * 78)
    fail = donor_fail > 0
    finest = {}
    for fam_attr, label in (("r2", "R2 cube<->wedge"), ("r1", "R1 wedge<->wedge")):
        pts = [(r.N, getattr(r, fam_attr).rms()) for r in resolutions]
        print(f"  {label}:")
        orders = []
        for a, b in zip(pts, pts[1:]):
            o = order(a[1], a[0], b[1], b[0])
            orders.append(o)
            print(f"    N {a[0]:>3d}->{b[0]:<3d}  E {a[1]:.4e} -> {b[1]:.4e}  "
                  f"order = {o:.3f}")
        slope = lsq_slope(pts)
        print(f"    least-squares slope over all N: {slope:.3f}")
        finest[fam_attr] = orders[-1]
        if orders[-1] < args.order_floor:
            fail = True
        print()

    print("-" * 78)
    if fail:
        reasons = []
        if donor_fail:
            reasons.append(f"{donor_fail} R2 donor failure(s)")
        for k, lab in (("r2", "R2"), ("r1", "R1")):
            if finest[k] < args.order_floor:
                reasons.append(f"{lab} finest-pair order {finest[k]:.3f} "
                               f"< floor {args.order_floor}")
        print("RESULT: FAIL -- " + "; ".join(reasons))
        sys.exit(1)
    print(f"RESULT: PASS -- R2 order {finest['r2']:.3f}, R1 order "
          f"{finest['r1']:.3f}, both >= {args.order_floor}; "
          f"all 6 face-centre axes have donors.")
    sys.exit(0)


if __name__ == "__main__":
    main()
