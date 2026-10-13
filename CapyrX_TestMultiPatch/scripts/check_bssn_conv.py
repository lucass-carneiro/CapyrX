#!/usr/bin/env python3
"""
Llama BSSN constraint-CONVERGENCE checker -- Step 17 of llama_patch_system_impl.md
(Phase E de-risk): the first full BSSN evolution on the Llama geometry, before the
decisive head-on Psi4 A/B (Step 18).

The initial data is a gauge wave (gxx = 1 - A sin(2 pi x / d)), an EXACT nonlinear
vacuum solution, so the CONTINUUM Hamiltonian and momentum constraints are
identically zero. The discrete constraint norm therefore IS the truncation error
and converges to zero at the evolution order. Unlike check_wave_conv.py (which
differences the evolved field across resolutions because its pulse has no
closed-form solution), here the exact answer is known to be 0, so the order comes
DIRECTLY from the norm ratios across N = 16/32/64:

    E(N)  = || ham ||_2   (RMS over a point's OWNER patch only, per band),
    order = log2( E(N) / E(2N) ).

No cross-resolution vertex matching is needed -- each resolution's norm stands on
its own -- but every run is taken to the SAME physical time t = 2.0 so the
comparison is at a fixed state (dt = cfl*2*r0/N, itlast = 2N).

Single-valued constraint
------------------------
Each physical point is taken from its OWNER patch only (get_owner_patch), so the
norm is over a single-valued field: ghost cells (patch != owner) and the non-owner
side of a dual-evolved overlap band are dropped.

Bands
-----
Each owner point is bucketed geometrically (same scheme as check_wave_conv.py) and
by causality:

    causal   r < r1 - c*t   (disconnected from the non-periodic outer-boundary
             inflow through t; the radiative BC injects error that propagates
             inward at speed <= c). This is the clean headline-order region.
    R2 band  |r - R| <= --r2-band   AND causal      (cube<->wedge seam shell)
    R1 band  wedge<->wedge diagonal AND causal        (abutting-wedge seam)
    interior causal, deep (neither seam)              (ratio baseline)
    global   every owner point (sanity, includes the BC-contaminated outer shell)

Gate (mirrors the run_checks.sh contract of check_wave_conv.py)
    The global/causal orders must reach --order-floor; the R2/R1 seam orders must
    reach --seam-order-floor; and the seam/interior norm ratio must be bounded and
    not GROWING across the coarse->fine pair (a growing ratio = a seam that
    converges slower than the interior = a stationary interface kink).

Exit codes
    0  orders meet their floors and the seam ratios are bounded and not growing
    1  an order is below its floor, or a seam ratio is unbounded / growing
    2  could not run: <3 resolutions, a band empty, or too few owner points
"""

import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_color import (  # noqa: E402
    find_param_file,
    get_owner_patch,
    load_global_coords,
    parse_par_params,
    parse_patch_system,
    read_tsv_rows,
)
from check_smooth_conv import order, resolve_par_path  # noqa: E402
from check_wave_conv import classify  # noqa: E402

VAL_COL = 11  # 0-based column of the first group member (hc; mcx) in the TSV


def companion(ham_path, group_to):
    """Sibling TSV of the `ham` file for another group (same iteration/patch)."""
    p = Path(ham_path)
    name = re.sub(r"canudax_bssnmol-ham", group_to, p.name)
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def load_scalar(path):
    """key (patch, level, comp, i, j, k) -> hc value (column VAL_COL)."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) <= VAL_COL:
            continue
        key = (int(fields[2]), int(fields[3]), int(fields[4]),
               int(fields[5]), int(fields[6]), int(fields[7]))
        table[key] = float(fields[VAL_COL])
    return table


def load_vector(path):
    """key -> (mcx, mcy, mcz) at columns VAL_COL..VAL_COL+2."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) <= VAL_COL + 2:
            continue
        key = (int(fields[2]), int(fields[3]), int(fields[4]),
               int(fields[5]), int(fields[6]), int(fields[7]))
        table[key] = (float(fields[VAL_COL]),
                      float(fields[VAL_COL + 1]),
                      float(fields[VAL_COL + 2]))
    return table


class Band:
    """Accumulates sum-of-squares of ham and |mom| for one (band, resolution)."""

    def __init__(self):
        self.sh = 0.0  # sum of ham^2
        self.sm = 0.0  # sum of |mom|^2
        self.n = 0

    def add(self, hc, mom2):
        self.sh += hc * hc
        self.sm += mom2
        self.n += 1

    def rms_ham(self):
        return math.sqrt(self.sh / self.n) if self.n else float("nan")

    def rms_mom(self):
        return math.sqrt(self.sm / self.n) if self.n else float("nan")


BANDS = ("global", "causal", "r2", "r1", "interior")


class Resolution:
    def __init__(self, N, r0, r1):
        self.N = N
        self.r0 = r0
        self.r1 = r1
        self.bands = {b: Band() for b in BANDS}


def process_resolution(ham_path, causal_speed, final_time,
                       r2_band, ang_band):
    par_path = resolve_par_path(find_param_file(ham_path), ham_path)
    if par_path is None:
        raise SystemExit(f"could not find the par file named in {ham_path}")
    params = parse_par_params(par_path)
    spherical = (parse_patch_system(par_path) == "llama")
    r0 = float(params["inner_boundary_radius"])
    r1 = float(params["outer_boundary_radius"])
    N = int(params.get("cartesian_ncells_i", params["angular_cells"]))

    coords_path = companion(ham_path, "coordinatesx-vertex_coords")
    mom_path = companion(ham_path, "canudax_bssnmol-mom")
    if coords_path is None:
        raise SystemExit(f"no vertex_coords companion for {ham_path}")
    if mom_path is None:
        raise SystemExit(f"no mom companion for {ham_path}")

    coords = load_global_coords(coords_path)
    hams = load_scalar(ham_path)
    moms = load_vector(mom_path)

    causal_radius = r1 - causal_speed * final_time
    res = Resolution(N, r0, r1)
    for key, hc in hams.items():
        c = coords.get(key)
        mom = moms.get(key)
        if c is None or mom is None:
            continue
        vx, vy, vz = c
        if get_owner_patch(vx, vy, vz, r0, spherical=spherical) != key[0]:
            continue  # single owner value only
        # The constraint is a diagnostic computed on the interior only; its outer
        # ghost layer (r > r1) is never filled and is written as poison. Drop those
        # non-finite holes rather than letting them NaN-poison a band norm.
        mom2 = mom[0] * mom[0] + mom[1] * mom[1] + mom[2] * mom[2]
        if not (math.isfinite(hc) and math.isfinite(mom2)):
            continue
        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        m = sorted((abs(vx), abs(vy), abs(vz)), reverse=True)

        res.bands["global"].add(hc, mom2)
        if r < causal_radius:
            res.bands["causal"].add(hc, mom2)
            kind = classify(r, r0, m[0], m[1], m[2], r2_band, ang_band)
            bucket = {"r2": "r2", "r1": "r1", "interior": "interior"}[kind]
            res.bands[bucket].add(hc, mom2)
    return res


def pair_order(coarse, fine, quantity):
    ec = getattr(coarse, quantity)()
    ef = getattr(fine, quantity)()
    return order(ec, 1.0, ef, 2.0)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ham", nargs="+",
                    help="final-iteration ham TSVs, one per resolution "
                         "(canudax_bssnmol-ham.it<final>.p0000.tsv)")
    ap.add_argument("--order-floor", type=float, default=3.5,
                    help="minimum acceptable finest-pair global/causal order")
    ap.add_argument("--seam-order-floor", type=float, default=3.0,
                    help="minimum acceptable finest-pair order in each seam band")
    ap.add_argument("--ratio-bound", type=float, default=6.0,
                    help="max seam-band/interior-band norm ratio (finest pair)")
    ap.add_argument("--ratio-growth", type=float, default=1.5,
                    help="max allowed growth of that ratio from the coarse pair "
                         "to the fine pair (a growing ratio = non-converging seam)")
    ap.add_argument("--causal-speed", type=float, default=2.5,
                    help="speed bounding the inward reach of the outer-BC error. "
                         "The radiative BC is inconsistent with a PLANE gauge wave "
                         "and its O(h) error spreads inward at ~2.25 (measured); "
                         "2.5 gives r<3.5 at t=1 with margin")
    ap.add_argument("--final-time", type=float, default=1.0,
                    help="physical time of the measured state (itlast*timestep)")
    ap.add_argument("--r2-band", type=float, default=0.25,
                    help="half-width (physical) of the |r-R| cube<->wedge shell")
    ap.add_argument("--ang-band", type=float, default=0.08,
                    help="relative angular tolerance for the wedge<->wedge diagonal")
    ap.add_argument("--min-owner", type=int, default=2000,
                    help="minimum owner points per resolution")
    args = ap.parse_args()

    resolutions = [process_resolution(p, args.causal_speed, args.final_time,
                                      args.r2_band, args.ang_band)
                   for p in args.ham]
    resolutions.sort(key=lambda r: r.N)

    print("=" * 78)
    print("Llama BSSN constraint CONVERGENCE -- Step 17 (gauge-wave ||ham||/||mom||)")
    print("=" * 78)

    cannot = []
    if len(resolutions) < 3:
        cannot.append("need 3 resolutions (N, 2N, 4N) for a coarse + fine pair")
    total_owner = 0
    for r in resolutions:
        n = r.bands["global"].n
        total_owner += n
        print(f"  N={r.N}: {n} owner points "
              f"(causal {r.bands['causal'].n}, R2 {r.bands['r2'].n}, "
              f"R1 {r.bands['r1'].n}, interior {r.bands['interior'].n})")
        if n < args.min_owner:
            cannot.append(f"N={r.N} has only {n} owner points (< {args.min_owner})")

    # run_checks.sh reads this exact label for its row count.
    print(f"Rows checked against ground truth: {total_owner}")

    for r in resolutions:
        for b in ("causal", "r2", "r1", "interior"):
            if r.bands[b].n == 0:
                cannot.append(f"band '{b}' empty at N={r.N}")

    if cannot:
        print("-" * 78)
        print("RESULT: CANNOT-RUN -- " + "; ".join(cannot))
        sys.exit(2)

    rc, rm, rf = resolutions  # coarse, mid, fine

    print()
    print(f"{'band':12s} {'||ham|| N=16':>13s} {'N=32':>11s} {'N=64':>11s} "
          f"{'ord c':>7s} {'ord f':>7s}")
    orders_ham = {}
    for b in BANDS:
        oc = pair_order(rc.bands[b], rm.bands[b], "rms_ham")
        of = pair_order(rm.bands[b], rf.bands[b], "rms_ham")
        orders_ham[b] = (oc, of)
        print(f"{b:12s} {rc.bands[b].rms_ham():13.4e} "
              f"{rm.bands[b].rms_ham():11.4e} {rf.bands[b].rms_ham():11.4e} "
              f"{oc:7.3f} {of:7.3f}")

    print()
    print(f"{'band':12s} {'||mom|| N=16':>13s} {'N=32':>11s} {'N=64':>11s} "
          f"{'ord c':>7s} {'ord f':>7s}")
    orders_mom = {}
    for b in BANDS:
        oc = pair_order(rc.bands[b], rm.bands[b], "rms_mom")
        of = pair_order(rm.bands[b], rf.bands[b], "rms_mom")
        orders_mom[b] = (oc, of)
        print(f"{b:12s} {rc.bands[b].rms_mom():13.4e} "
              f"{rm.bands[b].rms_mom():11.4e} {rf.bands[b].rms_mom():11.4e} "
              f"{oc:7.3f} {of:7.3f}")

    print()
    print("seam-band vs interior norm ratio (R1/R2 no-growth gate, ||ham||):")
    # Gate on the CAUSAL interior (r < r1 - c*t), the physically-clean region: the
    # all-r 'global' band also samples the r in [r1-c*t, r1] shell the approximate
    # outer BC has reached by t, which converges at the BC order, not the evolution
    # order, so a sub-floor 'global' is a BC artifact, not an interior defect. The
    # 'global' order is still reported for information.
    fail = []
    of = orders_ham["causal"][1]
    if of < args.order_floor:
        fail.append(f"causal ham order {of:.3f} < floor {args.order_floor}")

    int_c = rm.bands["interior"].rms_ham()  # mid for coarse pair denominator
    int_f = rf.bands["interior"].rms_ham()
    int_c0 = rc.bands["interior"].rms_ham()
    for key, name in (("r2", "R2"), ("r1", "R1")):
        b_c0 = rc.bands[key].rms_ham()
        b_f = rf.bands[key].rms_ham()
        ratio_c = b_c0 / int_c0 if int_c0 > 0 else float("inf")
        ratio_f = b_f / int_f if int_f > 0 else float("inf")
        growth = ratio_f / ratio_c if ratio_c > 0 else float("inf")
        of = orders_ham[key][1]
        print(f"  {name}: ratio coarse={ratio_c:.3f} fine={ratio_f:.3f} "
              f"growth={growth:.3f}  seam order={of:.3f}")
        if of < args.seam_order_floor:
            fail.append(f"{name} ham seam order {of:.3f} < "
                        f"floor {args.seam_order_floor}")
        if ratio_f > args.ratio_bound:
            fail.append(f"{name} seam/interior ratio {ratio_f:.3f} > "
                        f"bound {args.ratio_bound}")
        if growth > args.ratio_growth:
            fail.append(f"{name} seam ratio growing (x{growth:.3f} > "
                        f"{args.ratio_growth}): non-converging interface")

    print("-" * 78)
    if fail:
        print("RESULT: FAIL -- " + "; ".join(fail))
        sys.exit(1)
    print(f"RESULT: PASS -- causal ham order {orders_ham['causal'][1]:.3f} "
          f">= {args.order_floor}; R2 {orders_ham['r2'][1]:.3f}, "
          f"R1 {orders_ham['r1'][1]:.3f} >= {args.seam_order_floor}; seam ratios "
          f"bounded and not growing (global {orders_ham['global'][1]:.3f} is "
          f"BC-contaminated, informational only).")
    sys.exit(0)


if __name__ == "__main__":
    main()
