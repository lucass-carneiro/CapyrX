#!/usr/bin/env python3
"""
Llama WaveToy interpatch SELF-CONVERGENCE checker -- Step 11 of
llama_patch_system_impl.md (Phase C, first evolution on the Llama geometry).

Unlike check_smooth_conv.py (a single interpatch SYNC against an analytic field),
this reads the EVOLVED scalar-wave field phi after a Gaussian pulse has crossed
the cube<->wedge (R2) and wedge<->wedge (R1) seams, and measures convergence by
RICHARDSON SELF-CONVERGENCE across N = 16/32/64 (angular = radial = cube ncells,
doubled). The pulse has no closed-form solution on this geometry, so there is no
exact field to difference against; instead the order is read from

    E_pair(N) = || phi_N - phi_2N ||_2   (over coincident vertices),
    order     = log2( E(16,32) / E(32,64) ).

Because every resolution runs to the SAME physical time t = 2.0 (dt = cfl*(2/N),
itlast = 4N) and the cells double 2:1, every coarse NOMINAL-INTERIOR vertex
coincides with a fine one and phi can be differenced with no interpolation.

Single-valued physical field
----------------------------
Each physical point is taken from its OWNER patch only (get_owner_patch), so the
field is single-valued: ghost cells (patch != owner) and the non-owner side of a
dual-evolved overlap band are dropped. Overlap/ghost vertices sit at
resolution-dependent positions (the overlap is a fixed 2 cells, so its physical
width halves as N doubles) and simply fail to join across resolutions; only the
nominal-interior owner vertices -- which DO coincide -- enter the fit.

Seam gate (R1 / no growing interface noise)
-------------------------------------------
The primary R1 guard is the seam-band vs interior-band residual ratio. Each
coincident point is bucketed geometrically:

    R2 band   |r - R| <= --r2-band                      (cube<->wedge seam shell)
    R1 band   r > R + --r2-band  AND  near a wedge<->wedge diagonal
              (the top two |coords| nearly equal, the third clearly smaller)
    interior  everything else (deep cube + mid-wedge away from the diagonals)

A stationary interface kink shows up as a seam band that converges SLOWER than
the interior, i.e. a seam/interior residual ratio that GROWS from the coarse pair
(16->32) to the fine pair (32->64). The gate therefore requires, at the finest
pair, both seam families to (a) converge at or above --seam-order-floor and
(b) have a seam/interior ratio that is bounded and not growing across the two
pairs (within --ratio-growth). The global order must reach --order-floor.

Exit codes (match the run_checks.sh contract)
    0  global + both seam families converge and the seam ratios are bounded
    1  an order is below its floor, or a seam ratio is unbounded / growing
    2  could not run: <3 resolutions, too few coincident points, or a family empty
"""

import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_color import (  # noqa: E402
    CARTESIAN,
    find_param_file,
    get_owner_patch,
    load_global_coords,
    parse_par_params,
    parse_patch_system,
    read_tsv_rows,
)
from check_smooth_conv import order, resolve_par_path  # noqa: E402

PHI_COL = 11  # 0-based column of `phi` in the group-`state` TSV (12th field)
COORD_QUANT = 1e7  # round global coords to 1e-7 to key coincident vertices


def guess_coords_path(field_path):
    """The vertex_coords companion of a group-`state` TSV written next to it."""
    p = Path(field_path)
    name = re.sub(r"^capyrx_wavetoy-state",
                  "coordinatesx-vertex_coords", p.name)
    if name == p.name:
        return None
    candidate = p.with_name(name)
    return candidate if candidate.exists() else None


def load_phi(path):
    """key (patch, level, comp, i, j, k) -> phi value (column 12)."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) <= PHI_COL:
            continue
        key = (int(fields[2]), int(fields[3]), int(fields[4]),
               int(fields[5]), int(fields[6]), int(fields[7]))
        table[key] = float(fields[PHI_COL])
    return table


def ckey(vx, vy, vz):
    return (int(round(vx * COORD_QUANT)),
            int(round(vy * COORD_QUANT)),
            int(round(vz * COORD_QUANT)))


class Resolution:
    def __init__(self, N, r0):
        self.N = N
        self.r0 = r0
        self.phys = {}  # ckey -> (phi, r, m1, m2, m3) at owner points


def process_resolution(field_path):
    par_path = resolve_par_path(find_param_file(field_path), field_path)
    if par_path is None:
        raise SystemExit(f"could not find the par file named in {field_path}")
    params = parse_par_params(par_path)
    spherical = (parse_patch_system(par_path) == "llama")
    r0 = float(params["inner_boundary_radius"])
    N = int(params.get("cartesian_ncells_i", params["angular_cells"]))

    coords_path = guess_coords_path(field_path)
    if coords_path is None or not Path(coords_path).exists():
        raise SystemExit(f"no vertex_coords companion for {field_path}")
    coords = load_global_coords(coords_path)
    phis = load_phi(field_path)

    res = Resolution(N, r0)
    for key, phi in phis.items():
        c = coords.get(key)
        if c is None:
            continue
        vx, vy, vz = c
        if get_owner_patch(vx, vy, vz, r0, spherical=spherical) != key[0]:
            continue  # keep only each point's single owner value
        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        m = sorted((abs(vx), abs(vy), abs(vz)), reverse=True)
        res.phys[ckey(vx, vy, vz)] = (phi, r, m[0], m[1], m[2])
    return res


class Band:
    def __init__(self, name):
        self.name = name
        self.sc = 0.0  # sum of coarse-pair squared residuals
        self.sf = 0.0  # sum of fine-pair squared residuals
        self.n = 0

    def add(self, dc, df):
        self.sc += dc * dc
        self.sf += df * df
        self.n += 1

    def l2c(self):
        return math.sqrt(self.sc / self.n) if self.n else float("nan")

    def l2f(self):
        return math.sqrt(self.sf / self.n) if self.n else float("nan")

    def order(self):
        return order(self.l2c(), 1.0, self.l2f(), 2.0)


def classify(r, r0, m1, m2, m3, r2_band, ang):
    """'r2', 'r1', or 'interior' for a coincident point."""
    if abs(r - r0) <= r2_band:
        return "r2"
    if r > r0 + r2_band and m1 > 0:
        # near a wedge<->wedge diagonal: top two |coords| nearly equal, third
        # clearly smaller (a face-to-face seam, not a triple corner).
        if (m1 - m2) <= ang * m1 and (m2 - m3) > ang * m1:
            return "r1"
    return "interior"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("field", nargs="+",
                    help="group-`state` TSVs, one per resolution "
                         "(capyrx_wavetoy-state.it<final>.p0000.tsv)")
    ap.add_argument("--order-floor", type=float, default=3.5,
                    help="minimum acceptable global finest-pair self-conv order")
    ap.add_argument("--seam-order-floor", type=float, default=3.0,
                    help="minimum acceptable finest-pair order in each seam band")
    ap.add_argument("--ratio-bound", type=float, default=3.0,
                    help="max seam-band/interior-band residual ratio (finest pair)")
    ap.add_argument("--ratio-growth", type=float, default=1.5,
                    help="max allowed growth of that ratio from the coarse pair "
                         "to the fine pair (a growing ratio = non-converging seam)")
    ap.add_argument("--r2-band", type=float, default=0.25,
                    help="half-width (physical) of the |r-R| cube<->wedge shell")
    ap.add_argument("--ang-band", type=float, default=0.08,
                    help="relative angular tolerance for the wedge<->wedge diagonal")
    ap.add_argument("--min-common", type=int, default=2000,
                    help="minimum coincident points across all resolutions")
    args = ap.parse_args()

    resolutions = [process_resolution(p) for p in args.field]
    resolutions.sort(key=lambda r: r.N)

    print("=" * 78)
    print("Llama WaveToy interpatch SELF-CONVERGENCE -- Step 11 (R1/R2 under evolution)")
    print("=" * 78)

    cannot = []
    if len(resolutions) < 3:
        cannot.append("need 3 resolutions (N, 2N, 4N) for a coarse + fine pair")
    for r in resolutions:
        print(f"  N={r.N}: {len(r.phys)} owner points")

    if cannot:
        print("-" * 78)
        print("RESULT: CANNOT-RUN -- " + "; ".join(cannot))
        sys.exit(2)

    rc, rm, rf = resolutions  # coarse, mid, fine
    r0 = rc.r0
    common = set(rc.phys) & set(rm.phys) & set(rf.phys)
    print(f"coincident points across all three: {len(common)}")
    # run_checks.sh reads this exact label for its row count.
    print(f"Rows checked against ground truth: {len(common)}")
    if len(common) < args.min_common:
        print("-" * 78)
        print(f"RESULT: CANNOT-RUN -- only {len(common)} coincident points "
              f"(< {args.min_common}); the 2:1 vertex join is broken")
        sys.exit(2)

    glob = Band("global")
    bands = {"r2": Band("R2 cube<->wedge"),
             "r1": Band("R1 wedge<->wedge"),
             "interior": Band("interior")}
    for k in common:
        pc = rc.phys[k][0]
        pm = rm.phys[k][0]
        pf = rf.phys[k][0]
        dc = pc - pm
        df = pm - pf
        glob.add(dc, df)
        _phi, r, m1, m2, m3 = rf.phys[k]
        bands[classify(r, r0, m1, m2, m3, args.r2_band, args.ang_band)].add(dc, df)

    print()
    print(f"{'band':20s} {'points':>8s} {'E(16,32)':>12s} {'E(32,64)':>12s} "
          f"{'order':>7s}")
    for b in (glob, bands["r2"], bands["r1"], bands["interior"]):
        print(f"{b.name:20s} {b.n:8d} {b.l2c():12.4e} {b.l2f():12.4e} "
              f"{b.order():7.3f}")

    for key in ("r2", "r1", "interior"):
        if bands[key].n == 0:
            print("-" * 78)
            print(f"RESULT: CANNOT-RUN -- band '{bands[key].name}' is empty "
                  f"(adjust --r2-band/--ang-band or the geometry)")
            sys.exit(2)

    interior_f = bands["interior"].l2f()
    interior_c = bands["interior"].l2c()
    print()
    print("seam-band vs interior residual ratio (R1 no-growth gate):")
    fail = []
    if glob.order() < args.order_floor:
        fail.append(f"global order {glob.order():.3f} < floor {args.order_floor}")
    for key, lab in (("r2", "R2"), ("r1", "R1")):
        b = bands[key]
        ratio_f = b.l2f() / interior_f if interior_f > 0 else float("inf")
        ratio_c = b.l2c() / interior_c if interior_c > 0 else float("inf")
        growth = ratio_f / ratio_c if ratio_c > 0 else float("inf")
        print(f"  {lab}: ratio coarse={ratio_c:.3f} fine={ratio_f:.3f} "
              f"growth={growth:.3f}  seam order={b.order():.3f}")
        if b.order() < args.seam_order_floor:
            fail.append(f"{lab} seam order {b.order():.3f} < "
                        f"floor {args.seam_order_floor}")
        if ratio_f > args.ratio_bound:
            fail.append(f"{lab} seam/interior ratio {ratio_f:.3f} > "
                        f"bound {args.ratio_bound}")
        if growth > args.ratio_growth:
            fail.append(f"{lab} seam ratio growing (x{growth:.3f} > "
                        f"{args.ratio_growth}): non-converging interface")

    print("-" * 78)
    if fail:
        print("RESULT: FAIL -- " + "; ".join(fail))
        sys.exit(1)
    print(f"RESULT: PASS -- global order {glob.order():.3f} >= {args.order_floor}; "
          f"R2 order {bands['r2'].order():.3f}, R1 order {bands['r1'].order():.3f} "
          f">= {args.seam_order_floor}; seam ratios bounded and not growing.")
    sys.exit(0)


if __name__ == "__main__":
    main()
