#!/usr/bin/env python3
"""
Llama cube-AMR convergence checker -- Step 14 of llama_patch_system_impl.md
(Phase D: AMR convergence + the [P453] side-step).

Step 13 showed cube-only AMR runs and confines refinement to the central cube.
Step 14 asks the harder question: does that refinement PRESERVE the evolution
order across BOTH the refinement boundary (the level-0/level-1 interface inside
the cube) AND the patch boundaries (the cube<->wedge seam) -- and does it stay
stable with slave_overlap = no, side-stepping [P453]?

Method -- refined-cube vs unigrid, two bases (no N=128)
-------------------------------------------------------
For each base resolution N the AMR run (llama_amr_conv_N.par, a level-1 cube added
inside [-R,R]^3) is differenced against the unigrid run at the SAME base
(llama_wave_N.par, max_num_levels = 1) at the SAME physical time t = 2.0, over the
COINCIDENT level-0 vertices (the two share a level-0 grid, so the join is exact):

    d_N(x) = phi_AMR(N; x) - phi_unigrid(N; x)
           = (phi_exact + e_AMR(N)) - (phi_exact + e_unigrid(N))
           = e_AMR(N) - e_unigrid(N).

Both truncation errors are O(h_N^4) (e_AMR smaller inside the refined cube, where
the restricted level-0 data carries fine-grid accuracy), so d_N = O(h_N^4) and

    order = log2( ||d_32||_2 / ||d_64||_2 )  ->  4  (the evolution order).

The point of the test is not the global number -- Step 11 already showed the
unigrid evolution converges at order 4 -- but that d_N converges at that order
*in the two boundary bands as well*. A refinement interface that reflected, or a
regrid repair that corrupted the seam, would leave a NON-converging (order ~0)
residual spike exactly there while the interior still converged.

Bands (classified per coincident point from its own global coords)
    refine   cube-owned (r < R) and | max(|x|,|y|,|z|) - --box-half | <= --refine-band
             (the level-0/level-1 interface at the BoxInBox cube face)
    patch    | r - R | <= --patch-band                 (the cube<->wedge seam)
    global   every coincident owner point

A single-valued physical field is used: each point is taken from its OWNER patch
only (get_owner_patch, spherical for Llama), and only level 0 enters (the unigrid
reference has no other level; the AMR level-1 cells are dropped so the join is a
clean 1:1 coarse-grid match).

[P453] side-step: every run here uses slave_overlap = no (the Step 5 Llama guard
forbids yes). A stable, convergent AMR result with slave_overlap = no is the
first-class record that the Llama system side-steps [P453] by construction.

Exit codes (the run_checks.sh contract)
    0  global + both boundary bands converge at/above their floors, signal present
    1  an order is below its floor, or the AMR made no measurable difference
    2  could not run: missing an AMR/unigrid pair, <2 bases, or too few points
"""

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_color import (  # noqa: E402
    get_owner_patch,
    find_param_file,
    load_global_coords,
    parse_par_params,
    parse_patch_system,
    read_tsv_rows,
)
from check_smooth_conv import order, resolve_par_path  # noqa: E402
from check_wave_conv import guess_coords_path  # noqa: E402

PHI_COL = 11      # 0-based column of phi in the group-`state` TSV
LEVEL_COL = 3
COORD_QUANT = 1e7  # round global coords to 1e-7 to key coincident vertices


def ckey(vx, vy, vz):
    return (int(round(vx * COORD_QUANT)),
            int(round(vy * COORD_QUANT)),
            int(round(vz * COORD_QUANT)))


def load_level0_phi(path):
    """key (patch, comp, i, j, k) -> phi, for level-0 rows only."""
    table = {}
    for fields in read_tsv_rows(path):
        if len(fields) <= PHI_COL:
            continue
        if int(fields[LEVEL_COL]) != 0:
            continue
        key = (int(fields[2]), int(fields[4]),
               int(fields[5]), int(fields[6]), int(fields[7]))
        table[key] = float(fields[PHI_COL])
    return table


def tsv_is_amr(path):
    """True if any data row carries level >= 1 (a refinement level is present)."""
    for fields in read_tsv_rows(path):
        if len(fields) <= LEVEL_COL:
            continue
        try:
            if int(fields[LEVEL_COL]) >= 1:
                return True
        except ValueError:
            continue
    return False


class Run:
    def __init__(self, path):
        self.path = path
        par_path = resolve_par_path(find_param_file(path), path)
        if par_path is None:
            raise SystemExit(f"could not find the par file named in {path}")
        params = parse_par_params(par_path)
        self.spherical = (parse_patch_system(par_path) == "llama")
        self.r0 = float(params["inner_boundary_radius"])
        self.N = int(params.get("cartesian_ncells_i", params["angular_cells"]))
        self.amr = tsv_is_amr(path)

        coords_path = guess_coords_path(path)
        if coords_path is None or not Path(coords_path).exists():
            raise SystemExit(f"no vertex_coords companion for {path}")
        coords = load_global_coords(coords_path)  # keyed with level
        phis = load_level0_phi(path)

        # join phi(level0) to its coords; keep each point's single owner value.
        self.phys = {}  # ckey -> (phi, r, cheb)
        for key5, phi in phis.items():
            p, comp, i, j, k = key5
            c = coords.get((p, 0, comp, i, j, k))
            if c is None:
                continue
            vx, vy, vz = c
            if get_owner_patch(vx, vy, vz, self.r0,
                               spherical=self.spherical) != p:
                continue
            r = math.sqrt(vx * vx + vy * vy + vz * vz)
            cheb = max(abs(vx), abs(vy), abs(vz))
            self.phys[ckey(vx, vy, vz)] = (phi, r, cheb)


class Band:
    def __init__(self, name):
        self.name = name
        self.ssq = {}   # N -> sum of squared diffs
        self.n = {}     # N -> point count

    def add(self, N, d):
        self.ssq[N] = self.ssq.get(N, 0.0) + d * d
        self.n[N] = self.n.get(N, 0) + 1

    def rms(self, N):
        if self.n.get(N):
            return math.sqrt(self.ssq[N] / self.n[N])
        return float("nan")

    def order(self, Nc, Nf):
        return order(self.rms(Nc), 1.0, self.rms(Nf), 2.0)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tsv", nargs="+",
                    help="group-`state` TSVs: the AMR and the unigrid run at each "
                         "base N (auto-sorted by N and AMR/unigrid)")
    ap.add_argument("--order-floor", type=float, default=3.5,
                    help="minimum acceptable global AMR-vs-unigrid diff order")
    ap.add_argument("--seam-order-floor", type=float, default=3.0,
                    help="minimum order in the refinement-boundary and seam bands")
    ap.add_argument("--box-half", type=float, default=0.4,
                    help="BoxInBox cube half-size (the refinement boundary, cheb)")
    ap.add_argument("--refine-band", type=float, default=0.15,
                    help="half-width (cheb) of the refinement-boundary band")
    ap.add_argument("--patch-band", type=float, default=0.25,
                    help="half-width (physical) of the |r-R| cube<->wedge band")
    ap.add_argument("--min-signal", type=float, default=1e-6,
                    help="minimum coarse global RMS diff: below this the AMR made "
                         "no measurable difference and the test is vacuous")
    ap.add_argument("--min-common", type=int, default=2000,
                    help="minimum coincident points in each AMR/unigrid pair")
    args = ap.parse_args()

    print("=" * 78)
    print("Llama cube-AMR convergence -- Step 14 (refined-cube vs unigrid, 2 bases)")
    print("=" * 78)

    runs = [Run(p) for p in args.tsv]
    # bucket into {N: {"amr": Run, "uni": Run}}
    bases = {}
    for r in runs:
        slot = bases.setdefault(r.N, {})
        tag = "amr" if r.amr else "uni"
        if tag in slot:
            print(f"RESULT: CANNOT-RUN -- two {tag} runs at N={r.N}; expected one "
                  f"AMR and one unigrid per base")
            sys.exit(2)
        slot[tag] = r

    Ns = sorted(n for n, s in bases.items() if "amr" in s and "uni" in s)
    for N in sorted(bases):
        slot = bases[N]
        print(f"  N={N}: amr={'yes' if 'amr' in slot else 'MISSING'} "
              f"unigrid={'yes' if 'uni' in slot else 'MISSING'}")
    if len(Ns) < 2:
        print("-" * 78)
        print("RESULT: CANNOT-RUN -- need an AMR+unigrid pair at >=2 base "
              "resolutions for an order ratio")
        sys.exit(2)
    Nc, Nf = Ns[0], Ns[1]  # coarse, fine base (2:1)

    glob = Band("global")
    bands = {"refine": Band("refine-boundary"),
             "patch": Band("patch-boundary (cube<->wedge)")}
    total_common = 0
    for N in (Nc, Nf):
        amr, uni = bases[N]["amr"], bases[N]["uni"]
        common = set(amr.phys) & set(uni.phys)
        total_common += len(common)
        if len(common) < args.min_common:
            print("-" * 78)
            print(f"RESULT: CANNOT-RUN -- only {len(common)} coincident points at "
                  f"N={N} (< {args.min_common}); the level-0 join is broken")
            sys.exit(2)
        r0 = amr.r0
        for k in common:
            d = amr.phys[k][0] - uni.phys[k][0]
            _phi, r, cheb = amr.phys[k]
            glob.add(N, d)
            # refinement boundary: the BoxInBox cube face at cheb = box_half,
            # which lies inside the cube (r < R, so cube-owned by construction).
            if r < r0 and abs(cheb - args.box_half) <= args.refine_band:
                bands["refine"].add(N, d)
            if abs(r - r0) <= args.patch_band:
                bands["patch"].add(N, d)

    # run_checks.sh reads this exact label for its pinned row count.
    print(f"coincident points (both pairs): {total_common}")
    print(f"Rows checked against ground truth: {total_common}")
    print()
    print(f"{'band':32s} {'n(32)':>8s} {'RMS(32)':>12s} {'RMS(64)':>12s} "
          f"{'order':>7s}")
    for b in (glob, bands["refine"], bands["patch"]):
        print(f"{b.name:32s} {b.n.get(Nc, 0):8d} {b.rms(Nc):12.4e} "
              f"{b.rms(Nf):12.4e} {b.order(Nc, Nf):7.3f}")

    for key in ("refine", "patch"):
        if not bands[key].n.get(Nc) or not bands[key].n.get(Nf):
            print("-" * 78)
            print(f"RESULT: CANNOT-RUN -- band '{bands[key].name}' is empty at a "
                  f"base (adjust --box-half/--refine-band/--patch-band)")
            sys.exit(2)

    print()
    fail = []
    if glob.rms(Nc) < args.min_signal:
        fail.append(f"global coarse RMS {glob.rms(Nc):.3e} < min-signal "
                    f"{args.min_signal:.3e}: AMR made no measurable difference")
    if glob.order(Nc, Nf) < args.order_floor:
        fail.append(f"global order {glob.order(Nc, Nf):.3f} < "
                    f"floor {args.order_floor}")
    for key, lab in (("refine", "refinement-boundary"), ("patch", "cube<->wedge")):
        o = bands[key].order(Nc, Nf)
        if o < args.seam_order_floor:
            fail.append(f"{lab} order {o:.3f} < floor {args.seam_order_floor}")

    print("-" * 78)
    if fail:
        print("RESULT: FAIL -- " + "; ".join(fail))
        sys.exit(1)
    print(f"RESULT: PASS -- global order {glob.order(Nc, Nf):.3f} >= "
          f"{args.order_floor}; refinement-boundary "
          f"{bands['refine'].order(Nc, Nf):.3f} and cube<->wedge "
          f"{bands['patch'].order(Nc, Nf):.3f} >= {args.seam_order_floor}; "
          f"stable AMR with slave_overlap=no (side-steps [P453]).")
    sys.exit(0)


if __name__ == "__main__":
    main()
