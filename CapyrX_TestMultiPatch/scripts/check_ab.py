#!/usr/bin/env python3
"""
Llama vs Thornburg06 A/B checker -- Step 12 of llama_patch_system_impl.md
(Phase C: isolate the central-cube COUPLING).

Both codes evolve the SAME six Thornburg wedges (same R, outer, angular/radial
cells, patch_overlap, interpolation order, ghost width, outer BC, Gaussian-shell
initial data and PINNED time step). The only structural difference is the wedge
inner radial face at r = R:

    Llama        inner face is interpatch, fed by the central cube,
    Thornburg06  inner face is a physical outer boundary (no cube; a hole at r0).

So phi_Llama - phi_Thornburg06 at a shared wedge vertex isolates exactly what the
cube coupling changes relative to a plain inner BC.

Causal-band gate
----------------
The Gaussian shell (r_c = 2.5) splits into in/out-going waves; the inward half
reaches r0, where the difference between the two inner-face treatments is born
and then propagates OUTWARD at the wave speed c = 1. In the continuum, at
t_final = 2.0 nothing from r0 can have travelled past the causal front
r0 + c*t_final = 3.0, so the A/B difference is identically zero beyond it. The
numerical scheme (wide FD stencil x RK substeps) spreads a small precursor past
that front, decaying steeply with radius. Each shared wedge vertex is bucketed:

    inner band  r0 <= r <= --inner-hi          (deep in the domain of dep. of r0)
    causal cone --inner-hi < r < --outer-lo    (the difference in transit)
    far band    --outer-lo <= r <= r1          (beyond the causal front; reported
                                                isolated)

The A/B difference between a reflecting inner BC and the cube coupling is a
PHYSICAL, O(1) difference inside the cone -- it does NOT converge away. So the
test is not "the difference converges" but:

  * inner band: || phi_L - phi_T ||_2 is O(1) and plateaus under refinement ->
    the cube coupling IS exercised (teeth; a near-zero inner band = vacuous A/B);
  * far band: beyond the continuum causal front the difference is a negligible
    residual (--outer-max) that is STILL decreasing with resolution
    (--outer-order-floor) -> NO cube-only signal leaks past the causal horizon.
    The causal cone is reported for context but not gated (it legitimately holds
    the O(1) difference in transit).

A far band that is large, or that does NOT shrink with resolution, is a spurious
cube-only signal and the defect this test guards against.

Single-valued field
-------------------
Each physical point is taken from its OWNER wedge only (get_owner_patch), so the
field is single-valued; ghosts and the non-owner overlap side are dropped. The
two systems number their wedges differently (Llama: cube=0, wedges 1..6;
Thornburg06: wedges 0..5), so points are joined by rounded GLOBAL COORDINATE,
not by patch index. Only wedge vertices present in BOTH systems enter the A/B
(the Llama cube has no Thornburg06 counterpart).

Exit codes (match the run_checks.sh contract)
    0  inner band exercises the coupling AND the far band is isolated (negligible
       and decreasing)
    1  far band holds a cube-only signal (too large, or not decreasing), or the
       inner band is too small (coupling not exercised)
    2  could not run: <2 resolutions, missing companion coords, or a band empty
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
)
from check_smooth_conv import order, resolve_par_path  # noqa: E402
from check_wave_conv import ckey, guess_coords_path, load_phi  # noqa: E402


def detect_system(par_path):
    """"llama" or "thornburg06" from the raw patch_system keyword. (check_color's
    parse_patch_system only distinguishes llama vs cubed_sphere, so the A/B reads
    the keyword directly.)"""
    m = re.search(r'CapyrX_MultiPatch::patch_system\s*=\s*"([^"]*)"',
                  Path(par_path).read_text())
    if m is None:
        raise SystemExit(f"no patch_system in {par_path}")
    kw = m.group(1).strip().lower()
    if kw == "llama":
        return "llama"
    if kw == "thornburg06":
        return "thornburg06"
    raise SystemExit(f"check_ab.py only compares Llama vs Thornburg06, got '{kw}'")


def iteration_of(path):
    m = re.search(r"\.it(\d+)\.", Path(path).name)
    return int(m.group(1)) if m else -1


class Run:
    def __init__(self, system, N, it, field_path):
        self.system = system      # "llama" or "thornburg06"
        self.N = N
        self.it = it
        self.field_path = field_path
        self.phys = {}            # ckey -> (phi, r)


def load_run(field_path):
    par_path = resolve_par_path(find_param_file(field_path), field_path)
    if par_path is None:
        raise SystemExit(f"could not find the par file named in {field_path}")
    params = parse_par_params(par_path)
    system = detect_system(par_path)
    r0 = float(params["inner_boundary_radius"])
    N = int(params.get("cartesian_ncells_i", params["angular_cells"]))

    coords_path = guess_coords_path(field_path)
    if coords_path is None or not Path(coords_path).exists():
        raise SystemExit(f"no vertex_coords companion for {field_path}")
    coords = load_global_coords(coords_path)
    phis = load_phi(field_path)

    # Llama numbers wedges 1..6 (cube = 0); Thornburg06 numbers them 0..5. The
    # owner classifier returns the Llama numbering, so offset Thornburg06 by -1.
    offset = 0 if system == "llama" else 1

    run = Run(system, N, iteration_of(field_path), field_path)
    for key, phi in phis.items():
        c = coords.get(key)
        if c is None:
            continue
        vx, vy, vz = c
        owner = get_owner_patch(vx, vy, vz, r0, spherical=True)
        if owner == CARTESIAN:
            continue  # cube points have no Thornburg06 counterpart
        if owner - offset != key[0]:
            continue  # keep only each wedge point's single owner value
        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        run.phys[ckey(vx, vy, vz)] = (phi, r)
    run.r0 = r0
    return run


class Band:
    def __init__(self, name):
        self.name = name
        self.ss = 0.0
        self.n = 0

    def add(self, d):
        self.ss += d * d
        self.n += 1

    def l2(self):
        return math.sqrt(self.ss / self.n) if self.n else float("nan")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("field", nargs="+",
                    help="group-`state` TSVs: the final-iteration file for each "
                         "Llama and Thornburg06 resolution")
    ap.add_argument("--inner-hi", type=float, default=1.6,
                    help="outer radius of the inner (coupling) band")
    ap.add_argument("--outer-lo", type=float, default=3.25,
                    help="inner radius of the far (causally isolated) band -- must "
                         "sit BEYOND the continuum causal front r0 + c*t_final "
                         "(= 3.0 here) by a margin for numerical dispersion")
    ap.add_argument("--inner-min", type=float, default=1e-3,
                    help="minimum finest-pair inner-band ||diff|| (test teeth: "
                         "the coupling must actually be exercised)")
    ap.add_argument("--outer-max", type=float, default=1e-6,
                    help="far-band ||diff|| at the finest resolution must be below "
                         "this (a negligible residual, not a cube-only signal)")
    ap.add_argument("--outer-order-floor", type=float, default=1.5,
                    help="far-band ||diff|| must also be decreasing at least this "
                         "fast (a persistent cube-only signal would not converge)")
    ap.add_argument("--min-common", type=int, default=1000,
                    help="minimum joined wedge vertices per resolution")
    args = ap.parse_args()

    # Load every TSV, keep the highest-iteration file per (system, N).
    runs = {}
    for p in args.field:
        r = load_run(p)
        prev = runs.get((r.system, r.N))
        if prev is None or r.it > prev.it:
            runs[(r.system, r.N)] = r

    systems = sorted({k[0] for k in runs})
    Ns = sorted({k[1] for k in runs})

    print("=" * 78)
    print("Llama vs Thornburg06 A/B -- Step 12 (isolate central-cube coupling)")
    print("=" * 78)
    for (sys_, N), r in sorted(runs.items()):
        print(f"  {sys_:12s} N={N:3d} it={r.it:4d}: {len(r.phys)} wedge owner points")

    cannot = []
    if "llama" not in systems or "thornburg06" not in systems:
        cannot.append("need both a Llama and a Thornburg06 run")
    paired_Ns = [N for N in Ns
                 if ("llama", N) in runs and ("thornburg06", N) in runs]
    if len(paired_Ns) < 2:
        cannot.append("need >=2 resolutions present for BOTH systems")
    if cannot:
        print("-" * 78)
        print("RESULT: CANNOT-RUN -- " + "; ".join(cannot))
        sys.exit(2)

    r0 = runs[("llama", paired_Ns[0])].r0
    r1_ = None  # outer radius from the par geometry (outer_boundary_radius)
    # read it from any par
    pp = resolve_par_path(find_param_file(runs[("llama", paired_Ns[0])].field_path),
                          runs[("llama", paired_Ns[0])].field_path)
    r1_ = float(parse_par_params(pp)["outer_boundary_radius"])

    print()
    print(f"geometry: r0={r0}  r1={r1_}   "
          f"inner [{r0}, {args.inner_hi}]   "
          f"cone ({args.inner_hi}, {args.outer_lo})   "
          f"far [{args.outer_lo}, {r1_}]")
    print()
    print(f"{'N':>4s} {'joined':>9s} {'inner L2':>12s} {'cone L2':>12s} "
          f"{'far L2':>12s} {'far pts':>9s}")

    inner_l2 = {}
    far_l2 = {}
    total_joined = 0
    for N in paired_Ns:
        L = runs[("llama", N)].phys
        T = runs[("thornburg06", N)].phys
        common = set(L) & set(T)
        total_joined += len(common)
        inner = Band("inner")
        cone = Band("cone")
        far = Band("far")
        for k in common:
            (pl, r) = L[k]
            (pt, _) = T[k]
            d = pl - pt
            if r0 - 1e-9 <= r <= args.inner_hi:
                inner.add(d)
            elif r < args.outer_lo:
                cone.add(d)
            elif r <= r1_ + 1e-9:
                far.add(d)
        inner_l2[N] = (inner.l2(), inner.n)
        far_l2[N] = (far.l2(), far.n)
        print(f"{N:4d} {len(common):9d} {inner.l2():12.4e} {cone.l2():12.4e} "
              f"{far.l2():12.4e} {far.n:9d}")
        if len(common) < args.min_common:
            print("-" * 78)
            print(f"RESULT: CANNOT-RUN -- only {len(common)} joined points at N={N}")
            sys.exit(2)
        if inner.n == 0 or far.n == 0:
            print("-" * 78)
            print(f"RESULT: CANNOT-RUN -- the inner or far band is empty at N={N} "
                  f"(adjust --inner-hi/--outer-lo or the geometry)")
            sys.exit(2)

    # run_checks.sh reads this exact label for its row count.
    print(f"Rows checked against ground truth: {total_joined}")

    Nf = paired_Ns[-1]
    Nc = paired_Ns[-2]
    inner_order = order(inner_l2[Nc][0], 1.0, inner_l2[Nf][0], 2.0)
    far_order = order(far_l2[Nc][0], 1.0, far_l2[Nf][0], 2.0)

    print()
    print(f"finest pair N={Nc}->{Nf}:")
    print(f"  inner band L2: {inner_l2[Nc][0]:.4e} -> {inner_l2[Nf][0]:.4e} "
          f"(order {inner_order:.3f})  [cube coupling: O(1), plateauing]")
    print(f"  far   band L2: {far_l2[Nc][0]:.4e} -> {far_l2[Nf][0]:.4e} "
          f"(order {far_order:.3f})  [isolated: negligible + decreasing]")

    fail = []
    # Teeth: the coupling must be exercised in the inner band.
    if inner_l2[Nf][0] < args.inner_min:
        fail.append(f"inner-band ||diff|| {inner_l2[Nf][0]:.3e} < teeth floor "
                    f"{args.inner_min:.1e}: cube coupling not exercised")
    # Isolation: beyond the continuum causal front the difference must be a
    # negligible residual (tiny) AND still decreasing (a persistent cube-only
    # signal would be neither).
    far_fine = far_l2[Nf][0]
    if far_fine > args.outer_max:
        fail.append(f"far-band ||diff|| {far_fine:.3e} > {args.outer_max:.1e}: "
                    f"a cube-only signal leaked past the causal front")
    if far_order < args.outer_order_floor:
        fail.append(f"far-band order {far_order:.3f} < {args.outer_order_floor}: "
                    f"far-field difference is not decreasing (persistent signal)")

    print("-" * 78)
    if fail:
        print("RESULT: FAIL -- " + "; ".join(fail))
        sys.exit(1)
    print(f"RESULT: PASS -- inner band O(1) ({inner_l2[Nf][0]:.3e}, cube coupling "
          f"exercised); far band isolated ({far_fine:.3e}, order {far_order:.3f}) "
          f"beyond the causal front -- no cube-only signal.")
    sys.exit(0)


if __name__ == "__main__":
    main()
