#!/usr/bin/env python3
"""
Plane-wave A/B vs real Llama -- Step 16 of llama_patch_system_impl.md (Phase E).

Both codes evolve the SAME exact travelling plane wave
    phi = A cos(2*pi*(k.(x-x0) + |k|*(t-t0)))
(a Minkowski solution) on a matched central-cube + 6-wedge grid and store the
pointwise error = numerical - exact at the comparison time t = 1.0:

  CapyrX  (this checkout)  -> CapyrX_WaveToy `error` group, 3d TSV + vertex_coords
  real Llama (AxiDil exe)  -> LlamaWaveToy  `error` GF, CarpetIOASCII 3d .asc with
                             use_grid_coordinates=yes (global Cartesian columns)

Because each code is differenced against its OWN analytic solution, the error
norms stand alone per resolution -- no cross-code or cross-resolution node
matching is needed. The order is the plain Richardson ratio

    order(N->2N) = log2( ||error||_2(N) / ||error||_2(2N) ).

Single-valued field (owner only)
--------------------------------
Each physical point is kept from its OWNER patch only, so overlap/donor points are
not counted twice:
    r <= R            -> cube owner
    r >  R            -> wedge owner = dominant |coord| axis and sign
CapyrX uses get_owner_patch (spherical); the Llama side applies the same rule to
the global coords, dropping cube-file points with r>R and wedge-file points not in
their own angular sector.

Comparison region
-----------------
Only r < --r-max (default 2.5) enters the fit. At t=1.0 the outer-boundary signal
(speed 1 from r=4) has reached r=3, but its numerical precursor leaks a little
further in, and the two codes' native BCs (CapyrX linear-extrapolation vs Llama
radiative) contaminate to slightly different depths -- the radiative BC, fed a
plane wave it cannot treat as outgoing, reaches deeper. r-max=2.5 keeps a 0.5
margin inside the r=3 front so both codes show their clean FD order; at r-max=3.0
the Llama interior order drops to ~3.5 purely from this BC precursor (CapyrX,
contaminating less, still looks order 4). Points are bucketed into the
cube<->wedge seam shell (|r-R| <= --seam-band) vs the interior.

Green gate (matches the run_checks.sh contract)
    0  both codes reach the FD order (>= --order-floor) in interior and seam, AND
       the CapyrX seam/interior error ratio is bounded and not growing
       (no CapyrX-only interface signal)
    1  an order is below floor, or the CapyrX seam ratio is unbounded / growing
    2  could not run: missing resolutions/files, or a band empty
"""

import argparse
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_color import (  # noqa: E402
    get_owner_patch,
    load_global_coords,
    read_tsv_rows,
)

PHI_ERR_COL = 11  # phi_error is the 1st member of group `error` (col 12, 0-based 11)


def order(coarse, fine):
    if not (coarse > 0 and fine > 0):
        return float("nan")
    return math.log2(coarse / fine)


def dir_N(path):
    m = re.search(r"(\d+)\D*$", Path(path).name)
    if not m:
        raise SystemExit(f"cannot parse resolution N from {path}")
    return int(m.group(1))


def owner_is_cube(x, y, z, R):
    return (x * x + y * y + z * z) <= R * R


def band_of(r, R, seam):
    return "seam" if abs(r - R) <= seam else "interior"


class Bands:
    """L2 accumulators for interior / seam / (interior+seam)=global."""

    def __init__(self):
        self.s = {"interior": 0.0, "seam": 0.0}
        self.n = {"interior": 0, "seam": 0}

    def add(self, band, e):
        self.s[band] += e * e
        self.n[band] += 1

    def l2(self, band):
        return math.sqrt(self.s[band] / self.n[band]) if self.n[band] else float("nan")

    def l2_global(self):
        tot = self.s["interior"] + self.s["seam"]
        cnt = self.n["interior"] + self.n["seam"]
        return math.sqrt(tot / cnt) if cnt else float("nan")

    def count(self):
        return self.n["interior"] + self.n["seam"]


def capyrx_bands(run_dir, R, r_max, seam):
    """Owner-only banded error norms from a CapyrX run directory."""
    run_dir = Path(run_dir)
    errs = sorted(run_dir.glob("capyrx_wavetoy-error.it*.p0000.tsv"))
    if not errs:
        raise SystemExit(f"no CapyrX error TSV in {run_dir}")
    it = max(int(re.search(r"\.it(\d+)\.", p.name).group(1)) for p in errs)
    err_path = run_dir / f"capyrx_wavetoy-error.it{it:06d}.p0000.tsv"
    crd_path = run_dir / f"coordinatesx-vertex_coords.it{it:06d}.p0000.tsv"
    coords = load_global_coords(str(crd_path))
    b = Bands()
    for f in read_tsv_rows(str(err_path)):
        if len(f) <= PHI_ERR_COL:
            continue
        key = (int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]))
        c = coords.get(key)
        if c is None:
            continue
        vx, vy, vz = c
        if get_owner_patch(vx, vy, vz, R, spherical=True) != key[0]:
            continue
        r = math.sqrt(vx * vx + vy * vy + vz * vz)
        if r >= r_max:
            continue
        b.add(band_of(r, R, seam), abs(float(f[PHI_ERR_COL])))
    return b


def _asc_rows(path, want_it):
    """Yield (x, y, z, value) for rows of a CarpetIOASCII compact 3d .asc at
    iteration want_it. Row: it <ix iy iz> time <x y z> value."""
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line[0] == "#":
                continue
            parts = line.split()
            # compact: it ix iy iz time x y z value  -> 9 fields
            if len(parts) < 9:
                continue
            if int(parts[0]) != want_it:
                continue
            x, y, z = float(parts[5]), float(parts[6]), float(parts[7])
            yield x, y, z, float(parts[8])


def _wedge_axis(path, want_it, R):
    """Majority (axis, sign) of a wedge file's clearly-interior points."""
    votes = {}
    for x, y, z, _ in _asc_rows(path, want_it):
        r = math.sqrt(x * x + y * y + z * z)
        if r <= R * 1.3:
            continue
        a = (abs(x), abs(y), abs(z))
        ax = a.index(max(a))
        sign = 1 if (x, y, z)[ax] >= 0 else -1
        votes[(ax, sign)] = votes.get((ax, sign), 0) + 1
    return max(votes, key=votes.get) if votes else None


def llama_bands(run_dir, R, r_max, seam):
    """Owner-only banded error norms from a Llama run directory (per-patch .asc)."""
    run_dir = Path(run_dir)
    files = sorted(run_dir.glob("error.*.xyz.asc"),
                   key=lambda p: int(p.name.split(".")[1]))
    if not files:
        raise SystemExit(f"no Llama error .xyz.asc in {run_dir}")
    # final iteration = max it seen in the cube file (map 0)
    cube_file = run_dir / "error.0.xyz.asc"
    want_it = 0
    with open(cube_file) as fh:
        for line in fh:
            if not line.strip() or line[0] == "#":
                continue
            parts = line.split()
            if len(parts) >= 9:
                want_it = max(want_it, int(parts[0]))
    b = Bands()
    for path in files:
        pidx = int(path.name.split(".")[1])
        if pidx == 0:  # cube
            for x, y, z, e in _asc_rows(str(path), want_it):
                if not owner_is_cube(x, y, z, R):
                    continue  # donor-zone / corner point: wedge owns it
                r = math.sqrt(x * x + y * y + z * z)
                if r >= r_max:
                    continue
                b.add(band_of(r, R, seam), abs(e))
        else:  # wedge
            axis = _wedge_axis(str(path), want_it, R)
            if axis is None:
                continue
            fax, fsign = axis
            for x, y, z, e in _asc_rows(str(path), want_it):
                r = math.sqrt(x * x + y * y + z * z)
                if r <= R or r >= r_max:
                    continue  # inner overlap into cube, or outside region
                a = (abs(x), abs(y), abs(z))
                ax = a.index(max(a))
                sign = 1 if (x, y, z)[ax] >= 0 else -1
                if (ax, sign) != (fax, fsign):
                    continue  # angular-overlap donor: another wedge owns it
                b.add(band_of(r, R, seam), abs(e))
    return b


def report(name, runs, R, r_max, seam, loader):
    print(f"\n{'='*78}\n{name} plane-wave error vs analytic (owner points, r<{r_max})\n{'='*78}")
    series = []
    for d in runs:
        N = dir_N(d)
        b = loader(d, R, r_max, seam)
        series.append((N, b))
        print(f"  N={N:3d}  interior {b.l2('interior'):.4e} (n={b.n['interior']:8d})"
              f"   seam {b.l2('seam'):.4e} (n={b.n['seam']:7d})"
              f"   global {b.l2_global():.4e}")
    series.sort(key=lambda t: t[0])
    print(f"  {'pair':>11s} {'interior':>9s} {'seam':>9s} {'global':>9s}")
    orders = []
    for (Nc, bc), (Nf, bf) in zip(series, series[1:]):
        oi = order(bc.l2("interior"), bf.l2("interior"))
        os = order(bc.l2("seam"), bf.l2("seam"))
        og = order(bc.l2_global(), bf.l2_global())
        orders.append((Nc, Nf, oi, os, og))
        print(f"  {Nc:4d}->{Nf:<4d} {oi:9.3f} {os:9.3f} {og:9.3f}")
    ratios = [(N, b.l2("seam") / b.l2("interior")
               if b.l2("interior") > 0 else float("inf")) for N, b in series]
    print("  seam/interior ratio: " +
          "  ".join(f"N{N}={rt:.3f}" for N, rt in ratios))
    return series, orders, ratios


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--capyrx", nargs="+", required=True,
                    help="CapyrX run dirs, one per resolution (llama_pwave_*)")
    ap.add_argument("--llama", nargs="+", required=True,
                    help="Llama run dirs, one per resolution (lwt_pwave_*)")
    ap.add_argument("--inner-radius", type=float, default=1.0)
    ap.add_argument("--r-max", type=float, default=2.5)
    ap.add_argument("--seam-band", type=float, default=0.25)
    ap.add_argument("--order-floor", type=float, default=3.5,
                    help="minimum acceptable finest-pair order (both codes)")
    ap.add_argument("--seam-order-floor", type=float, default=3.0)
    ap.add_argument("--ratio-bound", type=float, default=3.0,
                    help="max CapyrX seam/interior error ratio")
    ap.add_argument("--ratio-growth", type=float, default=1.5,
                    help="max growth of the CapyrX seam/interior ratio N->2N->4N")
    args = ap.parse_args()

    R = args.inner_radius
    print("#" * 78)
    print("# Step 16 plane-wave A/B: CapyrX Llama vs real Llama (Thornburg04)")
    print("#" * 78)

    cx, cx_ord, cx_ratio = report("CapyrX (this checkout)", args.capyrx, R,
                                  args.r_max, args.seam_band, capyrx_bands)
    ll, ll_ord, ll_ratio = report("real Llama (AxiDil LlamaWaveToy)", args.llama, R,
                                  args.r_max, args.seam_band, llama_bands)

    # run_checks.sh reads this exact label for its row count.
    total_rows = sum(b.count() for _, b in cx) + sum(b.count() for _, b in ll)
    print(f"\nRows checked against ground truth: {total_rows}")

    fail = []
    if len(cx_ord) < 1 or len(ll_ord) < 1:
        print("RESULT: CANNOT-RUN -- need >= 2 resolutions per code")
        sys.exit(2)
    for _, b in cx + ll:
        if b.n["interior"] == 0 or b.n["seam"] == 0:
            print("RESULT: CANNOT-RUN -- an interior/seam band is empty")
            sys.exit(2)

    # finest pair of each code
    cNc, cNf, coi, cos_, cog = cx_ord[-1]
    lNc, lNf, loi, los, log_ = ll_ord[-1]
    if coi < args.order_floor:
        fail.append(f"CapyrX interior order {coi:.3f} < {args.order_floor}")
    if cos_ < args.seam_order_floor:
        fail.append(f"CapyrX seam order {cos_:.3f} < {args.seam_order_floor}")
    if loi < args.order_floor:
        fail.append(f"Llama interior order {loi:.3f} < {args.order_floor}")
    if los < args.seam_order_floor:
        fail.append(f"Llama seam order {los:.3f} < {args.seam_order_floor}")

    # CapyrX-only interface signal: seam/interior ratio bounded and not growing
    rt_fine = cx_ratio[-1][1]
    rt_coarse = cx_ratio[0][1]
    growth = rt_fine / rt_coarse if rt_coarse > 0 else float("inf")
    if rt_fine > args.ratio_bound:
        fail.append(f"CapyrX seam/interior ratio {rt_fine:.3f} > {args.ratio_bound}")
    if growth > args.ratio_growth:
        fail.append(f"CapyrX seam ratio growing (x{growth:.3f} > {args.ratio_growth})")

    print("-" * 78)
    if fail:
        print("RESULT: FAIL -- " + "; ".join(fail))
        sys.exit(1)
    print(f"RESULT: PASS -- CapyrX finest-pair order interior {coi:.3f} / seam "
          f"{cos_:.3f}; Llama {loi:.3f} / {los:.3f} (both >= {args.order_floor}/"
          f"{args.seam_order_floor}); CapyrX seam/interior ratio {rt_fine:.3f} "
          f"bounded, growth x{growth:.3f}. No CapyrX-only interface signal.")
    sys.exit(0)


if __name__ == "__main__":
    main()
