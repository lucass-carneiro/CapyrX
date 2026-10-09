#!/usr/bin/env python3
"""
Llama cube-only AMR structural checker -- Step 13 of llama_patch_system_impl.md
(Phase D: mesh refinement confined to the central cube).

The parfile llama_amr_32.par adds one refinement level to a Llama WaveToy run with
a BoxInBox region placed inside the cube. This checker reads the `level` column of
the state TSV and asserts the two structural guarantees Step 13 is about:

  * the CUBE (patch 0) carries a refined (level >= 1) band -- AMR actually
    happened inside the cube, and
  * every WEDGE (patches 1-6) holds level 0 ONLY -- the six Thornburg wedges
    stayed unigrid, i.e. the box-in-box tags were confined to [-R,R]^3.

"Ghosts re-filled after each regrid" (R4) is gated by the run itself, not here:
the parfile sets poison_undefined_values = yes, so any interpatch ghost the
regrid repair failed to re-fill would have aborted the run before this checker
ever sees a TSV. A clean completion (checker reached) is that gate.

The TSV columns are, per CarpetX io_tsv.cxx:
  0:iteration 1:time 2:patch 3:level 4:component 5:i 6:j 7:k ...
The cube is patch 0 in the Llama make_system (llama.cxx), matching BoxInBox_Setup,
which refines patch 0 only (boxinbox.cxx:109-115).

Exit codes (the run_checks.sh contract): 0 = pass, 1 = structural failure,
2 = cannot-run (no usable rows).
"""

import argparse
import sys

PATCH_COL = 2
LEVEL_COL = 3
CUBE_PATCH = 0


def load_rows(paths):
    """Return (max_iteration, [(patch, level), ...]) over the final iteration."""
    by_iter = {}
    for path in paths:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) <= LEVEL_COL:
                    continue
                try:
                    it = int(fields[0])
                    patch = int(fields[PATCH_COL])
                    level = int(fields[LEVEL_COL])
                except ValueError:
                    continue
                by_iter.setdefault(it, []).append((patch, level))
    if not by_iter:
        return None, []
    last = max(by_iter)
    return last, by_iter[last]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("tsv", nargs="+", help="state TSV file(s) with a level column")
    ap.add_argument("--max-examples", type=int, default=20,
                    help="cap on listed example cells per heading")
    args = ap.parse_args()

    last_it, rows = load_rows(args.tsv)
    if not rows:
        print("Rows checked against ground truth: 0")
        print("RESULT: CANNOT-RUN -- no data rows found in the state TSV(s); "
              "the run produced no level column to check.")
        return 2

    print(f"iteration checked: {last_it}")
    print(f"Rows checked against ground truth: {len(rows)}")

    cube_refined = [(p, l) for (p, l) in rows if p == CUBE_PATCH and l >= 1]
    wedge_refined = [(p, l) for (p, l) in rows if p != CUBE_PATCH and l >= 1]
    levels_seen = sorted({l for (_, l) in rows})
    patches_seen = sorted({p for (p, _) in rows})

    print(f"patches present: {patches_seen}")
    print(f"levels present:  {levels_seen}")
    print()

    print(f"=== CUBE-REFINED ({len(cube_refined)}) ===")
    print(f"  cube (patch {CUBE_PATCH}) cells at level >= 1: {len(cube_refined)}")

    print(f"=== WEDGE-OVER-LEVEL0 ({len(wedge_refined)}) ===")
    for (p, l) in wedge_refined[:args.max_examples]:
        print(f"  patch {p} has a cell at level {l} -- wedge is NOT unigrid")
    if len(wedge_refined) > args.max_examples:
        print(f"  ... and {len(wedge_refined) - args.max_examples} more")
    print()

    if wedge_refined:
        print(f"RESULT: {len(wedge_refined)} wedge cell(s) are refined above level 0 -- "
              "refinement leaked out of the cube; confinement is broken.")
        return 1
    if not cube_refined:
        print("RESULT: the cube carries no level-1 cells -- AMR did not refine the "
              "cube at all; the test is vacuous.")
        return 1

    print("RESULT: cube-only AMR confirmed -- the cube holds a refined level-1 band "
          "and all six wedges are unigrid (level 0).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
