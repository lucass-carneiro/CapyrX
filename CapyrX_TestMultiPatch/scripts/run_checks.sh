#!/usr/bin/env bash
#
# run_checks.sh -- the CapyrX_TestMultiPatch interpatch regression battery.
#
# Runs an ENUMERATED list of par-file cells with the built binary and feeds each
# one to the analytic checker that can judge it. Hand-invokable and
# CI-invokable; no Cactus-testsuite coupling. It does NOT build -- an existing
# binary is required (building is the operator's responsibility).
#
# ---------------------------------------------------------------------------
# THE MATRIX -- 31 runs plus 3 triplet checks = 34 verdicts, enumerated below
# in MATRIX[], never computed as a product. A product is how a matrix comes to
# run 12 of its 300 cells and exit 0; every cell here is written out, and every
# cell that does not reach a verdict is REPORTED, not skipped.
#
#   colour  3 pars x slave_overlap {no,yes} x {1 rank, 2 ranks}   = 12  -> check_color.py
#   nan     10 nan_*.par, 1 rank                                 = 10  -> check_nan_test.py
#   smooth  smooth_{z,p,o}_{neumann,linextrap,none}, 1 rank       =  9  -> check_smooth_test.py
#                                                                          (3 triplet checks)
#
# The 2-rank cells run, and they are not the 1-box-per-patch layout that aborts
# in MPI_Alltoallv. That abort (CarpetX/src/interpolate.cxx, the
# `recvbuf.size() != (nvars+1)*npoints` check after the Alltoallv) needs ONE
# AMReX box per patch, and none of the three colour pars produces that at
# nranks > 1: CarpetX::refine_grid_layout (default yes) re-chops each patch's
# box until there is one per rank, and a par file's own
# CarpetX::amrex_parameters = "amr.refine_grid_layout=0" does NOT override the
# parameter CarpetX sets itself -- it is inert. Therefore:
#
#   * DO NOT add `CarpetX::refine_grid_layout = no` to any multi-rank cell.
#     That one keyword is the reproducer for the upstream bug, so a cell
#     carrying it is a regression test for that bug, not a cell of this
#     battery. Such a cell is REFUSED below rather than run.
#   * A cell's box layout is READ OUT OF ITS LOG ("level L: N boxes", one line
#     per patch per regrid) and printed in the summary. It is never asserted
#     from the par file -- asserting it from the par file is what made the
#     2-rank blocker look fixed for four weeks.
#
# ---------------------------------------------------------------------------
# VERDICTS AND THE RUNNER'S EXIT CODE
#
#   PASS        the cell ran and matched the expectation recorded next to it
#   DRIFT       the cell ran and disagreed with that expectation. NOT
#               automatically a bug: the expectations are recorded
#               measurements, so a drift means "a number this battery pins has
#               moved -- decide whether that is the fix working or a
#               regression". It is never green.
#   REGRESSION  a checker reported incorrect interpatch data where correct data
#               is expected
#   CANNOT-RUN  the run aborted, the TSV is missing, the checker could not run,
#               or the checker validated ZERO rows
#   REFUSED     the cell's configuration is outside this battery (see the
#               refine_grid_layout rule above); nothing was launched
#   DROPPED     the cell was excluded by RANKS= or CELLS=. Named in the summary
#               and counted in the OVERALL line, never silent.
#
#   exit 0  every planned cell ran and passed
#   exit 1  at least one REGRESSION or DRIFT
#   exit 2  no regression, but at least one cell COULD NOT RUN or was REFUSED
#           -- the battery's coverage is not what it claims, so its green cells
#           are not a result
#   1 takes precedence over 2 as the headline because it is the strongest,
#   most actionable signal; both are CI failures. A run narrowed with RANKS= or
#   CELLS= can still exit 0, but the OVERALL line always states how many cells
#   it covered -- CI must run the default matrix.
#
# ---------------------------------------------------------------------------
# WHAT THE CHECKERS CANNOT DO. Stated here so it is not discovered later.
#
#  check_color.py
#   * It ports the production get_owner_patch, so a tie-break bug in that
#     function passes on both sides of the comparison. The circularity is real
#     and is conceded, not papered over.
#   * Its vacuous-pass guard is `checked == 0 => exit 2` plus a threshold on the
#     FRACTION of rows with no coords match. A join that silently drops, say,
#     30 % of its rows still passes. That is why this runner pins the
#     checked-row COUNT per cell (the `rows=` column) instead of trusting the
#     guard.
#   * It has no --slave-overlap mode. With slave_overlap = yes a slaved
#     interior cell legitimately holds another patch's marker, which its
#     own-interior test flags by construction -- so the slave-on colour cells
#     EXPECT exit 1, and their gate is the flagged-cell count, not the exit
#     code. The set-level differential against the pre-fix column lives in
#     evidence/fix/c1/, not here.
#
#  check_nan_test.py
#   * It has NO could-not-run exit code at all: an empty join reads as "no leak"
#     and exits 0. This runner therefore treats `rows checked == 0` as
#     CANNOT-RUN for every kind, and pins the row count per cell.
#   * The leak threshold is interpolation order 1 -- not 0, and not 2. Its old
#     corner-priority bug, which made order 2 look like the threshold, is
#     fixed.
#
#  check_smooth_test.py
#   * exit 1 = "Channel 1 is LIVE" is the EXPECTED POSITIVE, not a failure.
#     Only exit 2 (strict-interior invariant violated, i.e. the A/B premise is
#     broken, or the checker could not run) is a real failure.
#   * Its interior-invariant guard SKIPS ITSELF, printing one line, if it cannot
#     recover angular_cells/radial_cells from the par file -- and the exit code
#     does not change. This runner reads that line and calls a SKIPPED guard
#     CANNOT-RUN.
#   * The one_over_r field is singular at the origin, so the `o` triplet carries
#     a small number of non-finite cells the checker skips rather than flags.
#
# ---------------------------------------------------------------------------
# STREAMS: the two are kept SEPARATE (<id>.run.out / <id>.run.err) and must
# stay that way. A CCTK warning can land on BOTH: the flesh sends it to stderr
# when its level <= warning_level (default 1) and to stdout when its level <=
# logging_level -- which the -L3 below raises to 3 -- or, on any non-root rank,
# whenever level <= warning_level (Cactus/src/main/WarnLevel.c). So under -L3
# levels 0 and 1 appear twice while 2 and 3 are stdout-only, and the two copies
# are not even byte-identical (mpiexec pty-forwards rank 0's stdout, so Cactus
# bolds the stdout copy and not the stderr one). A merged log can therefore be
# neither stream-attributed nor de-duplicated, and its duplicated lines read as
# twice the evidence. (The Python checkers below keep 2>&1 on purpose: they are
# single-stream, and only their exit code and their stdout RESULT line are ever
# consumed.)
#
# LAUNCHER: derived from the binary's own resolved libmpi
# (<libmpi dir>/../bin/mpiexec), not from $PATH. On the machine this battery was
# written on, $PATH's mpiexec is Open MPI 4.1.6 while the binary links 5.0.5;
# that mismatch is harmless-looking at one rank and is the first thing a 2-rank
# cell hits. Set MPIEXEC=... to override, MPIEXEC=none to exec directly (1 rank
# only).
#
# ---------------------------------------------------------------------------
# Environment overrides:
#   CACTUS_SIM   path to the cactus binary   (default: <repo>/../../../exe/cactus_sim)
#   EXE_DIR      dir to run in / write output (default: dir of CACTUS_SIM)
#   RANKS        rank axis for the colour cells   (default: "1 2")
#   NRANKS       deprecated alias: NRANKS=n means RANKS="n"
#   CELLS        space-separated cell ids to run  (default: all)
#   MPIEXEC      launcher, or "none" for a direct exec
#   MPI_ARGS     launcher placement flags   (default: --map-by socket --bind-to socket)
#   SKIP_RUN=1   skip the cactus runs and check whatever output already exists
#                (use to re-check, or to point perturbed TSVs at the checkers).
#                A cell with no output is CANNOT-RUN, not skipped -- narrow the
#                matrix with CELLS= if that is what you mean.
#
set -uo pipefail
shopt -s nullglob

# =============================== THE RUN LIST ================================
# MATRIX rows: id | kind | par (basename, no .par) | slave | nranks | expect
# (the array is MATRIX, not CELLS: CELLS is the environment override that
#  selects a subset of it, and an array of that name would shadow it.)
#
# expect := <checker exit code>[:<KIND>=<count>[,<KIND>=<count>...]][:rows=<n>]
#   KIND is matched against the checker's own "=== KIND (N) ===" headings.
#   Every count here is a RECORDED MEASUREMENT with a source, not a derivation:
#     colour rows= / OWN-*-MISMATCH   : measured 2026-08-25 at CarpetX 4e50b591,
#                                       CapyrX 58e2060, evidence/fix/c1/. The
#                                       1-rank slave-on numbers are unchanged
#                                       from B10 and the 2086 / 2728 slaved sets
#                                       are the ones [P6] measured
#     nan LEAK-CROSS-PATCH            : mp_noise.md Part C, re-measured at every
#                                       Phase-B commit
#     smooth cross-patch cell count   : mp_noise.md's per-field table (Part B) --
#                                       see the note above SMOOTH_TRIPLETS
# A count that moves is a DRIFT, and the runner says so rather than deciding for
# you which direction is good.
MATRIX=(
  # ---- colour, 1 rank, slave off: the correctness gate ----------------------
  "color|color|color|no|1|0:rows=23491"
  "color_ghost|color|color_ghost|no|1|0:rows=75816"
  "color_ghost_overlap|color|color_ghost_overlap|no|1|0:rows=161875"
  # ---- colour, 1 rank, slave on: exit 1 BY CONSTRUCTION (see the header) ----
  "color_slave|color|color|yes|1|1:OWN-INTERIOR-MISMATCH=2086:rows=23491"
  "color_ghost_slave|color|color_ghost|yes|1|1:OWN-INTERIOR-MISMATCH=2728:rows=75816"
  "color_ghost_overlap_slave|color|color_ghost_overlap|yes|1|1:OWN-INTERIOR-MISMATCH=8006,OWN-OVERLAP-MISMATCH=39736:rows=161875"
  # ---- colour, 2 ranks: chopped layout, so the Alltoallv abort cannot fire --
  # These counts are NOT the 1-rank ones and must not be copied from them: the
  # chopped layout gives each patch more boxes, a cell appears once per box that
  # holds it, and that multiplicity is part of the count. color_ghost is the
  # exception -- it is already 16 boxes on six patches and 8 on the seventh at
  # ONE rank, so two ranks add no boxes and its counts do not move.
  # Measured 2026-08-25, evidence/fix/c1/.
  "color_n2|color|color|no|2|0:rows=29406"
  "color_ghost_n2|color|color_ghost|no|2|0:rows=75816"
  "color_ghost_overlap_n2|color|color_ghost_overlap|no|2|0:rows=183750"
  "color_slave_n2|color|color|yes|2|1:OWN-INTERIOR-MISMATCH=2186:rows=29406"
  "color_ghost_slave_n2|color|color_ghost|yes|2|1:OWN-INTERIOR-MISMATCH=2728:rows=75816"
  "color_ghost_overlap_slave_n2|color|color_ghost_overlap|yes|2|1:OWN-INTERIOR-MISMATCH=8202,OWN-OVERLAP-MISMATCH=40800:rows=183750"
  # ---- NaN injection: the leak threshold is order 1 ------------------------
  "nan_dirichlet_ghost|nan|nan_dirichlet_ghost|no|1|0:rows=23491"
  "nan_neumann_ghost_and_interior|nan|nan_neumann_ghost_and_interior|no|1|0:rows=23491"
  "nan_neumann_interior|nan|nan_neumann_interior|no|1|0:rows=23491"
  "nan_neumann_ghost_overlap|nan|nan_neumann_ghost_overlap|no|1|0:rows=161875"
  "nan_neumann_ghost_overlap_order4|nan|nan_neumann_ghost_overlap_order4|no|1|0:rows=161875"
  "nan_neumann_interior_overlap_order0|nan|nan_neumann_interior_overlap_order0|no|1|0:LEAK-CROSS-PATCH=0:rows=161875"
  "nan_neumann_interior_overlap_order1|nan|nan_neumann_interior_overlap_order1|no|1|1:LEAK-CROSS-PATCH=1104:rows=161875"
  "nan_neumann_interior_overlap_order2|nan|nan_neumann_interior_overlap_order2|no|1|1:LEAK-CROSS-PATCH=1512:rows=161875"
  "nan_neumann_interior_overlap_order3|nan|nan_neumann_interior_overlap_order3|no|1|1:LEAK-CROSS-PATCH=2208:rows=161875"
  "nan_neumann_interior_overlap_order4|nan|nan_neumann_interior_overlap_order4|no|1|1:LEAK-CROSS-PATCH=3312:rows=161875"
  # ---- smooth-field Channel-1 A/B, three fields x three outer BCs ----------
  # Checked as triplets (see SMOOTH_TRIPLETS); the per-leg expect is unused.
  "smooth_z_neumann|smoothleg|smooth_z_neumann|no|1|-"
  "smooth_z_linextrap|smoothleg|smooth_z_linextrap|no|1|-"
  "smooth_z_none|smoothleg|smooth_z_none|no|1|-"
  "smooth_p_neumann|smoothleg|smooth_p_neumann|no|1|-"
  "smooth_p_linextrap|smoothleg|smooth_p_linextrap|no|1|-"
  "smooth_p_none|smoothleg|smooth_p_none|no|1|-"
  "smooth_o_neumann|smoothleg|smooth_o_neumann|no|1|-"
  "smooth_o_linextrap|smoothleg|smooth_o_linextrap|no|1|-"
  "smooth_o_none|smoothleg|smooth_o_none|no|1|-"
)
# field | neumann-cell | linextrap-cell | none-cell | expected exit | expected
# cross-patch cell count. Exit 0 (inert) and 1 (LIVE) are both PASS for the
# CHECKER; the pinned pair below is the recorded state, and a departure from it
# is a DRIFT.
# The cross-patch counts are per FIELD and are not all the same number: 1512
# for parabola and one_over_r, 1488 for z_global, which undercounts by 24 cells
# (~1.6 %) in an equatorial blind spot where the field's own gradient vanishes.
# They come from mp_noise.md's per-field table, not from its prose summary --
# "1488-1512" and "1488/1512" appear in that file's headline sentences and in a
# memory index, and pinning 1488 for all three off those sentences is what this
# battery's first run caught.
SMOOTH_TRIPLETS=(
  "smooth_z|smooth_z_neumann|smooth_z_linextrap|smooth_z_none|1|1488"
  "smooth_p|smooth_p_neumann|smooth_p_linextrap|smooth_p_none|1|1512"
  "smooth_o|smooth_o_neumann|smooth_o_linextrap|smooth_o_none|1|1512"
)
# TSV group names emitted by each kind of run. The glob covers the per-rank
# files a multi-rank run writes (.p0000, .p0001, ...); the checkers take several
# TSVs and find each one's own coords/pre companion.
COLOR_TSV_GLOB="capyrx_testmultipatch-color.it000000.p*.tsv"
NAN_TSV_GLOB="capyrx_testmultipatch-nan_test.it000000.p*.tsv"
SMOOTH_TSV_GLOB="capyrx_testmultipatch-smooth_test.it000000.p*.tsv"
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THORN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PAR_DIR="$THORN_DIR/par"
CHECK_COLOR="$SCRIPT_DIR/check_color.py"
CHECK_NAN="$SCRIPT_DIR/check_nan_test.py"
CHECK_SMOOTH="$SCRIPT_DIR/check_smooth_test.py"

# scripts/../../../../exe  ==  <cactus>/exe
CACTUS_SIM="${CACTUS_SIM:-$(cd "$SCRIPT_DIR/../../../.." 2>/dev/null && pwd)/exe/cactus_sim}"
EXE_DIR="${EXE_DIR:-$(dirname "$CACTUS_SIM")}"
SKIP_RUN="${SKIP_RUN:-0}"
RANKS="${RANKS:-${NRANKS:-1 2}}"
CELLS_FILTER="${CELLS:-}"
MPI_ARGS="${MPI_ARGS:---map-by socket --bind-to socket}"

: "${OMP_NUM_THREADS:=8}"
: "${OMP_PLACES:=cores}"
: "${OMP_PROC_BIND:=close}"
export OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND

LOG_DIR="$EXE_DIR/run_checks_logs"
GEN_PAR_DIR="$EXE_DIR/run_checks_pars"

fail_early() { echo "run_checks.sh: $*" >&2; exit 2; }

[[ -x "$CACTUS_SIM" ]] || fail_early "cactus binary not found or not executable: $CACTUS_SIM (build it first, or set CACTUS_SIM)"
[[ -d "$EXE_DIR"    ]] || fail_early "output/run dir does not exist: $EXE_DIR (set EXE_DIR)"
for c in "$CHECK_COLOR" "$CHECK_NAN" "$CHECK_SMOOTH"; do
  [[ -f "$c" ]] || fail_early "missing checker: $c"
done
mkdir -p "$LOG_DIR" "$GEN_PAR_DIR"

# ---- launcher ---------------------------------------------------------------
resolve_launcher() {
  if [[ "${MPIEXEC:-}" == none ]]; then LAUNCHER=none; return 0; fi
  if [[ -n "${MPIEXEC:-}" ]]; then
    command -v "$MPIEXEC" >/dev/null || fail_early "MPIEXEC='$MPIEXEC' not found"
    LAUNCHER="$MPIEXEC"; return 0
  fi
  local libmpi bindir
  libmpi=$(ldd "$CACTUS_SIM" 2>/dev/null | sed -n 's/.*libmpi\.so[^ ]* => \([^ ]*\).*/\1/p' | head -1)
  if [[ -n "$libmpi" ]]; then
    bindir=$(cd "$(dirname "$libmpi")/../bin" 2>/dev/null && pwd)
    if [[ -n "$bindir" && -x "$bindir/mpiexec" ]]; then LAUNCHER="$bindir/mpiexec"; return 0; fi
  fi
  LAUNCHER=$(command -v mpiexec) || fail_early "no mpiexec: cannot derive one from $CACTUS_SIM and none on \$PATH (set MPIEXEC)"
  echo "run_checks.sh: WARNING -- could not derive a launcher from the binary's libmpi;" >&2
  echo "               falling back to \$PATH's $LAUNCHER ($("$LAUNCHER" --version 2>&1 | head -1))." >&2
  echo "               A version mismatch with the binary's MPI shows up first in the multi-rank cells." >&2
}
LAUNCHER=""
[[ "$SKIP_RUN" == "1" ]] || resolve_launcher

# ---- result accumulation ----------------------------------------------------
declare -a SUMMARY=()
had_regression=0; had_cantrun=0
n_planned=0; n_ran=0; n_dropped=0; n_refused=0
declare -a DROPPED_IDS=()

record() {  # record <status> <id> <ranks> <boxes> <rc> <rows> <detail>
  local status="$1"
  case "$status" in
    REGRESSION|DRIFT)      had_regression=1 ;;
    CANNOT-RUN|REFUSED)    had_cantrun=1 ;;
  esac
  SUMMARY+=("$(printf '%-11s %-30s %-2s %-9s %-4s %-9s %s' \
    "$status" "$2" "$3" "$4" "$5" "$6" "$7")")
}

# ---- readers ----------------------------------------------------------------
# The box layout, read from the log and never from the par file. CarpetX prints
# one "level L: N boxes" line per patch per regrid.
boxes_of() {  # boxes_of <run.out> -> "l0:2x7" | "-"
  local f="$1"
  [[ -f "$f" ]] || { printf -- '-'; return; }
  sed 's/\x1b\[[0-9;]*m//g' "$f" \
    | sed -n 's/.*level \([0-9]\+\): \([0-9]\+\) boxes.*/\1 \2/p' \
    | sort | uniq -c \
    | awk '{printf "%sl%s:%sx%s", (NR>1?",":""), $2, $3, $1} END{if(NR==0) printf "-"}'
}
rows_of() {  # rows_of <check.log> -> checked-row count, or -1 if unreadable
  local n
  n=$(sed -n 's/^Rows checked against ground truth:[[:space:]]*\([0-9]\+\).*/\1/p;
              s/^cells compared across neumann\/linextrap:[[:space:]]*\([0-9]\+\).*/\1/p' \
        "$1" 2>/dev/null | head -1)
  case "$n" in ''|*[!0-9]*) printf -- '-1' ;; *) printf '%s' "$n" ;; esac
}
count_kind() {  # count_kind <check.log> <KIND> -> N (0 if the heading is absent)
  awk -v k="$2" '
    index($0, "=== " k) == 1 {
      line = $0
      sub(/ ===[[:space:]]*$/, "", line)
      if (match(line, /\([0-9]+\)$/)) { print substr(line, RSTART+1, RLENGTH-2); found=1; exit }
    }
    END { if (!found) print 0 }' "$1" 2>/dev/null
}
smooth_cross_of() {  # smooth_cross_of <check.log> -> cross-patch cell count
  local n
  n=$(sed -n 's/^RESULT: Channel 1 is LIVE\. The outer-BC choice changes \([0-9]\+\) cross-patch.*/\1/p' \
        "$1" 2>/dev/null | head -1)
  if [[ -z "$n" ]] && grep -q '^RESULT: Channel 1 is INERT' "$1" 2>/dev/null; then n=0; fi
  case "$n" in ''|*[!0-9]*) printf -- '-1' ;; *) printf '%s' "$n" ;; esac
}

# ---- derived par files ------------------------------------------------------
# A cell that is not (slave off, 1 rank) needs its own par file, because
# IO::out_dir = $parfile means the output directory IS the par-file basename --
# two cells sharing a par file would share, and overwrite, one output tree. The
# derived file is a copy plus AT MOST ONE functional line, regenerated on every
# invocation so it cannot drift from the original.
derive_par() {  # derive_par <id> <base-par> <slave> -> echoes the par path
  local id="$1" base="$2" slave="$3" src="$PAR_DIR/$2.par" dst="$GEN_PAR_DIR/$1.par"
  if [[ "$id" == "$base" ]]; then printf '%s' "$src"; return 0; fi
  { cat "$src"
    printf '\n# --- derived by run_checks.sh, do not edit -----------------------------\n'
    printf '# Cell "%s" of the battery: %s.par with slave_overlap = %s.\n' "$id" "$base" "$slave"
    printf '# It exists as a separate file only because IO::out_dir = $parfile, so the\n'
    printf '# output directory is the par-file name and two cells must not share one.\n'
    if [[ "$slave" == yes ]]; then
      printf 'CapyrX_MultiPatch::slave_overlap = yes\n'
    fi
  } > "$dst" || return 1
  printf '%s' "$dst"
}

# A cell carrying CarpetX::refine_grid_layout = no at more than one rank is the
# upstream Alltoallv bug's reproducer, not a cell of this battery. (The
# amrex_parameters form of the same key is inert and is deliberately NOT matched
# here.)
forces_one_box() {  # forces_one_box <par-path> -> 0 if it disables refine_grid_layout
  grep -Eiq '^[[:space:]]*CarpetX::refine_grid_layout[[:space:]]*=[[:space:]]*"?(no|false|0)"?' "$1"
}

# ---- the runner -------------------------------------------------------------
run_cell() {  # run_cell <id> <par-path> <nranks> -> cactus rc (127 = never launched)
  local id="$1" par="$2" nr="$3"
  local out="$LOG_DIR/$id.run.out" err="$LOG_DIR/$id.run.err"
  if [[ "$SKIP_RUN" == "1" ]]; then
    # Do NOT touch <id>.run.{out,err}. They are the archived run's own streams
    # and the ONLY record of that run's box layout, which the summary reads back
    # out of them -- writing a "reusing existing" note over them destroys the
    # evidence the re-check exists to re-examine, and leaves the layout column
    # blank for every cell.
    [[ -f "$out" ]] || echo "run_checks.sh: SKIP_RUN=1 and no archived log for cell $id" >&2
    return 0
  fi
  rm -rf "${EXE_DIR:?}/$id"
  if [[ "$LAUNCHER" == none ]]; then
    [[ "$nr" == 1 ]] || { echo "MPIEXEC=none cannot run $nr ranks" >"$err"; : >"$out"; return 127; }
    ( cd "$EXE_DIR" && "$CACTUS_SIM" -L3 "$par" ) >"$out" 2>"$err"
  else
    ( cd "$EXE_DIR" && "$LAUNCHER" $MPI_ARGS -n "$nr" "$CACTUS_SIM" -L3 "$par" ) >"$out" 2>"$err"
  fi
}

# expect_check <check.log> <expect> <rc> -> "" if it matches, else the drift text
expect_check() {
  local log="$1" expect="$2" rc="$3" detail="" body want_rc field k want got
  want_rc=${expect%%:*}
  [[ "$rc" == "$want_rc" ]] || detail="exit $rc != recorded $want_rc"
  body=${expect#"$want_rc"}
  while [[ -n "$body" ]]; do
    body=${body#:}
    field=${body%%:*}
    body=${body#"$field"}
    local IFS=','
    for kv in $field; do
      k=${kv%%=*}; want=${kv#*=}
      if [[ "$k" == rows ]]; then got=$(rows_of "$log"); else got=$(count_kind "$log" "$k"); fi
      [[ "$got" == "$want" ]] || detail="$detail${detail:+; }$k=$got != recorded $want"
    done
  done
  printf '%s' "$detail"
}

echo "== run_checks.sh =="
echo "binary   : $CACTUS_SIM"
echo "exe dir  : $EXE_DIR"
echo "launcher : ${LAUNCHER:-<not resolved: SKIP_RUN=1>}${LAUNCHER:+  ($MPI_ARGS)}"
echo "ranks    : $RANKS      omp: $OMP_NUM_THREADS      SKIP_RUN=$SKIP_RUN"
echo "cells    : ${CELLS_FILTER:-all ${#MATRIX[@]}}"
echo "logs     : $LOG_DIR"
echo

declare -A SMOOTH_OK=()

for cell in "${MATRIX[@]}"; do
  IFS='|' read -r id kind base slave nr expect <<<"$cell"
  n_planned=$((n_planned + 1))

  # -- the two ways a cell is deliberately excluded, both of them loud --------
  if [[ -n "$CELLS_FILTER" && " $CELLS_FILTER " != *" $id "* ]]; then
    n_dropped=$((n_dropped + 1)); DROPPED_IDS+=("$id")
    record "DROPPED" "$id" "$nr" "-" "-" "-" "not in CELLS=$CELLS_FILTER"; continue
  fi
  if [[ " $RANKS " != *" $nr "* ]]; then
    n_dropped=$((n_dropped + 1)); DROPPED_IDS+=("$id")
    record "DROPPED" "$id" "$nr" "-" "-" "-" "rank $nr not in RANKS=$RANKS"; continue
  fi

  if [[ ! -f "$PAR_DIR/$base.par" ]]; then
    record "CANNOT-RUN" "$id" "$nr" "-" "-" "-" "missing par file: $PAR_DIR/$base.par"; continue
  fi
  par=$(derive_par "$id" "$base" "$slave") || {
    record "CANNOT-RUN" "$id" "$nr" "-" "-" "-" "could not write the derived par under $GEN_PAR_DIR"; continue; }

  if [[ "$nr" != 1 ]] && forces_one_box "$par"; then
    n_refused=$((n_refused + 1))
    record "REFUSED" "$id" "$nr" "-" "-" "-" \
      "CarpetX::refine_grid_layout is disabled: one box per patch at $nr ranks is the upstream Alltoallv reproducer, not a cell of this battery"
    continue
  fi

  echo "-- $id  ($base.par, slave=$slave, $nr rank(s))"
  run_cell "$id" "$par" "$nr"; run_rc=$?
  boxes=$(boxes_of "$LOG_DIR/$id.run.out")
  if [[ "$run_rc" != 0 ]]; then
    record "CANNOT-RUN" "$id" "$nr" "$boxes" "$run_rc" "-" \
      "cactus exited $run_rc (see $LOG_DIR/$id.run.{out,err})"; continue
  fi
  n_ran=$((n_ran + 1))

  case "$kind" in
    color)  tsvs=( "$EXE_DIR/$id"/$COLOR_TSV_GLOB ) ;;
    nan)    tsvs=( "$EXE_DIR/$id"/$NAN_TSV_GLOB ) ;;
    smoothleg) tsvs=( "$EXE_DIR/$id"/$SMOOTH_TSV_GLOB ) ;;
  esac
  if [[ ${#tsvs[@]} -eq 0 ]]; then
    record "CANNOT-RUN" "$id" "$nr" "$boxes" 0 "-" "the run produced no $kind TSV in $EXE_DIR/$id"
    continue
  fi

  if [[ "$kind" == smoothleg ]]; then
    SMOOTH_OK[$id]=1
    record "PASS" "$id" "$nr" "$boxes" 0 "-" "ran; judged as part of its triplet"
    continue
  fi

  clog="$LOG_DIR/$id.check.log"
  case "$kind" in
    color) python3 "$CHECK_COLOR" --max-examples 20 "${tsvs[@]}" >"$clog" 2>&1 ;;
    nan)   python3 "$CHECK_NAN"   --max-examples 20 "${tsvs[@]}" >"$clog" 2>&1 ;;
  esac
  rc=$?
  rows=$(rows_of "$clog")
  result_line="$(grep -m1 '^RESULT:' "$clog" | cut -c1-120 || true)"

  if [[ "$rc" -gt 2 ]]; then
    record "CANNOT-RUN" "$id" "$nr" "$boxes" "$rc" "$rows" "checker crashed (see $clog)"; continue
  fi
  if [[ "$rc" == 2 ]]; then
    record "CANNOT-RUN" "$id" "$nr" "$boxes" "$rc" "$rows" "${result_line:-checker could not run; see $clog}"; continue
  fi
  # check_nan_test.py has no could-not-run exit code at all, so an empty join
  # would read as a clean pass. Close that here, for every kind.
  if [[ "$rows" -le 0 ]]; then
    record "CANNOT-RUN" "$id" "$nr" "$boxes" "$rc" "$rows" \
      "the checker validated zero rows -- broken join or truncated output (see $clog)"; continue
  fi
  drift=$(expect_check "$clog" "$expect" "$rc")
  if [[ -n "$drift" ]]; then
    record "DRIFT" "$id" "$nr" "$boxes" "$rc" "$rows" "$drift -- ${result_line:-}"
  elif [[ "$rc" == 1 && "${expect%%:*}" != 1 ]]; then
    record "REGRESSION" "$id" "$nr" "$boxes" "$rc" "$rows" "${result_line:-incorrect interpatch data}"
  else
    record "PASS" "$id" "$nr" "$boxes" "$rc" "$rows" "${result_line:-}"
  fi
done

# ---- the smooth triplets ----------------------------------------------------
for t in "${SMOOTH_TRIPLETS[@]}"; do
  IFS='|' read -r field neu_id lin_id non_id want_rc want_cross <<<"$t"
  n_planned=$((n_planned + 1))
  missing=""
  for leg in "$neu_id" "$lin_id" "$non_id"; do
    [[ -n "${SMOOTH_OK[$leg]:-}" ]] || missing="$missing $leg"
  done
  if [[ -n "$missing" ]]; then
    # If every leg was dropped on purpose, the triplet is dropped too; if only
    # some are missing, the triplet could not run.
    if [[ -n "$CELLS_FILTER" && " $CELLS_FILTER " != *" $neu_id "* ]]; then
      n_dropped=$((n_dropped + 1)); DROPPED_IDS+=("${field}_triplet")
      record "DROPPED" "${field}_triplet" 1 "-" "-" "-" "legs not in CELLS="
    else
      record "CANNOT-RUN" "${field}_triplet" 1 "-" "-" "-" "leg(s) never produced output:$missing"
    fi
    continue
  fi
  echo "-- ${field}_triplet  (neumann / linextrap / none)"
  # --neumann/--linextrap/--none each take ONE path, and the smooth cells are
  # 1-rank only, so each glob below resolves to exactly one file. A multi-rank
  # smooth cell would need the checker to accept several TSVs per leg first.
  clog="$LOG_DIR/${field}_triplet.check.log"
  python3 "$CHECK_SMOOTH" --max-examples 20 \
    --neumann   "$EXE_DIR/$neu_id"/$SMOOTH_TSV_GLOB \
    --linextrap "$EXE_DIR/$lin_id"/$SMOOTH_TSV_GLOB \
    --none      "$EXE_DIR/$non_id"/$SMOOTH_TSV_GLOB >"$clog" 2>&1
  rc=$?
  rows=$(rows_of "$clog")
  cross=$(smooth_cross_of "$clog")
  result_line="$(grep -m1 '^RESULT:' "$clog" | cut -c1-120 || true)"
  n_ran=$((n_ran + 1))
  if [[ "$rc" -gt 2 || "$rc" == 2 ]]; then
    record "CANNOT-RUN" "${field}_triplet" 1 "-" "$rc" "$rows" \
      "${result_line:-interior invariant broken, or the checker could not run; see $clog}"
    continue
  fi
  # The guard that silently disables itself. A skipped guard means the A/B
  # premise was never verified, which is exactly the "could not run" case.
  if grep -q '^interior-invariant guard: SKIPPED' "$clog"; then
    record "CANNOT-RUN" "${field}_triplet" 1 "-" "$rc" "$rows" \
      "the checker's interior-invariant guard SKIPPED itself (cell counts unrecoverable); the A/B premise is unverified"
    continue
  fi
  if [[ "$rows" -le 0 ]]; then
    record "CANNOT-RUN" "${field}_triplet" 1 "-" "$rc" "$rows" "zero cells compared (see $clog)"
    continue
  fi
  detail=""
  [[ "$rc" == "$want_rc" ]] || detail="exit $rc != recorded $want_rc"
  [[ "$cross" == "$want_cross" ]] || detail="$detail${detail:+; }cross-patch=$cross != recorded $want_cross"
  if [[ -n "$detail" ]]; then
    record "DRIFT" "${field}_triplet" 1 "-" "$rc" "$rows" "$detail -- ${result_line:-}"
  else
    record "PASS" "${field}_triplet" 1 "-" "$rc" "$rows" "${result_line:-}"
  fi
done

# ---- summary ----------------------------------------------------------------
echo
echo "== summary =="
printf '  %-11s %-30s %-2s %-9s %-4s %-9s %s\n' STATUS CELL NR BOXES RC ROWS DETAIL
for line in "${SUMMARY[@]}"; do echo "  $line"; done
echo
echo "  BOXES is read from each run's own log ('level L: N boxes'), as lN:<boxes>x<patches>."
echo "  ROWS is the checker's own 'rows checked against ground truth'."
echo
printf 'COVERAGE: %d cells planned, %d ran, %d dropped, %d refused.\n' \
  "$n_planned" "$n_ran" "$n_dropped" "$n_refused"
[[ "$n_dropped" == 0 ]] || printf 'DROPPED : %s\n' "${DROPPED_IDS[*]}"

if [[ "$n_ran" == 0 ]]; then
  echo "OVERALL: FAIL -- nothing ran. A battery that covered no cell is not a result."
  exit 2
elif [[ "$had_regression" == "1" ]]; then
  echo "OVERALL: FAIL -- a cell regressed, or a pinned number drifted."
  exit 1
elif [[ "$had_cantrun" == "1" ]]; then
  echo "OVERALL: FAIL -- a cell could not run or was refused; the coverage above is not what it claims."
  exit 2
else
  echo "OVERALL: PASS -- every cell that ran matched its recorded state."
  exit 0
fi
