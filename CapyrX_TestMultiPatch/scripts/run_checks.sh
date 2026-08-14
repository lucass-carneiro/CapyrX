#!/usr/bin/env bash
#
# run_checks.sh -- standalone interpatch-check runner (mp_noise sec-12.9 CI).
#
# Runs the representative CapyrX_TestMultiPatch par files with the built binary
# and feeds their output to the two analytic interpatch checkers. Hand-invokable
# and CI-invokable; no Cactus-testsuite coupling. It does NOT build -- an
# existing binary is required (building is the operator's responsibility, see
# instructions.md).
#
# Coverage:
#   * color.par                          -> check_color.py
#       Correctness gate. color.par leaves run_outer_bc_tests=no (its default),
#       so the NaN-injecting NanTest/SmoothTest groups are NOT scheduled and the
#       poison_undefined_values=yes run completes. check_color exits
#       0=clean, 1=incorrect ghost data, 2=could-not-run (broken/empty join).
#   * smooth_z_{neumann,linextrap,none}  -> check_smooth_test.py
#       Channel-1 A/B *measurement*, not a pass/fail gate. On a healthy build it
#       reports "Channel 1 is LIVE" and exits 1 -- that is the expected positive,
#       NOT a failure. Only exit 2 (strict-interior invariant violated, i.e. the
#       A/B premise is broken, or the checker could not run) is a real failure.
#       So here: exit 0 (inert) OR 1 (live) => PASS; exit 2 => FAIL.
#
# Runner exit code:
#   0  every check passed
#   1  at least one REGRESSION was found (check_color reported incorrect data)
#   2  no regression, but at least one check COULD NOT RUN (aborted run, broken
#      join, or smooth interior-invariant violated) -- results are untrustworthy
# Regression (1) takes precedence over could-not-run (2) as the headline code
# because it is the strongest, most actionable signal; both are CI failures.
#
# Environment overrides:
#   CACTUS_SIM   path to the cactus binary   (default: <repo>/../../../exe/cactus_sim)
#   EXE_DIR      dir to run in / write output (default: dir of CACTUS_SIM)
#   NRANKS       MPI ranks                    (default: 1 -- see note below)
#   SKIP_RUN=1   skip the cactus runs and check whatever TSVs already exist
#                (use to re-check, or to point a perturbed TSV at the checkers)
#
# NOTE: NRANKS is 1 on purpose. color_ghost_overlap-style 1-box-per-patch runs
# with -n 2 abort in MPI_Alltoallv at CarpetX interpolate.cxx:917 (mp_noise
# sec-11); that is a separate, unaddressed blocker.

set -uo pipefail

# ---- Coverage configuration (edit these arrays to extend) --------------------
# Par-file basenames (without .par); output lands in "$EXE_DIR/<basename>/" because
# every par sets IO::out_dir = $parfile.
COLOR_PARS=( "color" )
SMOOTH_TRIPLET=( "smooth_z_neumann" "smooth_z_linextrap" "smooth_z_none" )
# TSV group names emitted by each kind of run:
COLOR_TSV="capyrx_testmultipatch-color.it000000.p0000.tsv"
SMOOTH_TSV="capyrx_testmultipatch-smooth_test.it000000.p0000.tsv"
# -----------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THORN_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PAR_DIR="$THORN_DIR/par"
CHECK_COLOR="$SCRIPT_DIR/check_color.py"
CHECK_SMOOTH="$SCRIPT_DIR/check_smooth_test.py"

# scripts/../../../../exe  ==  <cactus>/exe
CACTUS_SIM="${CACTUS_SIM:-$(cd "$SCRIPT_DIR/../../../.." 2>/dev/null && pwd)/exe/cactus_sim}"
EXE_DIR="${EXE_DIR:-$(dirname "$CACTUS_SIM")}"
NRANKS="${NRANKS:-1}"
SKIP_RUN="${SKIP_RUN:-0}"

: "${OMP_NUM_THREADS:=8}"
: "${OMP_PLACES:=cores}"
: "${OMP_PROC_BIND:=close}"
export OMP_NUM_THREADS OMP_PLACES OMP_PROC_BIND

LOG_DIR="$EXE_DIR/run_checks_logs"

fail_early() { echo "run_checks.sh: $*" >&2; exit 2; }

[[ -x "$CACTUS_SIM" ]] || fail_early "cactus binary not found or not executable: $CACTUS_SIM (build it first, or set CACTUS_SIM)"
[[ -d "$EXE_DIR"    ]] || fail_early "output/run dir does not exist: $EXE_DIR (set EXE_DIR)"
[[ -f "$CHECK_COLOR"  ]] || fail_early "missing checker: $CHECK_COLOR"
[[ -f "$CHECK_SMOOTH" ]] || fail_early "missing checker: $CHECK_SMOOTH"
mkdir -p "$LOG_DIR"

# Result accumulation.
declare -a SUMMARY=()          # human-readable "STATUS  name  detail" lines
had_regression=0
had_cantrun=0

record() {  # record <status> <name> <detail>
  local status="$1" name="$2" detail="$3"
  case "$status" in
    REGRESSION)    had_regression=1 ;;
    CANNOT-RUN)    had_cantrun=1 ;;
  esac
  SUMMARY+=("$(printf '%-11s %-26s %s' "$status" "$name" "$detail")")
}

# run_par <basename> -> populates "$EXE_DIR/<basename>/"; returns cactus rc.
#
# The two streams are kept SEPARATE (<base>.run.out / <base>.run.err) and must
# stay that way. A CCTK warning can land on BOTH: the flesh sends it to stderr
# when its level <= warning_level (default 1) and to stdout when its level <=
# logging_level -- which the -L3 below raises to 3 -- or, on any non-root rank,
# whenever level <= warning_level (Cactus/src/main/WarnLevel.c:576-688). So under
# -L3 levels 0 and 1 appear twice while 2 and 3 are stdout-only, and the two
# copies are not even byte-identical (mpiexec pty-forwards rank 0's stdout, so
# Cactus bolds the stdout copy and not the stderr one). A merged log can
# therefore be neither stream-attributed nor de-duplicated, and its duplicated
# lines read as twice the evidence.
# (The Python checkers below keep 2>&1 on purpose: they are single-stream, and
# only their exit code and their stdout RESULT line are ever consumed.)
run_par() {
  local base="$1" par="$PAR_DIR/$1.par"
  local out="$LOG_DIR/$1.run.out" err="$LOG_DIR/$1.run.err"
  [[ -f "$par" ]] || { echo "MISSING PAR: $par" >"$err"; : >"$out"; return 127; }
  if [[ "$SKIP_RUN" == "1" ]]; then
    echo "SKIP_RUN=1: reusing existing $EXE_DIR/$base" >"$out"
    : >"$err"
    return 0
  fi
  rm -rf "${EXE_DIR:?}/$base"
  ( cd "$EXE_DIR" && mpiexec -n "$NRANKS" "$CACTUS_SIM" -L3 "$par" ) >"$out" 2>"$err"
}

echo "== run_checks.sh =="
echo "binary : $CACTUS_SIM"
echo "exe dir: $EXE_DIR"
echo "ranks  : $NRANKS   (SKIP_RUN=$SKIP_RUN)"
echo "logs   : $LOG_DIR"
echo

# ---- Color correctness gate -------------------------------------------------
for base in "${COLOR_PARS[@]}"; do
  name="color:$base"
  echo "-- running $base.par ..."
  if ! run_par "$base"; then
    record "CANNOT-RUN" "$name" "cactus run aborted (see $LOG_DIR/$base.run.{out,err})"
    continue
  fi
  tsv="$EXE_DIR/$base/$COLOR_TSV"
  if [[ ! -f "$tsv" ]]; then
    record "CANNOT-RUN" "$name" "no color TSV produced: $tsv"
    continue
  fi
  clog="$LOG_DIR/$base.check.log"
  python3 "$CHECK_COLOR" "$tsv" >"$clog" 2>&1
  rc=$?
  result_line="$(grep -m1 '^RESULT:' "$clog" || true)"
  case "$rc" in
    0) record "PASS"       "$name" "${result_line:-clean}" ;;
    1) record "REGRESSION" "$name" "${result_line:-incorrect interpatch data}" ;;
    2) record "CANNOT-RUN" "$name" "checker could not run (${result_line:-broken/empty join}; see $clog)" ;;
    *) record "CANNOT-RUN" "$name" "checker crashed rc=$rc (see $clog)" ;;
  esac
done

# ---- Smooth-field Channel-1 A/B measurement ---------------------------------
name="smooth:${SMOOTH_TRIPLET[0]%_*}"   # e.g. smooth:smooth_z
declare -a smooth_tsvs=()
smooth_ok=1
for leg in "${SMOOTH_TRIPLET[@]}"; do
  echo "-- running $leg.par ..."
  if ! run_par "$leg"; then
    record "CANNOT-RUN" "$name" "leg $leg aborted (see $LOG_DIR/$leg.run.{out,err})"
    smooth_ok=0; break
  fi
  tsv="$EXE_DIR/$leg/$SMOOTH_TSV"
  if [[ ! -f "$tsv" ]]; then
    record "CANNOT-RUN" "$name" "leg $leg produced no smooth TSV: $tsv"
    smooth_ok=0; break
  fi
  smooth_tsvs+=("$tsv")
done

if [[ "$smooth_ok" == "1" ]]; then
  clog="$LOG_DIR/${name#smooth:}.check.log"
  python3 "$CHECK_SMOOTH" \
    --neumann   "${smooth_tsvs[0]}" \
    --linextrap "${smooth_tsvs[1]}" \
    --none      "${smooth_tsvs[2]}" >"$clog" 2>&1
  rc=$?
  result_line="$(grep -m1 '^RESULT:' "$clog" || true)"
  case "$rc" in
    0) record "PASS" "$name" "Channel 1 inert -- ${result_line:-no cross-patch change}" ;;
    1) record "PASS" "$name" "Channel 1 live (expected positive) -- ${result_line:-}" ;;
    2) record "CANNOT-RUN" "$name" "interior-invariant broken / checker could not run (see $clog)" ;;
    *) record "CANNOT-RUN" "$name" "checker crashed rc=$rc (see $clog)" ;;
  esac
fi

# ---- Summary ----------------------------------------------------------------
echo
echo "== summary =="
for line in "${SUMMARY[@]}"; do echo "  $line"; done
echo

if [[ "$had_regression" == "1" ]]; then
  echo "OVERALL: FAIL -- interpatch regression detected."
  exit 1
elif [[ "$had_cantrun" == "1" ]]; then
  echo "OVERALL: FAIL -- a check could not run; results are untrustworthy."
  exit 2
else
  echo "OVERALL: PASS -- all interpatch checks green."
  exit 0
fi
