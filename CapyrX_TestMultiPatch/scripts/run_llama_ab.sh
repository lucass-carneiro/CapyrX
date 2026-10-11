#!/usr/bin/env bash
# Step 16 of llama_patch_system_impl.md (Phase E): the plane-wave A/B between the
# CapyrX `Llama` patch system and the real Llama reference (LlamaWaveToy /
# Thornburg04). This is a CROSS-CODE test and needs TWO binaries, so it is NOT part
# of the single-exe run_checks.sh battery -- run it with this script instead.
#
#   CapyrX side : par/llama_pwave_{16,32,64}.par   -> $CAPYRX_EXE (this checkout)
#   Llama side  : par/lwt_pwave_{16,32,64}.par     -> $LWT_EXE    (AxiDil checkout)
#
# Both evolve the same exact travelling plane wave to t=1.0 on a matched central-
# cube + 6-wedge grid; check_llama_ab.py fits the error-vs-analytic convergence
# order per band (owner points, r<2.5) and gates on both codes reaching the FD
# order and the CapyrX seam/interior ratio staying bounded (no CapyrX-only
# interface signal). See check_llama_ab.py for the full method.
#
# Env overrides:
#   CAPYRX_EXE  CapyrX binary     (default: <repo>/../../../exe/cactus_sim)
#   LWT_EXE     Llama binary      (default: /home/lsanches/etk/AxiDil/Cactus/exe/cactus_sim)
#   CAPYRX_DIR  where CapyrX runs/outputs live (default: dir of CAPYRX_EXE)
#   LWT_DIR     where Llama runs/outputs live  (default: dir of LWT_EXE)
#   RES         resolutions to run             (default: "16 32 64")
#   SKIP_RUN=1  do not run cactus; just re-check existing output
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAR_DIR="$(cd "$SCRIPT_DIR/../par" && pwd)"
CHECK="$SCRIPT_DIR/check_llama_ab.py"

CAPYRX_EXE="${CAPYRX_EXE:-$SCRIPT_DIR/../../../../exe/cactus_sim}"
LWT_EXE="${LWT_EXE:-/home/lsanches/etk/AxiDil/Cactus/exe/cactus_sim}"
CAPYRX_DIR="${CAPYRX_DIR:-$(cd "$(dirname "$CAPYRX_EXE")" && pwd)}"
LWT_DIR="${LWT_DIR:-$(cd "$(dirname "$LWT_EXE")" && pwd)}"
RES="${RES:-16 32 64}"
SKIP_RUN="${SKIP_RUN:-0}"

run_leg() { # exe rundir parfile_basename
  local exe="$1" rundir="$2" base="$3"
  if [[ "$SKIP_RUN" == 1 ]]; then
    echo "  SKIP_RUN: reusing $rundir/$base"
    return
  fi
  echo "  running $base -> $rundir"
  ( cd "$rundir" && "$exe" "$PAR_DIR/$base.par" >"$rundir/$base.ablog" 2>&1 )
}

cap_dirs=(); lwt_dirs=()
for N in $RES; do
  run_leg "$CAPYRX_EXE" "$CAPYRX_DIR" "llama_pwave_$N"
  run_leg "$LWT_EXE"    "$LWT_DIR"    "lwt_pwave_$N"
  cap_dirs+=("$CAPYRX_DIR/llama_pwave_$N")
  lwt_dirs+=("$LWT_DIR/lwt_pwave_$N")
done

python3 -I "$CHECK" --capyrx "${cap_dirs[@]}" --llama "${lwt_dirs[@]}" "$@"
