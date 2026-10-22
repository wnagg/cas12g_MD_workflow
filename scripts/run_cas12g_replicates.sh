#!/usr/bin/env bash
# Cas12g panel runner: WT/V623Q/V623S, three independent 400 ns replicas.
# One server: REPLICATE_IDS="1 2 3". Three servers: set one ID per server.
set -Eeuo pipefail

if [[ $# -lt 1 || ! -f "$1" ]]; then
  echo "Usage: bash scripts/run_cas12g_replicates.sh CONFIG.env [--build-only|--pilot-only|--run-only|--analyze-only]" >&2
  exit 2
fi
CONFIG="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
MODE=all
if [[ $# -ge 2 ]]; then MODE="$2"; fi
case "$MODE" in all|--build-only|--pilot-only|--run-only|--analyze-only) ;; *) echo "Unknown mode: $MODE" >&2; exit 2;; esac

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE="$ROOT/scripts/run_md_new_protein.sh"
BUILDER="$ROOT/scripts/build_cas12g_mutants.py"
ANALYZER="$ROOT/analysis/analyze_cas12g_md.py"
AGGREGATOR="$ROOT/scripts/compare_replicates.py"
# shellcheck disable=SC1090
source "$CONFIG"

for name in GMX PYMOL PYTHON WT_PDB CHAIN WORKDIR VARIANTS REPLICATE_IDS BASE_SEED FORCEFIELD WATER BOX_TYPE BOX_DIST_NM SALT_MOLAR TEMP_K PRESSURE_BAR PREPARE_PROTEIN_ONLY NVT_PS NPT_CA_PS NPT_NOCA_PS MD_NS DT_PS XTC_INTERVAL_PS ENERGY_INTERVAL_PS LOG_INTERVAL_PS CHECKPOINT_MIN MD_DEFFNM CA_FC MAXWARN ION_MAXWARN ANALYSIS_BEGIN_NS FEL_BINS; do
  declare -p "$name" >/dev/null 2>&1 || { echo "ERROR: missing config variable: $name" >&2; exit 1; }
done
if [[ "$GMX" != */* ]]; then GMX="$(command -v "$GMX")"; fi
if [[ "$PYMOL" != */* ]]; then PYMOL="$(command -v "$PYMOL")"; fi
if [[ "$PYTHON" != */* ]]; then PYTHON="$(command -v "$PYTHON")"; fi
[[ -x "$GMX" ]] || { echo "ERROR: GROMACS not executable: $GMX" >&2; exit 1; }
[[ -x "$PYMOL" ]] || { echo "ERROR: PyMOL not executable: $PYMOL" >&2; exit 1; }
[[ -x "$PYTHON" ]] || { echo "ERROR: Python not executable: $PYTHON" >&2; exit 1; }
if [[ "$WT_PDB" != /* ]]; then WT_PDB="$ROOT/$WT_PDB"; fi
[[ -s "$WT_PDB" ]] || { echo "ERROR: WT_PDB missing: $WT_PDB" >&2; exit 1; }
mkdir -p "$WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd)"
STRUCT_DIR="$WORKDIR/structures"
RUN_ROOT="$WORKDIR/runs"
ANALYSIS_DIR="$WORKDIR/analysis"
CONFIG_DIR="$WORKDIR/.configs"
mkdir -p "$STRUCT_DIR" "$RUN_ROOT" "$ANALYSIS_DIR" "$CONFIG_DIR"
exec 9>"$WORKDIR/.cas12g_replicates.lock"
flock -n 9 || { echo "ERROR: another workflow is using $WORKDIR." >&2; exit 1; }

seen=" "
for rep in $REPLICATE_IDS; do
  case "$rep" in 1|2|3) ;; *) echo "ERROR: REPLICATE_IDS may contain only 1, 2, 3." >&2; exit 1;; esac
  case "$seen" in *" $rep "*) echo "ERROR: duplicate replica ID: $rep" >&2; exit 1;; esac
  seen="$seen$rep "
done
[[ "$BASE_SEED" =~ ^[0-9]+$ ]] || { echo "ERROR: BASE_SEED must be a positive integer." >&2; exit 1; }
[[ "$VARIANTS" == "WT V623Q V623S" ]] || { echo "ERROR: this workflow expects VARIANTS=\"WT V623Q V623S\"." >&2; exit 1; }
"$PYTHON" -c 'import numpy,pandas,matplotlib' || { echo "ERROR: install requirements.txt first." >&2; exit 1; }
[[ -s "$PIPELINE" && -s "$BUILDER" && -s "$ANALYZER" && -s "$AGGREGATOR" ]] || {
  echo "ERROR: one or more workflow files are missing." >&2; exit 1;
}

build_structures() {
  if find "$RUN_ROOT" -name '*.tpr' -print -quit | grep -q .; then
    echo "ERROR: run TPRs already exist. Use a fresh WORKDIR before rebuilding structures." >&2
    exit 1
  fi
  export CAS12G_WT_PDB="$WT_PDB" CAS12G_STRUCT_DIR="$STRUCT_DIR" CAS12G_CHAIN="$CHAIN"
  (cd "$ROOT" && "$PYMOL" -cq -d "run scripts/build_cas12g_mutants.py" -d "quit") > "$WORKDIR/build_mutants.log" 2>&1 || {
    tail -n 80 "$WORKDIR/build_mutants.log" >&2; exit 1;
  }
  grep -q '\[DONE\] All mutants were built independently from WT' "$WORKDIR/build_mutants.log" || {
    echo "ERROR: builder did not report successful completion; see $WORKDIR/build_mutants.log" >&2; exit 1;
  }
  for tag in WT V623Q V623S; do
    [[ -s "$STRUCT_DIR/$tag.pdb" ]] || { echo "ERROR: missing built structure $STRUCT_DIR/$tag.pdb" >&2; exit 1; }
  done
  echo "[OK] Built WT-derived structures. QC: $STRUCT_DIR/initial_models_qc.csv"
}

write_run_config() {
  local input_pdb="$1" run_dir="$2" system="$3" ns="$4" seed="$5" deffnm="$6" out="$7"
  {
    printf 'GMX=%q\n' "$GMX"
    printf 'INPUT_PDB=%q\n' "$input_pdb"
    printf 'WORKDIR=%q\n' "$run_dir"
    printf 'SYSTEM_NAME=%q\n' "$system"
    printf 'FORCEFIELD=%q\n' "$FORCEFIELD"
    printf 'WATER=%q\n' "$WATER"
    printf 'BOX_TYPE=%q\n' "$BOX_TYPE"
    printf 'BOX_DIST_NM=%q\n' "$BOX_DIST_NM"
    printf 'SALT_MOLAR=%q\n' "$SALT_MOLAR"
    printf 'TEMP_K=%q\n' "$TEMP_K"
    printf 'PRESSURE_BAR=%q\n' "$PRESSURE_BAR"
    printf 'PREPARE_PROTEIN_ONLY=%q\n' "$PREPARE_PROTEIN_ONLY"
    printf 'NVT_PS=%q\n' "$NVT_PS"
    printf 'NPT_CA_PS=%q\n' "$NPT_CA_PS"
    printf 'NPT_NOCA_PS=%q\n' "$NPT_NOCA_PS"
    printf 'MD_NS=%q\n' "$ns"
    printf 'DT_PS=%q\n' "$DT_PS"
    printf 'XTC_INTERVAL_PS=%q\n' "$XTC_INTERVAL_PS"
    printf 'ENERGY_INTERVAL_PS=%q\n' "$ENERGY_INTERVAL_PS"
    printf 'LOG_INTERVAL_PS=%q\n' "$LOG_INTERVAL_PS"
    printf 'CHECKPOINT_MIN=%q\n' "$CHECKPOINT_MIN"
    printf 'VELOCITY_SEED=%q\n' "$seed"
    printf 'MD_DEFFNM=%q\n' "$deffnm"
    printf 'CA_FC=%q\n' "$CA_FC"
    printf 'MAXWARN=%q\n' "$MAXWARN"
    printf 'ION_MAXWARN=%q\n' "$ION_MAXWARN"
  } > "$out"
}

run_one() {
  local rep="$1" tag="$2" ns="$3" seed="$4" deffnm="$5"
  local pdb="$STRUCT_DIR/$tag.pdb"
  local run_dir="$RUN_ROOT/rep$rep/$tag"
  local cfg="$CONFIG_DIR/rep$rep"_"$tag".env
  [[ -s "$pdb" ]] || { echo "ERROR: missing $pdb. Run --build-only first." >&2; exit 1; }
  mkdir -p "$run_dir"
  write_run_config "$pdb" "$run_dir" "Cas12g_"$tag"_rep"$rep "$ns" "$seed" "$deffnm" "$cfg"
  cat > "$run_dir/run_metadata.txt" <<EOF
replicate_id=$rep
variant=$tag
velocity_seed=$seed
production_ns=$ns
timestep_ps=$DT_PS
xtc_interval_ps=$XTC_INTERVAL_PS
forcefield=$FORCEFIELD
water=$WATER
salt_molar=$SALT_MOLAR
temperature_K=$TEMP_K
pressure_bar=$PRESSURE_BAR
EOF
  echo "========== replica $rep / $tag: $ns ns; seed=$seed =========="
  bash "$PIPELINE" "$cfg"
  "$PYTHON" "$ANALYZER" --workdir "$run_dir" --tag "$tag" --gmx "$GMX" --deffnm "$deffnm" --begin-ns "$ANALYSIS_BEGIN_NS" --temp-k "$TEMP_K" --fel-bins "$FEL_BINS" --outdir "$run_dir/analysis"
}

if [[ "$MODE" != --run-only && "$MODE" != --analyze-only ]]; then
  build_needed=no
  for tag in WT V623Q V623S; do [[ -s "$STRUCT_DIR/$tag.pdb" ]] || build_needed=yes; done
  if [[ "$REBUILD_MUTANTS" == yes ]]; then build_needed=yes; fi
  if [[ "$build_needed" == yes ]]; then build_structures; else echo "[OK] Existing structure set found."; fi
fi
if [[ "$MODE" == --build-only ]]; then exit 0; fi

if [[ "$MODE" == --pilot-only ]]; then
  pilot_dir="$WORKDIR/pilot/WT"
  pilot_cfg="$CONFIG_DIR/pilot_WT.env"
  pilot_seed=$((BASE_SEED + 9999))
  write_run_config "$STRUCT_DIR/WT.pdb" "$pilot_dir" "Cas12g_WT_pilot" 10 "$pilot_seed" mdpilot "$pilot_cfg"
  bash "$PIPELINE" "$pilot_cfg"
  "$PYTHON" "$ANALYZER" --workdir "$pilot_dir" --tag WT_pilot --gmx "$GMX" --deffnm mdpilot --begin-ns 2 --temp-k "$TEMP_K" --fel-bins 15 --outdir "$pilot_dir/analysis"
  echo "[DONE] WT pilot finished. Review $pilot_dir/logs/09_production.out before launching production."
  exit 0
fi

for rep in $REPLICATE_IDS; do
  idx=0
  for tag in WT V623Q V623S; do
    idx=$((idx + 1))
    seed=$((BASE_SEED + rep * 1000 + idx))
    deffnm="$MD_DEFFNM"
    if [[ "$MODE" == --analyze-only ]]; then
      run_dir="$RUN_ROOT/rep$rep/$tag"
      "$PYTHON" "$ANALYZER" --workdir "$run_dir" --tag "$tag" --gmx "$GMX" --deffnm "$deffnm" --begin-ns "$ANALYSIS_BEGIN_NS" --temp-k "$TEMP_K" --fel-bins "$FEL_BINS" --outdir "$run_dir/analysis"
    else
      run_one "$rep" "$tag" "$MD_NS" "$seed" "$deffnm"
    fi
  done
done
"$PYTHON" "$AGGREGATOR" --roots "$WORKDIR" --replicate-ids $REPLICATE_IDS --outdir "$ANALYSIS_DIR/replicate_comparison" --variants WT V623Q V623S --begin-ns "$ANALYSIS_BEGIN_NS"
echo "[DONE] Workflow mode $MODE complete. Results: $ANALYSIS_DIR"
