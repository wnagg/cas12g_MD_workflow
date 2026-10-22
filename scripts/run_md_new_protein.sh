#!/usr/bin/env bash
# Reusable GROMACS pipeline for standard protein-only systems.
# Usage: bash scripts/run_md_new_protein.sh /path/to/run.env
set -Eeuo pipefail
if [[ $# -ne 1 || ! -f "$1" ]]; then
  echo "Usage: bash scripts/run_md_new_protein.sh /path/to/run.env" >&2
  exit 2
fi
CONFIG_FILE="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
# shellcheck disable=SC1090
source "$CONFIG_FILE"

[[ -n "$GMX" && -n "$INPUT_PDB" && -n "$WORKDIR" && -n "$SYSTEM_NAME" ]] || {
  echo "ERROR: GMX, INPUT_PDB, WORKDIR and SYSTEM_NAME are required in config." >&2; exit 1;
}
if [[ "$GMX" == */* ]]; then
  [[ -x "$GMX" ]] || { echo "ERROR: GROMACS not executable: $GMX" >&2; exit 1; }
else
  GMX="$(command -v "$GMX")" || { echo "ERROR: GROMACS command not found." >&2; exit 1; }
fi
if [[ "$INPUT_PDB" != /* ]]; then INPUT_PDB="$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$INPUT_PDB"; fi
[[ -s "$INPUT_PDB" ]] || { echo "ERROR: input PDB missing: $INPUT_PDB" >&2; exit 1; }

mkdir -p "$WORKDIR"
cd "$WORKDIR"
mkdir -p logs .done analysis
exec > >(tee -a pipeline.log) 2>&1

signature="$(
  { printf '%s\n' "$SYSTEM_NAME" "$FORCEFIELD" "$WATER" "$BOX_TYPE" "$BOX_DIST_NM" "$SALT_MOLAR" "$TEMP_K" "$PRESSURE_BAR" "$NVT_PS" "$NPT_CA_PS" "$NPT_NOCA_PS" "$MD_NS" "$DT_PS" "$XTC_INTERVAL_PS" "$ENERGY_INTERVAL_PS" "$LOG_INTERVAL_PS" "$VELOCITY_SEED" "$CA_FC" "$MD_DEFFNM"; sha256sum "$INPUT_PDB"; } | sha256sum | awk '{print $1}'
)"
if [[ -s .run_signature ]]; then
  [[ "$(cat .run_signature)" == "$signature" ]] || {
    echo "ERROR: input or MD settings differ from this WORKDIR's original signature." >&2
    echo "Use a fresh WORKDIR for changed settings; restore the original config to resume." >&2
    exit 1
  }
else
  if find .done -maxdepth 1 -name '*.done' -print -quit | grep -q . || [[ -s "$MD_DEFFNM.tpr" || -s "$MD_DEFFNM.cpt" ]]; then
    echo "ERROR: existing run files found without a signature; refusing to mix runs." >&2
    exit 1
  fi
  printf '%s\n' "$signature" > .run_signature
fi

calc_steps() { awk -v ps="$1" -v dt="$DT_PS" 'BEGIN { if (dt <= 0) exit 1; printf "%.0f\n", ps/dt }'; }
NVT_STEPS="$(calc_steps "$NVT_PS")"
NPT_CA_STEPS="$(calc_steps "$NPT_CA_PS")"
NPT_NOCA_STEPS="$(calc_steps "$NPT_NOCA_PS")"
MD_STEPS="$(awk -v ns="$MD_NS" -v dt="$DT_PS" 'BEGIN {printf "%.0f\n", ns*1000/dt}')"
XTC_NSTEPS="$(calc_steps "$XTC_INTERVAL_PS")"
ENERGY_NSTEPS="$(calc_steps "$ENERGY_INTERVAL_PS")"
LOG_NSTEPS="$(calc_steps "$LOG_INTERVAL_PS")"
EQUIL_XTC_NSTEPS="$(calc_steps 10)"

run_stage() {
  local name="$1"; shift
  if [[ -f ".done/$name.done" ]]; then echo "[SKIP] $name already completed"; return; fi
  echo "[START] $name at $(date -Is)"
  "$@" 2>&1 | tee "logs/$name.out"
  touch ".done/$name.done"
  echo "[DONE] $name at $(date -Is)"
}

write_mdp_files() {
  cat > ions.mdp <<EOF
integrator = steep
emtol = 1000.0
emstep = 0.01
nsteps = 5000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
pbc = xyz
EOF
  cat > minim.mdp <<EOF
integrator = steep
nsteps = 50000
emtol = 1000.0
emstep = 0.01
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
constraints = none
pbc = xyz
EOF
  cat > nvt.mdp <<EOF
define = -DPOSRES_CA
integrator = md
dt = $DT_PS
nsteps = $NVT_STEPS
continuation = no
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
pbc = xyz
constraints = h-bonds
constraint_algorithm = lincs
lincs_iter = 1
lincs_order = 4
tcoupl = V-rescale
tc-grps = Protein Water_and_ions
tau_t = 0.1 0.1
ref_t = $TEMP_K $TEMP_K
annealing = single single
annealing_npoints = 2 2
annealing_time = 0 $NVT_PS 0 $NVT_PS
annealing_temp = 50 $TEMP_K 50 $TEMP_K
pcoupl = no
gen_vel = yes
gen_temp = 50
gen_seed = $VELOCITY_SEED
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = $EQUIL_XTC_NSTEPS
nstenergy = $ENERGY_NSTEPS
nstlog = $LOG_NSTEPS
EOF
  cat > npt_ca.mdp <<EOF
define = -DPOSRES_CA
integrator = md
dt = $DT_PS
nsteps = $NPT_CA_STEPS
continuation = yes
gen_vel = no
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
pbc = xyz
constraints = h-bonds
constraint_algorithm = lincs
lincs_iter = 1
lincs_order = 4
tcoupl = V-rescale
tc-grps = Protein Water_and_ions
tau_t = 0.1 0.1
ref_t = $TEMP_K $TEMP_K
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = $PRESSURE_BAR
compressibility = 4.5e-5
refcoord_scaling = com
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = $EQUIL_XTC_NSTEPS
nstenergy = $ENERGY_NSTEPS
nstlog = $LOG_NSTEPS
EOF
  cat > npt_noca.mdp <<EOF
integrator = md
dt = $DT_PS
nsteps = $NPT_NOCA_STEPS
continuation = yes
gen_vel = no
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
pbc = xyz
constraints = h-bonds
constraint_algorithm = lincs
lincs_iter = 1
lincs_order = 4
tcoupl = V-rescale
tc-grps = Protein Water_and_ions
tau_t = 0.1 0.1
ref_t = $TEMP_K $TEMP_K
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 2.0
ref_p = $PRESSURE_BAR
compressibility = 4.5e-5
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = $EQUIL_XTC_NSTEPS
nstenergy = $ENERGY_NSTEPS
nstlog = $LOG_NSTEPS
EOF
  cat > production.mdp <<EOF
integrator = md
dt = $DT_PS
nsteps = $MD_STEPS
continuation = yes
gen_vel = no
cutoff-scheme = Verlet
nstlist = 100
rlist = 1.0
coulombtype = PME
rcoulomb = 1.0
vdwtype = Cut-off
rvdw = 1.0
DispCorr = EnerPres
pbc = xyz
constraints = h-bonds
constraint_algorithm = lincs
lincs_iter = 1
lincs_order = 4
tcoupl = V-rescale
tc-grps = Protein Water_and_ions
tau_t = 0.1 0.1
ref_t = $TEMP_K $TEMP_K
pcoupl = C-rescale
pcoupltype = isotropic
tau_p = 5.0
ref_p = $PRESSURE_BAR
compressibility = 4.5e-5
comm-mode = Linear
nstcomm = 100
nstxout = 0
nstvout = 0
nstfout = 0
nstxout-compressed = $XTC_NSTEPS
nstenergy = $ENERGY_NSTEPS
nstlog = $LOG_NSTEPS
EOF
}

prepare_input() {
  if [[ "$PREPARE_PROTEIN_ONLY" == yes ]]; then
    awk '/^ATOM  / || /^TER   / || /^END   /' "$INPUT_PDB" > input_prepared.pdb
    grep -q '^END' input_prepared.pdb || echo END >> input_prepared.pdb
  else
    cp "$INPUT_PDB" input_prepared.pdb
  fi
  grep -q '^ATOM  ' input_prepared.pdb || { echo "ERROR: no ATOM records in prepared PDB." >&2; return 1; }
  echo "Prepared input: $INPUT_PDB"
  sha256sum "$INPUT_PDB"
}

topology_solvate_ions() {
  "$GMX" pdb2gmx -f input_prepared.pdb -o processed.gro -p topol.top -ff "$FORCEFIELD" -water "$WATER" -ignh
  "$GMX" editconf -f processed.gro -o newbox.gro -c -d "$BOX_DIST_NM" -bt "$BOX_TYPE"
  "$GMX" solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top
  "$GMX" grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn "$ION_MAXWARN"
  printf "SOL\n" | "$GMX" genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc "$SALT_MOLAR"
}

run_em() {
  "$GMX" grompp -f minim.mdp -c solvated_ions.gro -p topol.top -o em.tpr -maxwarn "$MAXWARN"
  "$GMX" mdrun -deffnm em -v -nb cpu -pme cpu
  grep -E "Steepest Descents converged|Fmax|Potential Energy" em.log | tail -n 20 || true
}

insert_posres_block() {
  python3 - "$1" "$2" <<'PY'
from pathlib import Path
import sys
itp, macro = sys.argv[1], sys.argv[2]
p = Path("topol.top")
text = p.read_text()
if f'#include "{itp}"' in text:
    print(f"{itp} already included")
    raise SystemExit(0)
lines = text.splitlines(True)
insert_at = next((i for i, line in enumerate(lines)
                  if any(token in line.lower() for token in ("tip3p.itp", "spc.itp", "spce.itp"))), None)
if insert_at is None:
    insert_at = next((i for i, line in enumerate(lines)
                      if line.strip().lower() == "[ system ]"), len(lines))
block = f"#ifdef {macro}\n#include \"{itp}\"\n#endif\n\n"
lines.insert(insert_at, block)
p.write_text("".join(lines))
PY
}

make_ca_restraint() {
  "$GMX" select -s em.tpr -select 'name CA and group "Protein"' -on index_ca.ndx
  printf "0\n" | "$GMX" genrestr -f em.gro -n index_ca.ndx -o posre_Ca.itp -fc "$CA_FC" "$CA_FC" "$CA_FC"
  insert_posres_block posre_Ca.itp POSRES_CA
}

run_mdrun_resumable() {
  local name="$1"; shift
  if [[ -s "$name.cpt" ]]; then
    echo "[RESUME] $name from checkpoint"
    "$GMX" mdrun -deffnm "$name" -cpi "$name.cpt" -append -v -cpt "$CHECKPOINT_MIN" "$@"
  else
    if [[ -s "$name.log" || -s "$name.xtc" || -s "$name.edr" ]]; then
      echo "ERROR: partial outputs exist for $name but no checkpoint was found; inspect manually." >&2
      return 1
    fi
    "$GMX" mdrun -deffnm "$name" -v -cpt "$CHECKPOINT_MIN" "$@"
  fi
}
run_nvt() {
  "$GMX" grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn "$MAXWARN"
  run_mdrun_resumable nvt -nb cpu -pme cpu
}
run_npt_ca() {
  "$GMX" grompp -f npt_ca.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -o npt_ca.tpr -maxwarn "$MAXWARN"
  run_mdrun_resumable npt_ca -nb cpu -pme cpu
}
run_npt_noca() {
  "$GMX" grompp -f npt_noca.mdp -c npt_ca.gro -t npt_ca.cpt -p topol.top -o npt_noca.tpr -maxwarn "$MAXWARN"
  run_mdrun_resumable npt_noca -nb cpu -pme cpu
}
prepare_production() {
  "$GMX" grompp -f production.mdp -c npt_noca.gro -t npt_noca.cpt -p topol.top -o "$MD_DEFFNM.tpr" -maxwarn "$MAXWARN"
}
run_production() {
  if [[ -s "$MD_DEFFNM.cpt" ]]; then
    "$GMX" mdrun -deffnm "$MD_DEFFNM" -cpi "$MD_DEFFNM.cpt" -append -v -cpt "$CHECKPOINT_MIN" -nb gpu -pme gpu
  else
    if [[ -s "$MD_DEFFNM.log" || -s "$MD_DEFFNM.xtc" || -s "$MD_DEFFNM.edr" ]]; then
      echo "ERROR: partial production output exists without a checkpoint; refusing to overwrite." >&2
      return 1
    fi
    "$GMX" mdrun -deffnm "$MD_DEFFNM" -v -cpt "$CHECKPOINT_MIN" -nb gpu -pme gpu
  fi
}

write_mdp_files
run_stage 01_prepare_input prepare_input
run_stage 02_topology_solvation_ions topology_solvate_ions
run_stage 03_energy_minimization run_em
run_stage 04_ca_restraints make_ca_restraint
run_stage 05_nvt run_nvt
run_stage 06_npt_ca run_npt_ca
run_stage 07_npt_noca run_npt_noca
run_stage 08_prepare_production prepare_production
run_stage 09_production run_production
touch .done/complete.done
echo "[DONE] $SYSTEM_NAME: requested production run finished."
echo "Output directory: $WORKDIR"
echo "Production trajectory: $WORKDIR/$MD_DEFFNM.xtc"
