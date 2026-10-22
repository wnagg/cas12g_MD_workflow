# Cas12g MD Workflow

A reproducible, configuration-driven GROMACS workflow for the Cas12g WT/V623Q/V623S comparison.

## Study design

- Systems: WT, V623Q, and V623S.
- Replicates: three independent replicas per system (nine production runs total).
- Production target: 400 ns per replica (3.6 microseconds aggregate production time).
- Pilot: a separate 10 ns WT run to check stability and estimate throughput before launching production.
- Starting mutants are built independently from the same WT PDB; neither mutant is used as the template for the other.
- Intended baseline: protein-only/apo simulation, AMBER99SB-ILDN, TIP3P water, dodecahedral box, at least 1.0 nm solute-to-box margin, neutralization plus 0.15 M NaCl, 300 K and 1 bar.

This repository contains workflow code and documentation, not the source structure or simulation trajectories. Supply a structure that you are permitted to use locally. The default mutation check expects chain A residue 623 to be VAL. Confirm residue numbering and biological context before running.

## Pipeline overview

1. **Preflight**: check GROMACS, PyMOL, Python dependencies, input PDB, chain A and residue 623 identity, GPU visibility, and free disk space.
2. **Build structures**: preserve WT and independently generate V623Q and V623S from WT; verify residue identities.
3. **10 ns pilot**: prepare the protein-only system, add solvent and ions, minimize, equilibrate, and run a short WT production trajectory in a separate directory.
4. **Production**: run each selected replica/system with a unique velocity seed. The default schedule is 100 ps NVT, 1 ns C-alpha-restrained NPT, 500 ps unrestrained NPT, then 400 ns production.
5. **Analysis**: calculate per-trajectory RMSD, radius of gyration and C-alpha RMSF; generate an RMSD–Rg occupancy/FEL summary and cross-replica plots when the analysis module is installed.

## Quick start

Clone the repository on the compute server and create a local configuration:

~~~bash
git clone https://github.com/wnagg/cas12g_MD_workflow.git
cd cas12g_MD_workflow
cp env/Cas12g_3rep_400ns.env.example env/Cas12g_3rep_400ns.env
~~~

Edit the local config and set at least:
- GMX: full path to the GROMACS executable, or the command name if it is on PATH.
- PYMOL: executable used for structure mutation.
- WT_PDB: path to the local WT PDB (not stored in this repository).
- WORKDIR: persistent results directory with sufficient free space.
- REPLICATE_IDS: 1 2 3 for one server, or one ID per server for a three-server split.

Install analysis dependencies in the intended Python environment and check the setup:

~~~bash
python3 -m pip install -r requirements.txt
bash scripts/preflight.sh env/Cas12g_3rep_400ns.env
~~~

Build structures and inspect the QC output before MD:

~~~bash
bash scripts/run_cas12g_replicates.sh env/Cas12g_3rep_400ns.env --build-only
~~~

Run the 10 ns WT pilot separately:

~~~bash
bash scripts/run_cas12g_replicates.sh env/Cas12g_3rep_400ns.env --pilot-only
~~~

Inspect the pilot log and confirm that there are no fatal errors or persistent LINCS warnings, and that temperature, pressure, density, and performance are reasonable. Estimate total wall time from measured ns/day before starting production.

Launch the production workflow only after the pilot passes:

~~~bash
nohup bash scripts/run_cas12g_replicates.sh env/Cas12g_3rep_400ns.env --run-only > cas12g_md.nohup.out 2>&1 &
echo $! > cas12g_md.pid
tail -f cas12g_md.nohup.out
~~~

The workflow is sequential on a given server to avoid GPU oversubscription. For three servers, use the same WT structure, software/force-field settings, BASE_SEED and WORKDIR naming convention, but set REPLICATE_IDS to 1, 2, or 3 on the respective server. Use distinct output directories per server or merge the exported results afterward. Do not run the same replica ID twice in a shared output directory.

## Important operational rules

- Keep local configs, input structures, logs, checkpoints, trajectories and large result files out of Git.
- Do not change force field, starting structure, timestep, production length or seed scheme inside an existing run directory. Use a new WORKDIR for a changed protocol.
- A restart must use the same TPR, checkpoint and run settings. Do not delete or replace a checkpoint to force a restart.
- Do not delete checkpoint/energy/topology/trajectory files until the result has been validated and backed up.
- C-alpha restraints are used only during the specified equilibration stages; production is unrestrained.
- A 10 ns pilot is a technical check, not evidence of equilibration or biological convergence.
- RMSD–Rg FEL plots are descriptive occupancy summaries, not proof of thermodynamic convergence. Compare independent replicas and block trends.
- This is a protein-only workflow. It is not suitable for interpreting RNA-, ligand-, or metal-dependent mechanisms unless the starting model and topology are changed to represent those components correctly.

## Repository layout

~~~text
env/       portable configuration template (copy to a local .env file)
input/     instructions for local-only starting structures
scripts/   preflight, mutant construction, GROMACS runner and replica orchestration
analysis/  trajectory analysis and comparison
docs/      method notes and reproducibility checklist
~~~

## Validation status

The repository documents the intended workflow and reusable scripts. GROMACS, PyMOL and GPU execution must be validated on the target server; passing shell/Python syntax checks does not establish successful end-to-end simulation. Record the actual GROMACS version, GPU model, wall time, warnings and completed simulation length in each run manifest.
