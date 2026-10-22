# Cas12g Molecular Dynamics Workflow

This repository contains scripts for preparing and running molecular dynamics (MD) simulations of wild-type (WT) Cas12g and mutant proteins derived from the WT structure. The structure-building script generates the specified point mutants independently from the WT model.

## Workflow

1. Check the input structure and local software environment.
2. Prepare the WT structure and WT-derived mutant structures.
3. Generate GROMACS topologies, solvate the systems, and add ions.
4. Perform energy minimization and equilibration, followed by production MD.
5. Optionally calculate structural summaries from the resulting trajectories.

The current structure builder is configured for chain A, residue 623, with V623Q and V623S substitutions. The input PDB is not included; provide the WT structure used in the study and verify its chain and residue numbering before running the scripts.

## Default simulation settings

- Force field: AMBER99SB-ILDN
- Water model: TIP3P
- Periodic box: dodecahedron, with a 1.0 nm solute-to-box margin
- Ionic conditions: neutralization and 0.15 M NaCl
- Temperature and pressure: 300 K and 1 bar
- Integration time step: 2 fs
- Equilibration: 100 ps NVT, 1 ns C-alpha-restrained NPT, and 500 ps unrestrained NPT
- Production trajectory: 400 ns; coordinates saved every 100 ps

These are the settings in the supplied configuration template. Confirm that they match the protocol used to generate the data before citing them in a manuscript.

## Requirements

- GROMACS
- PyMOL (for building mutant structures)
- Python 3 with the packages listed in `requirements.txt`

## Run

Clone the repository and create a local configuration file:

```bash
git clone https://github.com/wnagg/cas12g_MD_workflow.git
cd cas12g_MD_workflow
cp env/Cas12g_MD.env.example env/Cas12g_MD.env
```

Edit `WT_PDB`, `WORKDIR`, and the executable paths in the local configuration. Then run:

```bash
python3 -m pip install -r requirements.txt
bash scripts/preflight.sh env/Cas12g_MD.env
bash scripts/run_cas12g_md.sh env/Cas12g_MD.env --build-only
bash scripts/run_cas12g_md.sh env/Cas12g_MD.env --pilot-only
```

After checking the pilot output and confirming that the settings match the intended protocol, launch the production workflow:

```bash
bash scripts/run_cas12g_md.sh env/Cas12g_MD.env --run-only
```

The workflow uses a protein-only preparation by default. A system containing RNA, a ligand, a metal, or another required non-protein component needs an appropriately prepared structure and compatible topology; do not use the default protein-only preparation for such a complex.

Record the software versions and any deviations from the configuration used for the simulations.
