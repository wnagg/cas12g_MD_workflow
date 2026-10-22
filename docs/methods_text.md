# Methods template: Cas12g molecular dynamics

The WT, V623Q, and V623S Cas12g systems are prepared independently from a common WT starting structure. V623Q and V623S are generated as separate single-point mutants from WT and checked for the expected residue identity before system preparation. The baseline workflow is protein-only/apo; systems containing RNA, ligands, cofactors, or metals require an explicitly prepared complex and compatible topology and are outside this default configuration.

Each system is solvated in a dodecahedral box with a minimum 1.0 nm solute-to-box margin. The default force field and water model are AMBER99SB-ILDN and TIP3P. The system is neutralized and NaCl is added to a nominal concentration of 0.15 M. The equilibration schedule is 100 ps NVT, 1 ns C-alpha-restrained NPT, and 500 ps unrestrained NPT at 300 K and 1 bar. Production uses a 2 fs integration step and the configured trajectory output interval. The intended production target for the comparison panel is 400 ns per independent replica, with three replicas per system and independent initial velocity seeds.

A separate 10 ns WT pilot should be completed first to validate the local GROMACS/PyMOL installation, check numerical stability, inspect temperature/pressure/density and LINCS messages, and estimate throughput. The pilot is a technical validation step, not a convergence test.

Report the actual software versions, hardware, system atom counts, random seeds, warnings, completed production time, and any deviations from this template. Analyze RMSD, radius of gyration, C-alpha RMSF and RMSD–Rg occupancy distributions per trajectory; evaluate replicate consistency and block trends before drawing mechanistic conclusions. An occupancy-derived free-energy-like surface is descriptive and should not be presented as proof of thermodynamic convergence.

This text is a template. Replace any statement that does not match the final executed protocol, input model, or saved run metadata before using it in a manuscript.
