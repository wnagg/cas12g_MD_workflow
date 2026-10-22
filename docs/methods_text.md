# Cas12g molecular dynamics protocol

The scripts prepare WT Cas12g and mutant protein structures derived from a common WT starting structure. The current structure-building script generates the V623Q and V623S point mutants separately from WT and checks the residue identity at the target position before system preparation.

The default GROMACS setup uses the AMBER99SB-ILDN force field and TIP3P water in a dodecahedral periodic box with a 1.0 nm solute-to-box margin. The system is neutralized and supplemented with NaCl to a nominal concentration of 0.15 M. After energy minimization, the default equilibration protocol consists of 100 ps NVT, 1 ns C-alpha-restrained NPT, and 500 ps unrestrained NPT at 300 K and 1 bar. Production MD uses a 2 fs time step and a configured duration of 400 ns, with coordinates written every 100 ps.

This description reflects the configuration template in the repository. Before using it in a manuscript, verify the settings against the actual input structures, software versions, and simulation parameters used to generate the reported results.
