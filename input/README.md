# Input structures

Place the locally approved WT PDB here, or set WT_PDB in your local configuration to an external path. The structure is intentionally not committed.

For the default V623Q/V623S workflow, the mutation builder expects:
- chain A present;
- residue 623 is VAL in the WT structure;
- a standard protein-only structure with consistent residue numbering;
- no ligand, RNA, metal, or other nonstandard component that requires a separate topology.

The runner's protein-only preparation removes HETATM records. Do not use this default path for a complex whose ligand, RNA, cofactor, or metal is mechanistically required. Check structure provenance, chain mapping, missing residues, alternate conformations, termini, and protonation assumptions before simulation.
