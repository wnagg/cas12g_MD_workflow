# Input structure

Provide the WT Cas12g PDB structure locally and set its path in the configuration file. The structure is not distributed in this repository.

The current mutation builder expects chain A residue 623 to be valine and generates V623Q and V623S from the WT structure. Check the chain assignment, residue numbering, missing residues, and structural components before running.

The default preparation retains protein ATOM records and is intended for protein-only systems. If a ligand, RNA, metal, or other non-protein component is required for the system being studied, prepare the complete complex and its compatible topology rather than using the default protein-only workflow.
