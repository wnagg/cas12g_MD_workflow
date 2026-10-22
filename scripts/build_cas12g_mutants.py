#!/usr/bin/env python3
"""Build WT, V623Q and V623S structures independently from one WT PDB inside PyMOL."""
from __future__ import annotations

import csv
import os
import shutil
from pathlib import Path
from pymol import cmd

MUTATIONS = {
    "V623Q": (623, "VAL", "GLN"),
    "V623S": (623, "VAL", "SER"),
}


def residue_at(path: Path, chain: str, resi: int) -> str | None:
    found = set()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM  ") or len(line) < 27:
            continue
        if line[12:16].strip() != "CA":
            continue
        if (line[21].strip() or "_") != (chain or "_"):
            continue
        if line[22:26].strip() != str(resi) or line[26].strip():
            continue
        found.add(line[17:20].strip().upper())
    if len(found) > 1:
        raise RuntimeError(f"Ambiguous residue identity at {chain}:{resi}: {found}")
    return next(iter(found)) if found else None


def main() -> None:
    wt_env = os.environ.get("CAS12G_WT_PDB", "").strip()
    out_env = os.environ.get("CAS12G_STRUCT_DIR", "").strip()
    chain = os.environ.get("CAS12G_CHAIN", "A").strip() or "A"
    if not wt_env or not out_env:
        raise RuntimeError("Set CAS12G_WT_PDB and CAS12G_STRUCT_DIR.")
    wt_path = Path(wt_env).expanduser().resolve()
    out_dir = Path(out_env).expanduser().resolve()
    if not wt_path.is_file():
        raise FileNotFoundError(wt_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    for tag, (resi, expected, _) in MUTATIONS.items():
        observed = residue_at(wt_path, chain, resi)
        if observed != expected:
            raise RuntimeError(
                f"{tag}: expected WT {expected} at chain {chain}, residue {resi}; found {observed!r}. "
                "Check chain and residue numbering."
            )

    shutil.copy2(wt_path, out_dir / "WT.pdb")
    rows = [{"variant": "WT", "chain": chain, "residue": 623,
             "expected_residue": "VAL", "observed_residue": residue_at(wt_path, chain, 623)}]
    cmd.feedback("disable", "all", "actions")

    for tag, (resi, expected, target) in MUTATIONS.items():
        output = out_dir / f"{tag}.pdb"
        cmd.reinitialize()
        cmd.load(str(wt_path), "wt")
        selection = f"chain {chain} and resi {resi} and name CA"
        model = cmd.get_model(selection)
        if len(model.atom) != 1 or model.atom[0].resn.upper() != expected:
            raise RuntimeError(f"PyMOL selection/identity check failed for {tag}: {selection}")
        cmd.wizard("mutagenesis")
        wizard = cmd.get_wizard()
        wizard.set_mode(target)
        wizard.do_select(f"chain {chain} and resi {resi}")
        wizard.apply()
        cmd.set_wizard()
        cmd.save(str(output), "wt", state=1)
        cmd.delete("all")

        observed = residue_at(output, chain, resi)
        if observed != target:
            raise RuntimeError(f"Saved {tag} has residue {observed!r}; expected {target}.")
        if output.stat().st_size < 1000:
            raise RuntimeError(f"Output structure is unexpectedly small: {output}")
        rows.append({"variant": tag, "chain": chain, "residue": resi,
                     "expected_residue": target, "observed_residue": observed})
        print(f"[OK] {tag}: {chain}:{expected}{resi}{target} -> {output}")

    qc_path = out_dir / "initial_models_qc.csv"
    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"[DONE] All mutants were built independently from WT. QC: {qc_path}")


main()
