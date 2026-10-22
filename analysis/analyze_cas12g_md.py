#!/usr/bin/env python3
"""Generate per-trajectory RMSD, Rg, C-alpha RMSF and RMSD-Rg occupancy/FEL outputs."""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def run(command: list[str], stdin: str | None = None) -> str:
    result = subprocess.run(command, input=stdin, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("Command failed:\n" + " ".join(command) + "\n" + result.stdout[-5000:])
    return result.stdout


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            continue
    if not rows:
        raise RuntimeError(f"No numeric rows found in {path}")
    return np.asarray(rows, dtype=float)


def group_ids(make_ndx_output: str) -> dict[str, int]:
    groups = {}
    for line in make_ndx_output.splitlines():
        match = re.match(r"^\s*(\d+)\s+(.+?)\s*:\s+\d+\s+atoms\s*$", line)
        if match:
            groups[match.group(2).strip()] = int(match.group(1))
    return groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--deffnm", default="md400")
    parser.add_argument("--begin-ns", type=float, default=50.0)
    parser.add_argument("--temp-k", type=float, default=300.0)
    parser.add_argument("--fel-bins", type=int, default=30)
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    outdir = (args.outdir or workdir / "analysis").resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    tpr = workdir / f"{args.deffnm}.tpr"
    xtc = workdir / f"{args.deffnm}.xtc"
    if not tpr.is_file() or not xtc.is_file():
        raise FileNotFoundError(f"Expected production files: {tpr} and {xtc}")
    if args.fel_bins < 5:
        raise ValueError("--fel-bins must be at least 5")

    index_file = outdir / "analysis_index.ndx"
    ndx_log = run([args.gmx, "make_ndx", "-f", str(tpr), "-o", str(index_file)], "q\n")
    groups = group_ids(ndx_log)
    ca_group = groups.get("C-alpha")
    protein_group = groups.get("Protein")
    if ca_group is None or protein_group is None:
        raise RuntimeError("Could not find default GROMACS groups 'C-alpha' and 'Protein'. "
                           "Inspect make_ndx output and adjust the group parser for your GROMACS version.")

    rmsd_xvg = outdir / "rmsd_ca.xvg"
    rg_xvg = outdir / "rg_protein.xvg"
    rmsf_xvg = outdir / "rmsf_ca.xvg"
    run([args.gmx, "rms", "-s", str(tpr), "-f", str(xtc), "-n", str(index_file),
         "-o", str(rmsd_xvg), "-tu", "ns"], f"{ca_group}\n{ca_group}\n")
    run([args.gmx, "gyrate", "-s", str(tpr), "-f", str(xtc), "-n", str(index_file),
         "-o", str(rg_xvg), "-tu", "ns"], f"{protein_group}\n")
    run([args.gmx, "rmsf", "-s", str(tpr), "-f", str(xtc), "-n", str(index_file),
         "-o", str(rmsf_xvg), "-res", "-b", str(args.begin_ns * 1000.0)], f"{ca_group}\n")

    rmsd_data = read_xvg(rmsd_xvg)
    rg_data = read_xvg(rg_xvg)
    rmsf_data = read_xvg(rmsf_xvg)
    time = rmsd_data[:, 0]
    rmsd = rmsd_data[:, 1]
    rg = np.interp(time, rg_data[:, 0], rg_data[:, 1])
    metrics = pd.DataFrame({"time_ns": time, "rmsd_ca_nm": rmsd, "rg_protein_nm": rg})
    metrics.to_csv(outdir / f"metrics_{args.tag}.csv", index=False)

    rmsf = pd.DataFrame({"residue_index": rmsf_data[:, 0].astype(int),
                         "rmsf_ca_nm": rmsf_data[:, 1]})
    rmsf.to_csv(outdir / f"rmsf_{args.tag}.csv", index=False)

    mask = (time >= args.begin_ns) & np.isfinite(rmsd) & np.isfinite(rg)
    if mask.sum() < args.fel_bins:
        raise RuntimeError(f"Only {int(mask.sum())} frames are available at/after {args.begin_ns} ns.")
    hist, xedges, yedges = np.histogram2d(rmsd[mask], rg[mask], bins=args.fel_bins)
    probability = hist / hist.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        fel = -0.008314462618 * args.temp_k * np.log(probability / np.nanmax(probability))
    fel[probability == 0] = np.nan
    xcenters = (xedges[:-1] + xedges[1:]) / 2
    ycenters = (yedges[:-1] + yedges[1:]) / 2
    xx, yy = np.meshgrid(xcenters, ycenters, indexing="ij")
    fel_table = pd.DataFrame({
        "rmsd_ca_bin_nm": xx.ravel(),
        "rg_protein_bin_nm": yy.ravel(),
        "probability": probability.ravel(),
        "relative_fel_kj_mol": fel.ravel(),
    })
    fel_table.to_csv(outdir / f"fel_rmsd_rg_{args.tag}.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    axes[0].plot(time, rmsd, linewidth=0.8)
    axes[0].set(xlabel="Time (ns)", ylabel="C-alpha RMSD (nm)", title=f"{args.tag}: RMSD")
    axes[1].plot(time, rg, linewidth=0.8)
    axes[1].set(xlabel="Time (ns)", ylabel="Protein Rg (nm)", title="Radius of gyration")
    mesh = axes[2].pcolormesh(xedges, yedges, fel.T, shading="auto", cmap="viridis")
    axes[2].set(xlabel="C-alpha RMSD (nm)", ylabel="Protein Rg (nm)",
                title=f"Relative occupancy FEL (from {args.begin_ns:g} ns)")
    fig.colorbar(mesh, ax=axes[2], label="Relative FEL (kJ/mol)")
    fig.suptitle(f"Cas12g MD summary: {args.tag}")
    fig.savefig(outdir / f"summary_{args.tag}.png", dpi=220)
    plt.close(fig)

    summary = [
        f"tag={args.tag}",
        f"workdir={workdir}",
        f"frames={len(time)}",
        f"analysis_begin_ns={args.begin_ns}",
        f"rmsd_mean_nm={np.nanmean(rmsd[mask]):.4f}",
        f"rmsd_sd_nm={np.nanstd(rmsd[mask], ddof=1):.4f}",
        f"rg_mean_nm={np.nanmean(rg[mask]):.4f}",
        f"rg_sd_nm={np.nanstd(rg[mask], ddof=1):.4f}",
        "interpretation=FEL is a descriptive occupancy-derived surface, not proof of convergence.",
    ]
    (outdir / f"summary_{args.tag}.txt").write_text("\n".join(summary) + "\n")
    print(f"[DONE] Analysis outputs written to {outdir}")


if __name__ == "__main__":
    main()
