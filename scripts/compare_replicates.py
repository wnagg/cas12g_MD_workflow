#!/usr/bin/env python3
"""Summarize trajectory metrics across independent Cas12g replicas."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--roots", nargs="+", type=Path, required=True,
                        help="One or more result roots, each containing runs/repN/VARIANT/analysis.")
    parser.add_argument("--replicate-ids", nargs="+", type=int, required=True)
    parser.add_argument("--variants", nargs="+", default=["WT", "V623Q", "V623S"])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--begin-ns", type=float, default=50.0)
    parser.add_argument("--target-ns", type=float, default=400.0)
    parser.add_argument("--final-window-ns", type=float, default=50.0)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    records = []
    time_series = {}
    for variant in args.variants:
        frames = []
        for root in args.roots:
            for rep in args.replicate_ids:
                path = root / "runs" / f"rep{rep}" / variant / "analysis" / f"metrics_{variant}.csv"
                if not path.is_file():
                    print(f"[WARNING] missing metrics: {path}")
                    continue
                df = pd.read_csv(path)
                df["replicate_id"] = rep
                df["variant"] = variant
                df["source_root"] = str(root)
                frames.append(df)
                end_time = float(df["time_ns"].max())
                final = df[df["time_ns"] >= max(args.begin_ns, end_time - args.final_window_ns)]
                complete = end_time >= args.target_ns - 0.5
                records.append({
                    "source_root": str(root),
                    "replicate_id": rep,
                    "variant": variant,
                    "last_time_ns": end_time,
                    "target_reached": complete,
                    "final_window_start_ns": float(final["time_ns"].min()) if len(final) else np.nan,
                    "rmsd_mean_nm": float(final["rmsd_ca_nm"].mean()) if len(final) else np.nan,
                    "rmsd_sd_nm": float(final["rmsd_ca_nm"].std()) if len(final) > 1 else np.nan,
                    "rg_mean_nm": float(final["rg_protein_nm"].mean()) if len(final) else np.nan,
                    "rg_sd_nm": float(final["rg_protein_nm"].std()) if len(final) > 1 else np.nan,
                })
        if not frames:
            continue
        combined = pd.concat(frames, ignore_index=True)
        combined = combined[combined["time_ns"] >= args.begin_ns]
        summary = combined.groupby("time_ns", as_index=False).agg(
            n_replicas=("replicate_id", "nunique"),
            rmsd_mean_nm=("rmsd_ca_nm", "mean"),
            rmsd_sd_nm=("rmsd_ca_nm", "std"),
            rg_mean_nm=("rg_protein_nm", "mean"),
            rg_sd_nm=("rg_protein_nm", "std"),
        )
        summary.to_csv(args.outdir / f"metrics_{variant}_mean_sd.csv", index=False)
        time_series[variant] = summary

    if records:
        pd.DataFrame(records).to_csv(args.outdir / "replicate_summary.csv", index=False)
    else:
        raise RuntimeError("No per-replica metrics were found. Check result-root paths and analysis output.")

    for metric, mean_col, sd_col, ylabel, filename in [
        ("RMSD", "rmsd_mean_nm", "rmsd_sd_nm", "C-alpha RMSD (nm)", "comparison_RMSD.png"),
        ("Rg", "rg_mean_nm", "rg_sd_nm", "Protein Rg (nm)", "comparison_Rg.png"),
    ]:
        fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True)
        for variant, df in time_series.items():
            x = df["time_ns"].to_numpy()
            mean = df[mean_col].to_numpy()
            sd = df[sd_col].fillna(0).to_numpy()
            ax.plot(x, mean, linewidth=1.4, label=variant)
            ax.fill_between(x, mean - sd, mean + sd, alpha=0.18)
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Cas12g replicate comparison: {metric}")
        ax.legend(frameon=False)
        fig.savefig(args.outdir / filename, dpi=220)
        plt.close(fig)
    print(f"[DONE] Cross-replica summaries written to {args.outdir}")


if __name__ == "__main__":
    main()
