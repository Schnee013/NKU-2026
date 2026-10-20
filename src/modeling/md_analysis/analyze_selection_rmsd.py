#!/usr/bin/env python
"""Compute backbone-fitted RMSD for a config-selected molecular readout.
"""
from pathlib import Path
import argparse
import sys

import numpy as np

import MDAnalysis as mda
from .common import add_system_argument, analysis_context, ensure_dir, equal_weight_kabsch_coordinates, time_ns, replicate_plot_path, save_csv
from .plotting import lineplot_by_rep

parser = argparse.ArgumentParser()
add_system_argument(parser)
args = parser.parse_args()
ctx = analysis_context(args.system)

metric_dir = ctx.config.get("selection_rmsd_dir", "selection_rmsd")
metric_label = ctx.config.get("selection_rmsd_label", "selection RMSD")
metric_selection_label = ctx.config.get("selection_rmsd_selection_label", "configured selection")

reps = ctx.replicates()
for rep in reps:
    mobile = mda.Universe(str(rep.topology), str(rep.trajectory))
    reference = mda.Universe(str(rep.topology), str(rep.trajectory))
    reference.trajectory[0]
    mobile_fit = mobile.select_atoms(ctx.backbone_selection)
    reference_fit = reference.select_atoms(ctx.backbone_selection)
    mobile_selection = mobile.select_atoms(ctx.selection_rmsd_selection)
    reference_selection = reference.select_atoms(ctx.selection_rmsd_selection)
    for label, mobile_atoms, reference_atoms in [
        ("fit", mobile_fit, reference_fit),
        ("readout", mobile_selection, reference_selection),
    ]:
        if len(mobile_atoms) == 0 or len(reference_atoms) == 0:
            raise RuntimeError(
                f"empty selection-RMSD {label} selection for {ctx.system} {rep.rep}: "
                f"mobile={len(mobile_atoms)} reference={len(reference_atoms)}"
            )
        if len(mobile_atoms) != len(reference_atoms) or not np.array_equal(
            mobile_atoms.indices, reference_atoms.indices
        ):
            raise RuntimeError(
                f"selection-RMSD {label} atom count/order mismatch for "
                f"{ctx.system} {rep.rep}"
            )

out = ensure_dir(ctx.analysis_base / metric_dir)
rows=[]
for rep in reps:
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    ref = mda.Universe(str(rep.topology), str(rep.trajectory))
    ref.trajectory[0]
    mob_fit = u.select_atoms(ctx.backbone_selection); ref_fit = ref.select_atoms(ctx.backbone_selection)
    mob_sel = u.select_atoms(ctx.selection_rmsd_selection); ref_sel = ref.select_atoms(ctx.selection_rmsd_selection)
    ref_fit_pos = ref_fit.positions.copy()
    ref_sel_pos = ref_sel.positions.copy()
    for ts in u.trajectory:
        sel, _ = equal_weight_kabsch_coordinates(
            mob_fit.positions,
            ref_fit_pos,
            mob_sel.positions,
        )
        diff = sel - ref_sel_pos
        val = float(np.sqrt((diff * diff).sum() / len(sel)))
        rows.append({'rep': rep.rep, 'time_ns': time_ns(ts), 'selection_rmsd_A': val})

df = save_csv(out / 'selection_rmsd.csv', rows)
lineplot_by_rep(
    df,
    'time_ns',
    'selection_rmsd_A',
    replicate_plot_path(out, 'selection_rmsd'),
    f'{ctx.label} {metric_label}',
    f'{metric_label} (A)',
)
for rep, sub in df.groupby('rep'):
    tail = ctx.summary(sub)
    print(
        f"{rep}: {metric_selection_label} {ctx.summary_label} "
        f"mean={tail.selection_rmsd_A.mean():.2f} A, "
        f"max={tail.selection_rmsd_A.max():.2f} A"
    )
