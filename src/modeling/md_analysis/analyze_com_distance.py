#!/usr/bin/env python
from pathlib import Path
import argparse
import sys, numpy as np
import MDAnalysis as mda
from .common import add_system_argument, analysis_context, ensure_dir, time_ns, linear_slope, replicate_plot_path, save_csv
from .plotting import lineplot_by_rep

parser = argparse.ArgumentParser()
add_system_argument(parser)
args = parser.parse_args()
ctx = analysis_context(args.system)

reps = ctx.replicates()
for rep in reps:
    universe = mda.Universe(str(rep.topology), str(rep.trajectory))
    target = universe.select_atoms(ctx.partner_com_selection)
    if len(target) == 0:
        raise RuntimeError(
            f"empty COM target selection for {ctx.system} {rep.rep}: "
            f"{ctx.partner_com_selection}"
        )
    for name, selection in ctx.receptor_interface_clusters.items():
        atoms = universe.select_atoms(selection)
        if len(atoms) == 0:
            raise RuntimeError(
                f"empty COM cluster selection for {ctx.system} {rep.rep} "
                f"cluster={name}: {selection}"
            )

out = ensure_dir(ctx.analysis_base / 'com_distance')
rows=[]
for rep in reps:
    u=mda.Universe(str(rep.topology), str(rep.trajectory))
    partner_atoms=u.select_atoms(ctx.partner_com_selection)
    clusters={name: u.select_atoms(spec) for name, spec in ctx.receptor_interface_clusters.items()}
    for ts in u.trajectory:
        c2=partner_atoms.center_of_mass()
        for name, sel in clusters.items():
            rows.append({'rep': rep.rep, 'time_ns': time_ns(ts), 'cluster': name, 'distance_A': float(np.linalg.norm(sel.center_of_mass()-c2))})

df=save_csv(out/'com.csv', rows)
for cluster, sub in df.groupby('cluster'):
    lineplot_by_rep(
        sub,
        'time_ns',
        'distance_A',
        replicate_plot_path(out, cluster),
        f'{ctx.label} COM distance {cluster} to {ctx.com_target_label}',
        'COM distance (A)',
    )
context_note = ctx.com_distance_context_note
for (rep, cluster), sub in df.groupby(['rep','cluster']):
    tail=ctx.summary(sub)
    slope=linear_slope(tail.time_ns.to_numpy(), tail.distance_A.to_numpy())
    print(
        f"{rep} {cluster}: {ctx.summary_label} slope={slope:.4g} A/ns; "
        f"mean={tail.distance_A.mean():.2f} A"
        + (f"; {context_note}" if context_note else "")
    )
