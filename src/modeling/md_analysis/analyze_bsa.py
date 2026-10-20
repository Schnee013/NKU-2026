#!/usr/bin/env python
"""Compute configured receptor-partner BSA from retained-complex trajectories.

Uses MDTraj Shrake-Rupley SASA on stride-sampled frames:
BSA = SASA(partner) + SASA(receptor) - SASA(complex)
Units are converted from nm^2 to A^2.
"""
from pathlib import Path
import argparse
import sys


import mdtraj as md
import MDAnalysis as mda
import pandas as pd

from .common import add_system_argument, analysis_context, ensure_dir, replicate_plot_path
from .plotting import lineplot_by_rep

parser = argparse.ArgumentParser()
add_system_argument(parser)
parser.add_argument('--stride', type=int, default=50, help='Frame stride; 50 samples every 500 ps from a 10 ps trajectory')
parser.add_argument('--n-sphere-points', type=int, default=240, help='Shrake-Rupley sphere points; lower is faster, relative BSA trend is the target')
args = parser.parse_args()
if args.stride <= 0:
    parser.error("--stride must be a positive integer")
ctx = analysis_context(args.system)

out = ensure_dir(ctx.analysis_base / 'bsa')
rows = []
for rep in ctx.replicates():
    print(f'[bsa] loading {rep.rep} stride={args.stride}', flush=True)
    u = mda.Universe(str(rep.topology))
    partner_idx = u.select_atoms(ctx.partner_selection).indices
    receptor_idx = u.select_atoms(ctx.receptor_selection).indices
    if len(partner_idx) == 0 or len(receptor_idx) == 0:
        raise RuntimeError(
            f'empty BSA selection for {ctx.system} {rep.rep}: '
            f'{ctx.partner_name}={len(partner_idx)} {ctx.receptor_name}={len(receptor_idx)}'
        )

    traj = md.load_xtc(str(rep.trajectory), top=str(rep.topology), stride=args.stride)
    time_ns = traj.time / 1000.0

    complex_sasa = md.shrake_rupley(traj, mode='atom', n_sphere_points=args.n_sphere_points).sum(axis=1)
    partner_sasa = md.shrake_rupley(traj.atom_slice(partner_idx), mode='atom', n_sphere_points=args.n_sphere_points).sum(axis=1)
    receptor_sasa = md.shrake_rupley(traj.atom_slice(receptor_idx), mode='atom', n_sphere_points=args.n_sphere_points).sum(axis=1)
    bsa_a2 = (partner_sasa + receptor_sasa - complex_sasa) * 100.0

    initial = float(bsa_a2[0]) if len(bsa_a2) else float('nan')
    for t, val in zip(time_ns, bsa_a2):
        rows.append({
            'rep': rep.rep,
            'time_ns': float(t),
            'bsa_A2': float(val),
            'bsa_fraction_of_initial': float(val / initial) if initial else float('nan'),
            'stride': args.stride,
            'n_sphere_points': args.n_sphere_points,
        })

df = pd.DataFrame(rows)
df.to_csv(out / 'bsa.csv', index=False)
lineplot_by_rep(df, 'time_ns', 'bsa_A2', replicate_plot_path(out, 'bsa'), f'{ctx.label} {ctx.receptor_name}-{ctx.partner_name} BSA', 'BSA (A^2)')
lineplot_by_rep(df, 'time_ns', 'bsa_fraction_of_initial', replicate_plot_path(out, 'bsa_fraction'), f'{ctx.label} BSA fraction of initial', 'BSA / initial')

for rep, sub in df.groupby('rep'):
    tail = ctx.summary(sub)
    print(
        f'{rep}: {ctx.summary_label} '
        f'mean fraction={tail.bsa_fraction_of_initial.mean():.3f}, '
        f'min fraction={tail.bsa_fraction_of_initial.min():.3f}'
    )
