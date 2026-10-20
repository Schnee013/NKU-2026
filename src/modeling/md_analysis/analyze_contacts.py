#!/usr/bin/env python
"""Compute initial receptor-partner interface contact persistence.

Initial residue contacts are configured receptor-interface residues vs all
configured partner residues with minimum heavy-atom distance <= cutoff in
frame 0. Persistence is the fraction of these initial residue pairs retained
over time.
"""
from pathlib import Path
import argparse
import sys


import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

from .common import add_system_argument, analysis_context, ensure_dir, time_ns, replicate_plot_path, in_time_window
from .plotting import lineplot_by_rep

parser = argparse.ArgumentParser()
add_system_argument(parser)
parser.add_argument('--stride', type=int, default=50, help='Frame stride; 50 samples every 500 ps from a 10 ps trajectory')
parser.add_argument('--cutoff', type=float, default=4.0, help='Heavy-atom residue contact cutoff in Angstrom')
args = parser.parse_args()
if args.stride <= 0:
    parser.error("--stride must be a positive integer")
ctx = analysis_context(args.system)

out = ensure_dir(ctx.analysis_base / 'contacts')
summary_rows = []
contact_rows = []

def heavy(sel):
    return sel.select_atoms('not name H* and not type H')

for rep in ctx.replicates():
    print(f'[contacts] {rep.rep} stride={args.stride} cutoff={args.cutoff} A', flush=True)
    u = mda.Universe(str(rep.topology), str(rep.trajectory))
    receptor = heavy(u.select_atoms(ctx.receptor_interface_selection))
    partner = heavy(u.select_atoms(ctx.partner_selection))
    if len(receptor) == 0 or len(partner) == 0:
        raise RuntimeError(
            f'empty selection for {ctx.system} {rep.rep}: '
            f'{ctx.receptor_name}_interface={len(receptor)} '
            f'{ctx.interaction_partner_role}={len(partner)}'
        )

    u.trajectory[0]
    d0 = distance_array(receptor.positions, partner.positions)
    pairs = np.argwhere(d0 <= args.cutoff)
    initial_pairs = set()
    for i, j in pairs:
        initial_pairs.add((receptor[i].resid, receptor[i].resname, partner[j].resid, partner[j].resname))
    if not initial_pairs:
        raise RuntimeError(f'no initial contacts for {rep.rep}')

    # Pre-group atoms by residue for repeated minimum-distance checks.
    receptor_res = {}
    partner_res = {}
    for atom in receptor:
        receptor_res.setdefault((atom.resid, atom.resname), []).append(atom.ix)
    for atom in partner:
        partner_res.setdefault((atom.resid, atom.resname), []).append(atom.ix)

    initial_pairs = sorted(initial_pairs)
    per_pair_counts = {p: 0 for p in initial_pairs}
    n_frames = 0
    for ts in u.trajectory[::args.stride]:
        count_for_occupancy = in_time_window(ts, ctx.window('occupancy_window_ns'))
        retained = 0
        for p in initial_pairs:
            receptor_resid, receptor_resname, partner_resid, partner_resname = p
            receptor_idx = receptor_res[(receptor_resid, receptor_resname)]
            partner_idx = partner_res[(partner_resid, partner_resname)]
            min_dist = distance_array(u.atoms[receptor_idx].positions, u.atoms[partner_idx].positions).min()
            if min_dist <= args.cutoff:
                retained += 1
                if count_for_occupancy:
                    per_pair_counts[p] += 1
        frac = retained / len(initial_pairs)
        summary_rows.append({
            'rep': rep.rep,
            'time_ns': time_ns(ts),
            'initial_contacts': len(initial_pairs),
            'retained_contacts': retained,
            'retained_fraction': frac,
            'cutoff_A': args.cutoff,
            'stride': args.stride,
        })
        if count_for_occupancy:
            n_frames += 1

    if n_frames == 0:
        raise ValueError('contact occupancy window has no sampled frames')

    for p, count in per_pair_counts.items():
        receptor_resid, receptor_resname, partner_resid, partner_resname = p
        row = {
            'rep': rep.rep,
            'receptor_name': ctx.receptor_name,
            'receptor_resid': receptor_resid,
            'receptor_resname': receptor_resname,
            'partner_role': ctx.interaction_partner_role,
            'partner_name': ctx.partner_name,
            'partner_resid': partner_resid,
            'partner_resname': partner_resname,
            'occupancy': count / n_frames,
            'n_frames': n_frames,
            'cutoff_A': args.cutoff,
            'stride': args.stride,
            'window_ns': str(ctx.window('occupancy_window_ns') or 'full'),
        }
        contact_rows.append(row)

summary = pd.DataFrame(summary_rows)
summary.to_csv(out / 'contact_persistence.csv', index=False)
contacts = pd.DataFrame(contact_rows).sort_values(['rep', 'occupancy'], ascending=[True, False])
contacts.to_csv(out / 'initial_contact_occupancy.csv', index=False)
lineplot_by_rep(summary, 'time_ns', 'retained_fraction', replicate_plot_path(out, 'contact_persistence'), f'{ctx.label} initial contact persistence', 'Retained initial-contact fraction')

for rep, sub in summary.groupby('rep'):
    tail = ctx.summary(sub)
    print(
        f'{rep}: initial residue contacts={int(sub.initial_contacts.iloc[0])}; '
        f'{ctx.summary_label} mean retained={tail.retained_fraction.mean():.3f}, '
        f'min={tail.retained_fraction.min():.3f}'
    )
