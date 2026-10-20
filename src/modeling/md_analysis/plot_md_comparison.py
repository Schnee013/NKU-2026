"""Plot all configured systems from saved basic and local numerical tables."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import REPO_ROOT, ensure_dir

STUDY = REPO_ROOT/'studies/receptor-engineering/molecular-dynamics'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', type=Path, default=STUDY)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = ensure_dir(args.output)
    cases = json.loads((args.study/'inputs.json').read_text())
    plt.rcParams.update({'font.size': 9, 'svg.fonttype': 'none', 'axes.spines.top': False,
                         'axes.spines.right': False, 'font.family': 'DejaVu Sans'})
    colors = plt.get_cmap('tab10').colors
    labels = [c['label'].replace(' + ', ' / ') for c in cases]
    definitions = [
        ('convergence/rmsd_summary.csv', 'backbone_rmsd_A', 'Receptor backbone RMSD (Å)'),
        ('convergence/rg_summary.csv', 'rg_A', 'Receptor Cα radius of gyration (Å)'),
        ('selection_rmsd/selection_rmsd.csv', 'selection_rmsd_A', 'Ligand heavy-atom RMSD (Å)'),
        ('com_distance/com.csv', 'distance_A', 'Ligand COM to own initial-contact cluster (Å)'),
        ('bsa/bsa.csv', 'bsa_A2', 'Buried surface area (Å²)'),
        ('contacts/contact_persistence.csv', 'retained_fraction', 'Initial-contact retention (fraction)'),
    ]

    def finish(fig, name):
        fig.legend(loc='lower center', ncol=2, bbox_to_anchor=(.5, .015), fontsize=8)
        fig.tight_layout(rect=(0, .18, 1, .95))
        for ext in ('png', 'svg'):
            fig.savefig(out/f'{name}.{ext}', dpi=200, facecolor='white')
        plt.close(fig)

    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    for ax, (filename, column, title) in zip(axes.flat, definitions):
        for i, case in enumerate(cases):
            table = pd.read_csv(args.study/'results/base'/case['dataset_id']/filename)
            if 'cluster' in table:
                table = table[table.cluster == 'own_initial_4A_contacts']
            ax.plot(table.time_ns, table[column], color=colors[i], lw=.9,
                    label=labels[i] if ax is axes[0, 0] else None)
        ax.set(xlabel='Time (ns)', ylabel=title, xlim=(0, 100))
        ax.axvspan(50, 100, color='#e8edf2', zorder=-1)
    fig.suptitle('Eight independently prepared Gpr2 systems: one 100 ns trajectory per system')
    finish(fig, 'eight_case_timeseries')

    fig, ax = plt.subplots(figsize=(12, 5.5))
    for i, case in enumerate(cases):
        table = pd.read_csv(args.study/'results/base'/case['dataset_id']/'rmsf/rmsf.csv')
        ax.plot(table.resid, table.rmsf_A, color=colors[i], lw=1, label=labels[i])
    ax.set(xlabel='Receptor residue', ylabel='Cα RMSF, 50–100 ns (Å)', xlim=(1, 406))
    fig.suptitle('Receptor fluctuations on backbone-fitted coordinates')
    finish(fig, 'eight_case_receptor_rmsf')

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for row, metric in enumerate(('ligand_rmsd', 'lid_rmsd')):
        for col, fit in enumerate(('global', 'core', 'pocket')):
            ax = axes[row, col]
            for i, case in enumerate(cases):
                table = pd.read_csv(args.study/'results/local'/case['dataset_id']/'frame_metrics.csv')
                ax.plot(table.time_ns, table[f'{metric}_{fit}_A'], color=colors[i], lw=.9,
                        label=labels[i] if row == col == 0 else None)
            ax.set(xlabel='Time (ns)', ylabel=('Ligand' if row == 0 else 'Lid')+' RMSD (Å)',
                   title=fit.capitalize()+' fit', xlim=(0, 100))
            ax.axvspan(50, 100, color='#e8edf2', zorder=-1)
    fig.suptitle('Local structural motion: fixed receptor, core and pocket selections')
    finish(fig, 'eight_case_local_fits')


if __name__ == '__main__':
    main()
