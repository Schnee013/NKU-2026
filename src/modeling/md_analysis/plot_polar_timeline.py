"""Plot saved 1-ns occupancies; no trajectory reanalysis or smoothing."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

import argparse
from .common import ensure_dir

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--table',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();ensure_dir(args.output)
    data = pd.read_csv(args.table)
    data = data[data.case_id == 'v368r_ascr3_arg189_start']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(10.1, 3.9))
    fig.subplots_adjust(left=.09, bottom=.18, right=.98, top=.73)
    fig.text(.09, .945, 'Asp89 anchoring is frequent for one stage, then largely lost',
             fontsize=16, weight='bold', color='#263945')
    fig.text(.09, .865, 'V368R–ascr3, Arg189-directed start', fontsize=11.5, color='#52636f')
    ax.axvspan(20, 50, color='#dff0ec', zorder=0)
    ax.axvspan(50, 100, color='#f1f4f6', zorder=0)
    ax.plot(data.time_ns, data.asp89_hbond_occupancy, color='#087d73', lw=2.15,
            label='Asp89–sugar hydrogen bond', zorder=3)
    ax.plot(data.time_ns, data.arg345_salt_occupancy, color='#295fb8', lw=1.8,
            label='Arg345–carboxylate salt bridge', zorder=2)
    ax.set(xlim=(0, 100), ylim=(-.035, 1.08), xlabel='Time (ns)', ylabel='Occupancy in 1 ns bins')
    ax.set_xticks(range(0, 101, 10))
    ax.set_yticks([0, .5, 1])
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    ax.grid(axis='y', color='#dce2e6', linewidth=.7)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.05), frameon=False, ncol=2, fontsize=10.3)
    for ext in ['png', 'svg', 'pdf']:
        fig.savefig(args.output / ('figure1_asp89_timeline.' + ext), dpi=260)
    print('Saved Asp89 / Arg345 time-series figure.')

if __name__ == '__main__':
    main()
