"""Plot a saved six-state model run without recomputing its trajectory."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',type=Path)
    parser.add_argument('--output',type=Path,required=True,
                        help='New PNG file outside packaged input/reference directories')
    args=parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    destination = args.output.resolve()
    protected = root/'studies/signal-model/pathway'
    if destination == protected or destination.is_relative_to(protected):
        parser.error('--output must be outside packaged input and reference directories')
    if destination == args.csv.resolve() or destination.exists():
        parser.error('--output already exists or equals input; choose a new PNG file')
    data=np.atleast_1d(np.genfromtxt(args.csv,delimiter=',',names=True))
    needed={'ligand_normalized','time_s','MAPK_fraction','color_au'}
    if not data.size or not needed.issubset(data.dtype.names or ()):
        parser.error('CSV must contain ligand_normalized, time_s, MAPK_fraction and color_au')
    if any(not np.isfinite(data[name]).all() for name in needed):
        parser.error('Plot columns must be finite')
    fig,axes=plt.subplots(1,3,figsize=(12,3.7))
    for dose in np.unique(data['ligand_normalized']):
        row=data[data['ligand_normalized']==dose]
        early=row['time_s']<=900
        axes[0].plot(row['time_s'][early]/60,row['MAPK_fraction'][early],label=f'L = {dose:g}')
        axes[1].plot(row['time_s']/3600,row['MAPK_fraction'])
        axes[2].plot(row['time_s']/3600,row['color_au'])
    for ax in axes:
        ax.grid(alpha=.2)
        ax.set_ylim(bottom=0)
    axes[0].set(xlabel='Time (min)',ylabel='Active Fus3 / total',title='Early core response')
    axes[1].set(xlabel='Time (h)',ylabel='Active Fus3 / total',title='Sustained-input scenario')
    axes[2].set(xlabel='Time (h)',ylabel='Color (arbitrary units)',title='Effective color output')
    axes[0].legend(fontsize=8)
    fig.suptitle('Six-state ODE prototype | normalized ligand | uncalibrated endpoints',fontsize=11)
    fig.tight_layout()
    out=args.output
    # Exclusive creation also prevents an overwrite if another process races us.
    with out.open('xb') as handle:
        fig.savefig(handle,dpi=180,format='png')
    plt.close(fig)
    print(out)


if __name__=='__main__':
    main()
