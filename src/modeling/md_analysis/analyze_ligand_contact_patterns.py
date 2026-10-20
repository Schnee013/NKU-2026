"""Targeted and all-residue polar-contact analysis of supplied trajectories."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from .polar_io import load_config

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for data in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def table(path, rows):
    import pandas as pd
    pd.DataFrame(rows).to_csv(path, index=False)

def runs(mask, times):
    """Continuous sampled spans; no gap filling or duration extrapolation."""
    edges = np.diff(np.r_[False, mask, False].astype(int))
    starts, ends = np.where(edges == 1)[0], np.where(edges == -1)[0] - 1
    return [(float(times[a]), float(times[b]), float(times[b] - times[a]))
            for a, b in zip(starts, ends)]

def frequent_itemsets(masks, minimum_count):
    """Yield every nonempty conjunction reaching minimum_count (boolean AND)."""
    eligible = np.flatnonzero(masks.sum(axis=0) >= minimum_count).tolist()

    def descend(prefix, current, remaining):
        for pos, feature in enumerate(remaining):
            joint = current & masks[:, feature]
            if int(joint.sum()) >= minimum_count:
                itemset = prefix + (feature,)
                yield itemset, joint
                yield from descend(itemset, joint, remaining[pos + 1:])

    yield from descend((), np.ones(len(masks), dtype=bool), eligible)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--case')
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--polar',action='store_true')
    modes.add_argument('--all-polar',action='store_true')
    modes.add_argument('--compare-all-polar',action='store_true')
    args=parser.parse_args();cfg=load_config(args.config)
    if args.compare_all_polar:
        if args.case:parser.error('--compare-all-polar uses the configured target and control sets')
        from .polar_contact_screen import compare_all_polar
        compare_all_polar(cfg)
        return
    cases={c['id']:c for c in cfg['cases']}
    if args.case not in cases:parser.error('Choose a configured --case')
    if args.polar:
        from .polar_contacts import run_polar_case
        run_polar_case(cfg,cases[args.case])
    else:
        from .polar_contact_screen import scan_all_polar
        scan_all_polar(cfg,cases[args.case])

if __name__=='__main__':main()
