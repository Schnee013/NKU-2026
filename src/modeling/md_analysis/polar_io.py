"""Resolve polar-analysis inputs and keep outputs outside scientific inputs."""
from pathlib import Path
import hashlib
import json
from .common import REPO_ROOT

PATH_KEYS = {'topology', 'reference', 'trajectory', 'alignment_reference',
             'polar_screen', 'screen_path', 'prior_polar_output', 'output',
             'output_root', 'shared_output', 'chemistry_contract', 'sdf',
             'representatives_root'}


def load_config(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    cfg = json.loads(data)
    cases = cfg['cases']
    ids = [case['id'] for case in cases]
    if not cases or not all(ids) or len(set(ids)) != len(ids):
        raise ValueError('Cases must have unique nonempty identities')

    def resolve(value):
        if isinstance(value, list):
            return [resolve(item) for item in value]
        if isinstance(value, dict):
            result = {}
            for key, item in value.items():
                if key in PATH_KEYS:
                    p = Path(item)
                    result[key] = str((p if p.is_absolute() else REPO_ROOT / p).resolve())
                else:
                    result[key] = resolve(item)
            return result
        return value

    cfg = resolve(cfg)
    inputs = [Path(c[key]) for c in cfg['cases'] for key in ('topology', 'reference', 'trajectory')]
    inputs += [Path(c['prior_polar_output']) for c in cfg['cases'] if 'prior_polar_output' in c]
    for key in ('alignment_reference',):
        if key in cfg:inputs.append(Path(cfg[key]))
    compatibility = cfg.get('compatibility', {})
    if 'representatives_root' in compatibility:
        inputs.append(Path(compatibility['representatives_root']))
    outputs = [Path(c['output']) for c in cfg['cases']]
    outputs += [Path(cfg[k]) for k in ('output_root', 'shared_output') if k in cfg]
    if len({c['output'] for c in cfg['cases']}) != len(cases):
        raise ValueError('Each case needs a separate output directory')
    for out in outputs:
        for name in ('studies', 'data', 'src', 'tests', 'environments'):
            protected = (REPO_ROOT / name).resolve()
            if out == protected or protected in out.parents:
                raise ValueError('Output overlaps protected repository data')
        for source in inputs:
            if out == source or out in source.parents or source.parent in out.parents:
                raise ValueError('Output overlaps a scientific input')
    cfg['_path'] = str(path)
    cfg['_config_sha256'] = hashlib.sha256(data).hexdigest()
    return cfg
