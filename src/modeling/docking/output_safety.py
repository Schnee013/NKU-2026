"""Protect shipped study inputs and reference results from runtime writes."""
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]

def protect_output(path, *, empty_directory=False):
    path=Path(path).expanduser().resolve()
    studies=(REPOSITORY_ROOT/'studies').resolve()
    if path==studies or studies in path.parents:
        raise ValueError('Runtime outputs must be outside the shipped studies directory')
    if empty_directory and path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise FileExistsError(f'Output directory must be new or empty: {path}')
    return path
