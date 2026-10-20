"""Protect scientific inputs and publish a complete new output directory."""
from contextlib import contextmanager
from pathlib import Path
import shutil
import tempfile


@contextmanager
def new_output(output, *input_paths):
    output = Path(output).resolve()
    protected = Path(__file__).resolve().parents[3] / 'studies'
    if output == protected or protected in output.parents:
        raise ValueError('output must be outside the read-only study inputs and references')
    for item in input_paths:
        item = Path(item).resolve()
        if output == item or output in item.parents or item in output.parents:
            raise ValueError('output overlaps scientific inputs')
    if output.exists():
        raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.' + output.name + '-', dir=output.parent))
    try:
        yield staging
        if output.exists():
            raise FileExistsError(output)
        staging.rename(output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
