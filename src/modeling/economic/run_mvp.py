"""Run the illustrative economic framework from a JSON scenario config."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path
import sys

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from src.modeling.economic.economic_model import evaluate_config
    from src.modeling.economic._paths import PROJECT_DIR, checked_directory
else:
    from .economic_model import evaluate_config
    from ._paths import PROJECT_DIR, checked_directory


def clean_value(value):
    if isinstance(value, float):
        return round(value, 6)
    return value


def clean_rows(rows: list[dict]) -> list[dict]:
    return [{key: clean_value(value) for key, value in row.items()} for row in rows]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run deterministic illustrative economic framework scenarios.")
    parser.add_argument(
        "--config",
        default=str(PROJECT_DIR / "studies/economic-model/configs/public_mvp_framework.json"),
        help="Scenario JSON config.",
    )
    parser.add_argument(
        "--out",
        default="economic-framework-results/results.csv",
        help="CSV output path.",
    )
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict]) -> None:
    checked_directory(path.parent)
    if not rows:
        raise ValueError("no rows to write")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", newline="", encoding="utf-8", dir=path.parent, prefix=".economic-csv-", delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main() -> None:
    args = parse_args()
    config_path = Path(args.config)
    out_path = Path(args.out)

    config = json.loads(config_path.read_text(encoding="utf-8"))
    rows = clean_rows(evaluate_config(config))
    write_csv(out_path, rows)

    print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
