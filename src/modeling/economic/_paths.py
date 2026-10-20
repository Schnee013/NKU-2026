"""Keep generated economic outputs separate from packaged scientific inputs."""
from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[3]
STUDIES_DIR = PROJECT_DIR / "studies"
REPORT_FILES = {
    "diagnostic_report.txt", "diagnostic_report.json",
    "diagnostic_report.html", "diagnostic_report_en.html",
}


def checked_directory(path: str | Path, *, must_be_empty: bool = True) -> Path:
    """Reject reference/input trees and occupied output directories."""
    directory = Path(path).expanduser().resolve()
    if directory == STUDIES_DIR or STUDIES_DIR in directory.parents:
        raise ValueError("Output must be outside the packaged studies/reference/input tree")
    if directory.exists() and not directory.is_dir():
        raise ValueError("Output directory is an existing file")
    if must_be_empty and directory.exists() and any(directory.iterdir()):
        raise ValueError("Output directory must be absent or empty; existing output is not overwritten")
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def report_targets(output_dir: str | Path, names: tuple[str, ...]) -> tuple[Path, ...]:
    """Allow companion report files, never replace any existing artifact."""
    directory = checked_directory(output_dir, must_be_empty=False)
    existing = {p.name for p in directory.iterdir()}
    if not existing.issubset(REPORT_FILES):
        raise ValueError("Output directory contains unrelated files")
    targets = tuple(directory / name for name in names)
    if any(path.exists() or path.is_symlink() for path in targets):
        raise ValueError("Report output already exists; choose a new output directory")
    return targets
