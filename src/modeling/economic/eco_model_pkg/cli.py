"""Generate the fixed Spain protected-tomato model-scenario report."""
from __future__ import annotations

import argparse
import os
import shutil
import tempfile
from pathlib import Path

from .._paths import PROJECT_DIR, checked_directory
from .break_even import break_even_analysis
from .config import Model, load_density, load_model
from .economics import evaluate
from .reporting import (
    build_diagnostic_report, open_report_in_browser,
    write_diagnostic_report, write_english_html_report,
)

CONFIG_DIR = PROJECT_DIR / "studies" / "economic-model" / "configs"


def print_report(model: Model, density: float | None) -> None:
    print(f"Model scenario: {model.country_en}; {model.crop_en}; {model.nematode_en}")
    print(f"Density: {model.density_unit}; money: {model.currency}/{model.area_unit}")
    for baseline, label in (("no_control", "Testing plus control vs no control"),
                            ("uniform_control", "Density-guided vs uniform control")):
        roots, jumps = break_even_analysis(model, baseline=baseline)
        print(f"{label}: strict zero-benefit densities = {roots}; tier jumps = {jumps}")
    if density is not None:
        result = evaluate(density, model)
        print(f"Input density = {density:g}; model-scenario strategy = {result.treatment.name_en}")
        print(f"Program net benefit vs no control = {result.program_net_benefit_vs_no_control:.8f}")
        print(f"Conditional testing value vs uniform control = {result.detection_net_benefit_vs_uniform:.8f}")
        print("Test cost is a team pricing assumption; input density is treated as known and accurate.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("density", nargs="?", type=float, help="Override input density (J2/100 cm^3 soil)")
    parser.add_argument("--config", type=Path, default=CONFIG_DIR / "eco_model_example_model.json")
    parser.add_argument("--input", type=Path, default=CONFIG_DIR / "eco_model_example_input.json")
    parser.add_argument("--break-even-only", action="store_true")
    parser.add_argument("--report-dir", type=Path, default=Path.cwd() / "economic-results")
    parser.add_argument("--no-report", action="store_true")
    parser.add_argument("--open", action="store_true", help="Open generated English HTML in a browser")
    args = parser.parse_args()
    try:
        model = load_model(args.config)
        density = None if args.break_even_only else args.density if args.density is not None else load_density(args.input)
        # Evaluate before creating output; failed input cannot leave a report.
        report = None if density is None else build_diagnostic_report(model, evaluate(density, model))
        if report is not None and not args.no_report:
            output = checked_directory(args.report_dir)
            # Publish a complete directory only after every report was rendered.
            # A failed/interrupted render cannot leave a partially valid report
            # at the requested final location.
            stage = Path(tempfile.mkdtemp(prefix=".economic-report-", dir=output.parent))
            try:
                staged_paths = write_diagnostic_report(report, stage)
                staged_english = write_english_html_report(report, stage)
                os.replace(stage, output)
            finally:
                if stage.exists():
                    shutil.rmtree(stage)
            paths = tuple(output / p.name for p in staged_paths)
            english_path = output / staged_english.name
            print("Reports: " + ", ".join(str(p) for p in (*paths, english_path)))
            if args.open:
                open_report_in_browser(english_path)
        print_report(model, density)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
