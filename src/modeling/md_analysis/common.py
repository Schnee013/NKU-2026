"""Common helpers for config-driven GPCR-complex MD analysis.

Each system defines its molecular roles, atom selections and analysis windows.
The default system is the WT-ascr18 reference specified in systems.json.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align

REPO_ROOT = Path(__file__).resolve().parents[3]
CONFIG_PATH = REPO_ROOT / "studies/receptor-engineering/molecular-dynamics/configs/systems.json"

@dataclass(frozen=True)
class Replicate:
    rep: str
    topology: Path
    trajectory: Path


@dataclass(frozen=True)
class AnalysisContext:
    system: str
    label: str
    config: dict
    base: Path
    payload_base: Path

    @property
    def receptor_name(self) -> str:
        return self.config["receptor"]["name"]

    @property
    def partner_name(self) -> str:
        return self.config["partner"]["name"]

    @property
    def receptor_selection(self) -> str:
        return self.config["receptor"]["selection"]

    @property
    def partner_selection(self) -> str:
        return self.config["partner"]["selection"]

    @property
    def protein_selection(self) -> str:
        return f"({self.partner_selection}) or ({self.receptor_selection})"

    @property
    def backbone_selection(self) -> str:
        return self.config.get("backbone_selection") or f"({self.protein_selection}) and backbone"

    @property
    def ca_selection(self) -> str:
        return self.config.get("ca_selection") or f"({self.protein_selection}) and name CA"

    def window(self, key: str) -> tuple[float, float] | None:
        value = self.config.get(key)
        if value is None:
            return None
        if len(value) != 2:
            raise ValueError(f"{key} requires [start_ns, end_ns]")
        start, end = map(float, value)
        if not np.isfinite([start, end]).all() or not 0 <= start < end:
            raise ValueError(f"invalid {key}: {value}")
        return start, end

    def summary(self, table: pd.DataFrame) -> pd.DataFrame:
        window = self.window("summary_window_ns")
        if window is None:
            return table[table.time_ns >= 100] if table.time_ns.max() >= 100 else table
        selected = table[table.time_ns.between(*window)]
        if len(selected) < 2:
            raise ValueError(f"summary window {window} contains fewer than two samples")
        return selected

    @property
    def summary_label(self) -> str:
        window = self.window("summary_window_ns")
        return f"{window[0]:g}-{window[1]:g} ns" if window else "after100"

    @property
    def selection_rmsd_selection(self) -> str:
        selection = self.config["selection_rmsd_selection"]
        if not selection:
            raise ValueError("selection_rmsd_selection must be explicitly configured")
        return selection

    @property
    def interaction_partner_role(self) -> str:
        """Return the scientific role of the interaction partner."""
        return self.config.get("interaction_partner_role", "partner")

    @property
    def com_target_label(self) -> str:
        return self.config.get("com_target_label") or self.partner_name

    @property
    def com_distance_context_note(self) -> str | None:
        return self.config.get("com_distance_context_note")

    @property
    def partner_com_selection(self) -> str:
        return self.config["partner"]["com_selection"]

    @property
    def receptor_interface_selection(self) -> str:
        return self.config["receptor"]["interface_selection"]

    @property
    def receptor_interface_clusters(self) -> dict[str, str]:
        clusters = self.config["receptor"].get("interface_clusters") or {}
        return clusters or {"receptor_interface": self.receptor_interface_selection}

    @property
    def plot_segids(self) -> list[str]:
        return list(self.config.get("plot_segids") or ["PROA", "PROB"])

    @property
    def analysis_base(self) -> Path:
        return self.base / "analysis"

    def replicates(self) -> list[Replicate]:
        out = []
        for i in configured_replicate_numbers():
            d = self.payload_base / "processed" / f"rep{i}"
            out.append(Replicate(f"rep{i}", d / "protein_ref.pdb", d / "prod_protein.xtc"))
        return out


def configured_replicate_numbers() -> list[int]:
    raw = os.environ.get("MD_REPS")
    if not raw:
        return [1]
    reps: list[int] = []
    for item in raw.replace(",", " ").split():
        token = item.strip().lower()
        if token.startswith("rep"):
            token = token[3:]
        try:
            rep = int(token)
        except ValueError as exc:
            raise ValueError(f"invalid MD_REPS token: {item!r}") from exc
        if rep <= 0:
            raise ValueError(f"MD_REPS entries must be positive: {item!r}")
        reps.append(rep)
    if not reps:
        raise ValueError("MD_REPS is set but no replicate numbers were parsed")
    duplicates = sorted({rep for rep in reps if reps.count(rep) > 1})
    if duplicates:
        raise ValueError(
            "MD_REPS contains duplicate replicate numbers: "
            + ", ".join(map(str, duplicates))
        )
    return reps


def replicate_output_tag() -> str:
    """Return a filename tag for the configured replicate numbers."""
    reps = configured_replicate_numbers()
    if len(reps) == 1:
        return f"rep{reps[0]}"
    return f"{len(reps)}rep"


def replicate_plot_path(directory: Path, stem: str) -> Path:
    return directory / f"{stem}_{replicate_output_tag()}.png"


def load_config(path: Path | None = None) -> dict:
    path = path or Path(os.environ.get("MD_CONFIG", CONFIG_PATH))
    return json.loads(path.read_text(encoding="utf-8"))


def available_systems() -> list[str]:
    return sorted(load_config()["systems"])


def system_config(
    system: str | None = None,
    config_path: Path | None = None,
) -> tuple[str, dict]:
    root = load_config(config_path)
    system = system or os.environ.get("MD_SYSTEM") or root["default"]
    try:
        return system, root["systems"][system]
    except KeyError as exc:
        choices = ", ".join(sorted(root["systems"]))
        raise KeyError(f"unknown MD analysis system '{system}'. Available: {choices}") from exc


def analysis_context(
    system: str | None = None,
    base: Path | None = None,
    payload_base: Path | None = None,
    config_path: Path | None = None,
) -> AnalysisContext:
    system, cfg = system_config(system, config_path=config_path)
    output_root = os.environ.get("MD_OUTPUT_ROOT")
    if base is None and not output_root:
        raise ValueError("Set MD_OUTPUT_ROOT to an independent output directory")
    result_base = base or Path(output_root) / system
    if payload_base is None:
        direct_payload_root = os.environ.get("MD_PAYLOAD_ROOT")
        if direct_payload_root:
            payload_base = Path(direct_payload_root)
        else:
            payload_base = Path(cfg["result_dir"])
    if not payload_base.is_absolute():
        payload_base = REPO_ROOT / payload_base
    result_resolved = result_base.resolve()
    payload_resolved = payload_base.resolve()
    if result_resolved == payload_resolved or payload_resolved in result_resolved.parents:
        raise ValueError("Output overlaps the trajectory input directory")
    return AnalysisContext(
        system=system,
        label=cfg.get("label", system),
        config=cfg,
        base=result_base,
        payload_base=payload_base,
    )


def add_system_argument(parser) -> None:
    parser.add_argument(
        "--system",
        default=None,
        choices=available_systems(),
        help="Configured MD system id; defaults to MD_SYSTEM or the fixed WT reference.",
    )


def ensure_dir(path: Path) -> Path:
    """Create an empty result directory and protect shipped inputs and references."""
    path = Path(path).resolve()
    for name in ("studies", "data", "src", "tests", "environments"):
        protected = (REPO_ROOT / name).resolve()
        if path == protected or protected in path.parents:
            raise ValueError(f"Output overlaps protected repository data: {path}")
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"Refuse nonempty output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def time_ns(ts) -> float:
    # GROMACS XTC stores time in ps; use frame index fallback if absent.
    try:
        return float(ts.time) / 1000.0
    except Exception:
        return float(ts.frame)


def in_time_window(ts, window: tuple[float, float] | None) -> bool:
    return window is None or window[0] <= time_ns(ts) <= window[1]


def trajectory_window_slice(universe, window: tuple[float, float] | None) -> slice:
    if window is None:
        return slice(None)
    times = np.array([time_ns(ts) for ts in universe.trajectory])
    indices = np.flatnonzero((times >= window[0]) & (times <= window[1]))
    if len(indices) < 2 or not np.all(np.diff(indices) == 1):
        raise ValueError(f"invalid/empty contiguous trajectory window: {window}")
    return slice(int(indices[0]), int(indices[-1]) + 1)


def equal_weight_kabsch_coordinates(
    mobile_fit_positions: np.ndarray,
    reference_fit_positions: np.ndarray,
    coordinates: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Transform row-vector coordinates with an equal-weight Kabsch fit.

    ``MDAnalysis.analysis.align.rotation_matrix`` returns a matrix acting on
    column vectors. MDAnalysis atom coordinates are stored as row vectors, so
    the returned matrix must be transposed when it is applied with ``@`` or
    ``np.dot``.
    """
    mobile_fit_positions = np.asarray(mobile_fit_positions, dtype=float)
    reference_fit_positions = np.asarray(reference_fit_positions, dtype=float)
    coordinates = np.asarray(coordinates, dtype=float)
    if mobile_fit_positions.shape != reference_fit_positions.shape:
        raise ValueError(
            "fit selection shape mismatch: "
            f"{mobile_fit_positions.shape} != {reference_fit_positions.shape}"
        )
    if mobile_fit_positions.ndim != 2 or mobile_fit_positions.shape[1] != 3 or len(mobile_fit_positions) == 0:
        raise ValueError(
            "fit selections must be non-empty N x 3 coordinate arrays: "
            f"got {mobile_fit_positions.shape}"
        )
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "coordinates must be an N x 3 array: "
            f"got {coordinates.shape}"
        )

    mobile_center = mobile_fit_positions.mean(axis=0)
    reference_center = reference_fit_positions.mean(axis=0)
    rotation, fit_rmsd = align.rotation_matrix(
        mobile_fit_positions - mobile_center,
        reference_fit_positions - reference_center,
    )
    aligned = (coordinates - mobile_center) @ rotation.T + reference_center
    return aligned, float(fit_rmsd)


def linear_slope(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return float("nan")
    return float(np.polyfit(x[mask], y[mask], 1)[0])


def save_csv(path: Path, rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    return df
