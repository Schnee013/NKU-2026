"""Read local AlphaFold 3 inputs and native or numbered prediction exports."""
from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re


def single_protein_sequence(document: dict | list) -> str:
    """Read one protein copy from a native input or a compatible older request."""
    if isinstance(document, list):
        if len(document) != 1:
            raise ValueError("Expected one prediction request")
        document = document[0]
    if not isinstance(document, dict):
        raise ValueError("Expected one prediction request object")
    entities = document.get("sequences", [])
    if len(entities) != 1:
        raise ValueError("Expected one protein chain")
    entity = entities[0]
    if "protein" in entity:
        if document.get("dialect") != "alphafold3":
            raise ValueError("Expected alphafold3 input dialect")
        protein = entity["protein"]
        ids = protein.get("id")
        ids = [ids] if isinstance(ids, str) else ids
        if not isinstance(ids, list) or len(ids) != 1 or not isinstance(ids[0], str) or not ids[0]:
            raise ValueError("Expected one explicit protein chain ID")
    elif "proteinChain" in entity:
        protein = entity["proteinChain"]
        if protein.get("count", 1) != 1:
            raise ValueError("Expected one protein chain copy")
    else:
        raise ValueError("Expected one protein chain")
    sequence = protein.get("sequence")
    if not isinstance(sequence, str) or not sequence:
        raise ValueError("Missing protein sequence")
    return sequence


def one_file(folder: Path, pattern: str, *, optional: bool = False) -> Path | None:
    files = sorted(folder.glob(pattern))
    if not files and optional:
        return None
    if len(files) != 1:
        raise ValueError(f"Expected one {pattern} in {folder}; found {len(files)}")
    return files[0]


@dataclass(frozen=True)
class Prediction:
    ranking_score: float
    summary_path: Path
    cif_path: Path | None
    confidence_path: Path | None
    seed: int | None = None
    sample: int | None = None
    model_index: int | None = None

    def summary(self) -> dict:
        return json.loads(self.summary_path.read_text())

    def confidence(self) -> dict:
        if self.confidence_path is None:
            raise ValueError(f"Missing full confidence file for {self.summary_path}")
        return json.loads(self.confidence_path.read_text())

    def identity(self) -> dict:
        return {"best_model": self.model_index if self.model_index is not None else "",
                "best_seed": self.seed if self.seed is not None else "",
                "best_sample": self.sample if self.sample is not None else ""}


def _score(value: object) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite AF3 ranking score")
    return number


def predictions(folder: Path) -> list[Prediction]:
    """Read native CSV order or the five retained numbered export indices."""
    ranking = one_file(folder, "*ranking_scores.csv", optional=True)
    if ranking is not None:
        result, seen = [], set()
        with ranking.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if not {"seed", "sample", "ranking_score"} <= set(reader.fieldnames or []):
                raise ValueError(f"Invalid ranking CSV columns: {ranking}")
            for row in reader:
                seed, sample = int(row["seed"]), int(row["sample"])
                if seed < 0 or sample < 0 or (seed, sample) in seen:
                    raise ValueError(f"Invalid or duplicate seed/sample in {ranking}")
                seen.add((seed, sample))
                directory = folder / f"seed-{seed}_sample-{sample}"
                summary = one_file(directory, "*summary_confidences.json")
                confidence_files = [p for p in directory.glob("*confidences.json")
                                    if not p.name.endswith("summary_confidences.json")]
                if len(confidence_files) != 1:
                    raise ValueError(f"Expected one full confidence file in {directory}")
                result.append(Prediction(_score(row["ranking_score"]), summary,
                                         one_file(directory, "*model.cif"), confidence_files[0],
                                         seed=seed, sample=sample))
        if not result:
            raise ValueError(f"Empty ranking CSV: {ranking}")
        # A partial CSV must not silently omit sampled predictions.
        directories = {p.name for p in folder.glob("seed-*_sample-*") if p.is_dir()}
        if directories != {f"seed-{s}_sample-{n}" for s, n in seen}:
            raise ValueError(f"Ranking CSV and sample directories differ: {folder}")
        data_path = one_file(folder, "*_data.json", optional=True)
        if data_path:
            data = json.loads(data_path.read_text())
            if data.get("dialect") != "alphafold3" or set(data.get("modelSeeds", [])) != {s for s, _ in seen}:
                raise ValueError(f"Input dialect/seeds and ranking CSV differ: {folder}")
        return result
    if any(folder.glob("seed-*_sample-*")) or any(folder.glob("*_data.json")):
        raise ValueError(f"Missing native ranking CSV in {folder}")
    result = []
    for index in range(5):
        summary = one_file(folder, f"*summary_confidences_{index}.json")
        data = json.loads(summary.read_text())
        result.append(Prediction(_score(data["ranking_score"]), summary,
                                 one_file(folder, f"*model_{index}.cif", optional=True),
                                 one_file(folder, f"*full_data_{index}.json", optional=True),
                                 model_index=index))
    return result


def best_prediction(models: list[Prediction]) -> Prediction:
    """Maximize full CSV precision; ties retain the original input order."""
    return max(models, key=lambda model: model.ranking_score)


def ca_plddts(path: Path | None) -> list[float]:
    """Read C-alpha pLDDT by mmCIF column names, preserving residue order."""
    if path is None:
        raise ValueError("Missing selected model CIF")
    columns, tokens, values = [], [], []
    active = False
    token_pattern = re.compile(r"'[^']*'(?=\s|$)|\"[^\"]*\"(?=\s|$)|\S+")
    for line in path.read_text().splitlines():
        line = line.strip()
        if line.startswith("_atom_site."):
            if tokens:
                raise ValueError(f"Malformed atom-site loop: {path}")
            columns.append(line.split()[0].removeprefix("_atom_site."))
            active = True
            continue
        if not active or not line:
            continue
        if line.startswith(("#", "_", "loop_", "data_")):
            if tokens:
                raise ValueError(f"Incomplete atom-site row: {path}")
            break
        tokens.extend(token[1:-1] if token[:1] in {"'", '"'} else token
                      for token in token_pattern.findall(line))
        while len(tokens) >= len(columns):
            row = dict(zip(columns, tokens[:len(columns)]))
            del tokens[:len(columns)]
            if row.get("group_PDB") == "ATOM" and row.get("label_atom_id", row.get("auth_atom_id")) == "CA":
                value = _score(row["B_iso_or_equiv"])
                if not 0 <= value <= 100:
                    raise ValueError(f"Invalid C-alpha pLDDT in {path}")
                values.append(value)
    if tokens or not values:
        raise ValueError(f"Missing or malformed C-alpha confidence records: {path}")
    return values
