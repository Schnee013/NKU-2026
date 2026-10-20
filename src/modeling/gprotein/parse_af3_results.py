"""Extract AF3 confidence tables using the retained scoring definitions."""
import argparse
import csv
import json
import math
from pathlib import Path

from ..af3 import best_prediction, ca_plddts, one_file, predictions

GRAFT_LENGTHS = {
    "Gpa1-WT": 0,
    "Gpa1-AoGpa2-C5": 5, "Gpa1-AoGpa2-C10": 10, "Gpa1-AoGpa2-C15": 15,
    "Gpa1-AoGpa2-C25": 25, "Gpa1-AoGpa2-C37": 37,
    "Gpa1-AfGasA-C5": 5, "Gpa1-AfGasA-C10": 10, "Gpa1-AfGasA-C15": 15,
    "Gpa1-AfGasA-C25": 25, "Gpa1-AfGasA-C37": 37,
}
COGNATE_PAIRS = {("Gpr2", "AoGpa2"), ("Gpr3", "AoGpa2"), ("GprC", "AfGasA")}
STUDY = Path(__file__).resolve().parents[3] / "studies" / "gprotein-coupling"


def job_name_to_folder(name: str) -> str:
    """Normalize the retained job IDs used by numbered exports."""
    return name.lower().replace("-", "_")


def discover_jobs(downloads: Path) -> list[dict]:
    """Match native input names, or numbered export folders, to the job table."""
    with (STUDY / "prediction-jobs.csv").open(newline="") as handle:
        expected = list(csv.DictReader(handle))
    by_name = {r["job_id"]: r for r in expected}
    by_folder = {job_name_to_folder(r[key]): r for r in expected for key in ("job_id", "score_id")}
    jobs = []
    for folder in sorted(downloads.iterdir()):
        if not folder.is_dir():
            continue
        data_path = one_file(folder, "*_data.json", optional=True)
        if data_path:
            data = json.loads(data_path.read_text())
            record = by_name.get(data.get("name"))
            if record is not None:
                request = json.loads((STUDY / record["input_path"]).read_text())
                chains = lambda doc: [(p["protein"]["id"], p["protein"]["sequence"]) for p in doc["sequences"]]
                if chains(data) != chains(request):
                    raise ValueError(f"Prediction sequences differ from {record['job_id']}")
        else:
            record = by_folder.get(folder.name)
        if record is None:
            raise ValueError(f"Unrecognized prediction job: {folder}")
        job = {"name": record["score_id"], "type": record["type"],
               "chimera": record["chimera"], "folder": folder}
        if job["type"] == "complex":
            receptor = record["score_id"].split("_", 1)[0]
            donor = "AoGpa2" if "AoGpa2" in job["chimera"] else "AfGasA" if "AfGasA" in job["chimera"] else "WT"
            job.update(receptor=receptor, pairing="cognate" if (receptor, donor) in COGNATE_PAIRS
                       else "control" if job["chimera"] == "Gpa1-WT" else "cross")
        jobs.append(job)
    observed = [j["name"] for j in jobs]
    wanted = {r["score_id"] for r in expected}
    if len(observed) != len(set(observed)) or set(observed) != wanted:
        raise ValueError(f"prediction jobs incomplete/mismatched: missing={sorted(wanted-set(observed))}, "
                         f"extra={sorted(set(observed)-wanted)}, duplicate jobs={len(observed)-len(set(observed))}")
    return jobs


def find_best_model(folder: Path):
    """Return the selected prediction with its seed/sample or export index."""
    return best_prediction(predictions(folder))


def compute_interface_pae(full_data: dict) -> float:
    """Mean over all directed token pairs belonging to different chains."""
    pae = full_data.get("pae", [])
    chains = full_data.get("token_chain_ids", [])
    if not pae or not chains:
        return float("inf")
    n = len(chains)
    if len(pae) != n or any(len(row) != n for row in pae):
        raise ValueError("PAE matrix dimensions differ from token_chain_ids")
    if any(not math.isfinite(value) or value < 0 for row in pae for value in row):
        raise ValueError("Invalid PAE value")
    values = [pae[i][j] for i in range(n) for j in range(n) if chains[i] != chains[j]]
    return sum(values) / len(values) if values else float("inf")


def analyze_monomers(jobs: list[dict]) -> list[dict]:
    results = []
    for job in jobs:
        if job["type"] != "monomer":
            continue
        models = predictions(job["folder"])
        best = best_prediction(models)
        summary = best.summary()
        values = ca_plddts(best.cif_path)
        graft_len = GRAFT_LENGTHS[job["chimera"]]
        if len(values) != 472:
            raise ValueError(f"Expected 472 Gpa1 C-alpha records for {job['name']}; found {len(values)}")
        scores = [model.summary()["ranking_score"] for model in models]
        row = {
            "name": job["name"], "chimera": job["chimera"], "graft_len": graft_len,
            **best.identity(),
            "plddt_global": round(sum(values) / len(values), 2),
            "plddt_cterminal": round(sum(values[-graft_len:]) / graft_len, 2) if graft_len else "N/A",
            "ptm": round(summary["ptm"], 4),
            "ranking_score": round(summary["ranking_score"], 4),
            "fraction_disordered": round(summary["fraction_disordered"], 3),
            "score_range": f"{min(scores):.4f}-{max(scores):.4f}",
        }
        results.append(row)
    return results


def analyze_complexes(jobs: list[dict]) -> list[dict]:
    results = []
    for job in jobs:
        if job["type"] != "complex":
            continue
        models = predictions(job["folder"])
        best = best_prediction(models)
        summary = best.summary()
        ipae = compute_interface_pae(best.confidence())
        if not math.isfinite(ipae):
            raise ValueError(f"Missing cross-chain PAE for {job['name']}")
        cppm = summary.get("chain_pair_pae_min", [])
        min_cross_pae = round(cppm[0][1], 2) if len(cppm) >= 2 and len(cppm[0]) >= 2 else "N/A"
        iptms = [s["iptm"] for model in models if (s := model.summary()).get("iptm") is not None]
        row = {
            "name": job["name"], "receptor": job["receptor"], "chimera": job["chimera"],
            "pairing": job["pairing"], **best.identity(),
            "iptm": round(summary.get("iptm") or 0, 4), "ptm": round(summary["ptm"], 4),
            "ipae_mean": round(ipae, 2), "chain_pair_pae_min": min_cross_pae,
            "ranking_score": round(summary["ranking_score"], 4),
            "fraction_disordered": round(summary["fraction_disordered"], 3),
            "has_clash": summary["has_clash"],
            "iptm_range": f"{min(iptms):.4f}-{max(iptms):.4f}" if iptms else "N/A",
        }
        results.append(row)
    return results


def write_csv(rows: list[dict], output_path: Path):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--downloads", type=Path, required=True, help="Directory containing the 42 local AF3 job outputs")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args(argv)
    downloads, out = args.downloads.resolve(), args.output.resolve()
    if not downloads.is_dir():
        ap.error("required 42-job AF3 downloads directory does not exist")
    if out.is_relative_to(STUDY) or out.is_relative_to(downloads):
        ap.error("output must be outside study references and prediction inputs")
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        ap.error("output must be a new or empty directory")
    try:
        jobs = discover_jobs(downloads)
        monomers, complexes = analyze_monomers(jobs), analyze_complexes(jobs)
    except (ValueError, OSError, KeyError, TypeError) as exc:
        ap.error(str(exc))
    if len(monomers) != 11 or len(complexes) != 31:
        ap.error("incomplete score extraction")
    out.mkdir(parents=True, exist_ok=True)
    write_csv(monomers, out / "monomer_scores.csv")
    write_csv(complexes, out / "complex_scores.csv")
    print(f"Extracted {len(monomers)} monomer and {len(complexes)} complex scores to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
