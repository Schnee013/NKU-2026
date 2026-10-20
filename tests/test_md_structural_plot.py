"""Exercise structural-summary figures with synthetic CSV measurements."""
from __future__ import annotations

import csv
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def write_summaries(study: Path, *, nonzero: bool) -> None:
    """Create independent, known-value inputs for every plotted measurement."""
    tables = study / "results/local-summary"
    tables.mkdir(parents=True)
    wt18, n18 = "native_wt_ascr18", "raw_n84v_ascr18"
    v18 = "raw_v368r_ascr18"
    wt3 = "native_wt_ascr3_lid_rebuilt_own_pose_v1"
    v3 = "v368r_ascr3_lid_rebuilt_own_pose_v1"

    def write(name, fields, rows):
        with (tables / name).open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(fields.split(","))
            writer.writerows(rows)

    write("local_pair_summary.csv", "case,site,neighbor,late_mean", [
        (n18, 84, 171, 4.25 if nonzero else 5.0),
        (wt18, 84, 171, 6.75 if nonzero else 5.0),
    ])
    write("anchor_hbond_union.csv", "case,resid,late_any_hbond_occupancy", [
        (n18, 189, 0.125 if nonzero else 0.0),
        (wt3, 345, 0.25 if nonzero else 0.0),
        (v3, 345, 0.5 if nonzero else 0.0),
    ])
    write("anchor_group_occupancy.csv", "case,resid,ligand_group,late_occupancy", [
        (n18, 189, "carboxylate", 0.6 if nonzero else 0.0),
        (wt3, 345, "carboxylate", 0.35 if nonzero else 0.0),
        (v3, 345, "carboxylate", 0.8 if nonzero else 0.0),
    ])
    write("Arg368_Asp70_summary.csv",
          "case,late_guanidinium_Asp70_mean_A,late_any_hbond_occupancy,first_1ns_charged_separation_ns", [
        (v18, 8.5, 0.2 if nonzero else 0.0, None),
        (v3, 9.75, 0.4 if nonzero else 0.0, 2.5 if nonzero else None),
    ])
    bridges = [
        (v3, 368, 364, 345, 0.45 if nonzero else 0.0),
        (wt3, 368, 364, 345, 0.15 if nonzero else 0.0),
    ]
    if nonzero:
        bridges.extend([(n18, 84, 150, 186, 0.3), (n18, 84, 151, 189, 0.4)])
    write("proximity_bridges.csv",
          "case,site,neighbor,anchor,late_bridge_and_anchor_ligand_occupancy", bridges)
    write("case_summary.csv", "case,ligand_rmsd_core_A", [(n18, 1.75)])
    write("event_timeline.csv", "case,event,first_1ns_start_ns", [
        (v3, "residue_345_ligand_contact", 12.5 if nonzero else None),
    ])


class StructuralSummaryPlotTests(unittest.TestCase):
    def render(self, root: Path, language: str) -> str:
        output = root / language
        env = dict(os.environ, MPLCONFIGDIR=str(root / "mpl"))
        result = subprocess.run(
            [sys.executable, "-m", "src.modeling.md_analysis.plot_structural_summary",
             "--study", str(root / "study"), "--output", str(output),
             "--language", language],
            cwd=ROOT, env=env, capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        png = output / "mechanism_evidence_map.png"
        self.assertEqual(png.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
        self.assertGreater(png.stat().st_size, 1000)
        svg = output / "mechanism_evidence_map.svg"
        self.assertEqual(ET.parse(svg).getroot().tag, "{http://www.w3.org/2000/svg}svg")
        # Matplotlib retains label text in SVG comments when glyphs are paths.
        return "\n".join(re.findall(r"<!--\s*(.*?)\s*-->", svg.read_text(), re.S))

    def test_zero_contacts_equal_distances_and_absent_events(self):
        with tempfile.TemporaryDirectory(prefix="md summary zero ") as tmp:
            root = Path(tmp)
            write_summaries(root / "study", nonzero=False)
            for language in ("en", "zh"):
                with self.subTest(language=language):
                    labels = self.render(root, language)
                    colon = ": " if language == "en" else "："
                    self.assertIn(f"N84V–ascr18{colon}5.00 Å", labels)
                    self.assertIn(f"WT–ascr18{colon}5.00 Å", labels)
                    if language == "en":
                        expected = ("WT–ascr3 contact: 0%", "V368R–ascr3 contact: 0%",
                                    "Arg345 contact onset: —", "Charged-end separation: —",
                                    "Two-step paths in input table: 0",
                                    "Mean distance over 50–100 ns",
                                    "Contact occupancy: 0%", "H-bond occupancy: 0.0%",
                                    "Core-fit ligand RMSD: 1.75 Å")
                    else:
                        expected = ("WT–ascr3：近接 0%", "V368R–ascr3：近接 0%",
                                    "Arg345 近接起点：—", "带电端分离起点：—",
                                    "输入表中的两步路径：0 条",
                                    "50–100 ns 的平均距离",
                                    "近接占有率：0%", "任一几何氢键：0.0%",
                                    "核心拟合 RMSD：1.75 Å")
                    for label in expected:
                        self.assertIn(label, labels)
                    self.assertNotIn("nan ns", labels)

    def test_nonzero_measurements_and_event_times_follow_the_tables(self):
        with tempfile.TemporaryDirectory(prefix="md summary values ") as tmp:
            root = Path(tmp)
            write_summaries(root / "study", nonzero=True)
            for language in ("en", "zh"):
                with self.subTest(language=language):
                    labels = self.render(root, language)
                    if language == "en":
                        expected = (
                            "N84V–ascr18: 4.25 Å", "WT–ascr18: 6.75 Å",
                            "Contact occupancy: 60%", "H-bond occupancy: 12.5%",
                            "Core-fit ligand RMSD: 1.75 Å",
                            "WT–ascr3 contact: 35%", "V368R–ascr3 contact: 80%",
                            "H-bond: WT 25.0%; V368R 50.0%",
                            "Joint contact: ascr3 45.0%; WT 15.0%",
                            "Charged-end separation: 2.5 ns", "Arg345 contact onset: 12.5 ns",
                            "Two-step paths in input table: 2",
                        )
                    else:
                        expected = (
                            "N84V–ascr18：4.25 Å", "WT–ascr18：6.75 Å",
                            "近接占有率：60%", "任一几何氢键：12.5%",
                            "核心拟合 RMSD：1.75 Å",
                            "WT–ascr3：近接 35%", "V368R–ascr3：近接 80%",
                            "几何氢键：WT 25.0%；V368R 50.0%",
                            "联合近接：ascr3 45.0%；WT 15.0%",
                            "带电端分离起点：2.5 ns", "Arg345 近接起点：12.5 ns",
                            "输入表中的两步路径：2 条",
                        )
                    for label in expected:
                        self.assertIn(label, labels)


if __name__ == "__main__":
    unittest.main()
