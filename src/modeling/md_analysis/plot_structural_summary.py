"""Plot local packing, ligand contacts and timing from supplied summary tables."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.patches import FancyBboxPatch
import pandas as pd

from .common import REPO_ROOT, ensure_dir

STUDY = REPO_ROOT / 'studies/receptor-engineering/molecular-dynamics'


def event_time(time_ns):
    """Format a sustained-event onset, using a dash for an absent event."""
    return '—' if pd.isna(time_ns) else f'{time_ns:.1f} ns'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, default=STUDY)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--font', type=Path, help='Chinese-capable font file')
    parser.add_argument('--language', choices=['zh', 'en'], default='zh', help='Figure label language')
    args = parser.parse_args()
    out = ensure_dir(args.output)
    if args.font:
        fontManager.addfont(args.font)
        family = FontProperties(fname=args.font).get_name()
    else:
        family = 'Noto Sans CJK SC' if args.language == 'zh' else 'DejaVu Sans'
    font = FontProperties(family=[family, 'DejaVu Sans'])
    tables = args.study / 'results/local-summary'

    def value(table, field, **keys):
        rows = pd.read_csv(tables / table)
        for key, expected in keys.items():
            rows = rows.loc[rows[key] == expected]
        if len(rows) != 1:
            raise ValueError(f'Expected one row in {table}: {keys}')
        return float(rows.iloc[0][field])

    wt18, n18 = 'native_wt_ascr18', 'raw_n84v_ascr18'
    v18 = 'raw_v368r_ascr18'
    v3 = 'v368r_ascr3_lid_rebuilt_own_pose_v1'
    wt3 = 'native_wt_ascr3_lid_rebuilt_own_pose_v1'
    pair = lambda case: value('local_pair_summary.csv', 'late_mean', case=case, site=84, neighbor=171)
    hb = lambda case, resid: value('anchor_hbond_union.csv', 'late_any_hbond_occupancy', case=case, resid=resid) * 100
    contact = lambda case, resid: value('anchor_group_occupancy.csv', 'late_occupancy', case=case, resid=resid, ligand_group='carboxylate') * 100
    charged = lambda case, field: value('Arg368_Asp70_summary.csv', field, case=case)
    bridges = pd.read_csv(tables / 'proximity_bridges.csv')
    npaths = len(bridges.loc[(bridges.case == n18) & (bridges.site == 84) & bridges.anchor.isin([186, 189])])
    joint = lambda case: value('proximity_bridges.csv', 'late_bridge_and_anchor_ligand_occupancy', case=case, site=368, neighbor=364, anchor=345) * 100
    rmsd = value('case_summary.csv', 'ligand_rmsd_core_A', case=n18)
    start = value('event_timeline.csv', 'first_1ns_start_ns', case=v3, event='residue_345_ligand_contact')

    fig = plt.figure(figsize=(16, 9.6), facecolor='white')
    ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis('off')
    teal, orange, gray = '#108681', '#B26128', '#576779'

    def text(x, y, label, size=15, color='#161C23', weight='normal'):
        ax.text(x, y, label, fontproperties=font, fontsize=size, color=color,
                weight=weight, va='top', linespacing=1.65)

    def box(x, y, title, body, color):
        ax.add_patch(FancyBboxPatch((x, y), .292, .265,
                     boxstyle='round,pad=0.008,rounding_size=0.012',
                     facecolor='#F5F8FA', edgecolor=color, linewidth=2))
        text(x + .014, y + .23, title, 16, color, 'bold')
        text(x + .014, y + .16, body, 13)

    if args.language == 'en':
        text(.025, .975, 'Local packing, ligand contacts and timing in six Gpr2 systems', 23, weight='bold')
        text(.025, .905, 'Six receptor–ligand systems; one 100 ns trajectory per system; separate starting conformations', 14, gray)
        text(.025, .84, 'N84V–ascr18: packing near residue 84 and Arg189 contacts', 18, teal)
        box(.03, .515, 'Residue 84–Gly171 distance',
            f'N84V–ascr18: {pair(n18):.2f} Å\nWT–ascr18: {pair(wt18):.2f} Å\nMean distance over 50–100 ns', teal)
        box(.355, .515, '84 to Phe186 / Arg189',
            f'Two-step paths in input table: {npaths}\nBoth protein edges ≤4.5 Å\nSequence gap >2 on each edge', teal)
        box(.68, .515, 'Arg189–carboxylate contacts',
            f'Contact occupancy: {contact(n18,189):.0f}%\nH-bond occupancy: {hb(n18,189):.1f}%\nCore-fit ligand RMSD: {rmsd:.2f} Å', teal)
        text(.025, .46, 'V368R: Arg368–Asp70 geometry, local proximity and Arg345 contacts', 18, orange)
        box(.03, .13, 'Arg368 charged end–Asp70',
            f'ascr18: {charged(v18,"late_guanidinium_Asp70_mean_A"):.2f} Å; H-bond {charged(v18,"late_any_hbond_occupancy")*100:.1f}%\n'
            f'ascr3: {charged(v3,"late_guanidinium_Asp70_mean_A"):.2f} Å; H-bond {charged(v3,"late_any_hbond_occupancy")*100:.0f}%\nCharged-end atoms: NE / NH1 / NH2', orange)
        box(.355, .13, '368–Lys364–Arg345–ligand',
            f'Joint contact: ascr3 {joint(v3):.1f}%; WT {joint(wt3):.1f}%\nCharged-end separation: {event_time(charged(v3,"first_1ns_charged_separation_ns"))}\nArg345 contact onset: {event_time(start)}', orange)
        box(.68, .13, 'Arg345–carboxylate contacts',
            f'WT–ascr3 contact: {contact(wt3,345):.0f}%\nV368R–ascr3 contact: {contact(v3,345):.0f}%\nH-bond: WT {hb(wt3,345):.1f}%; V368R {hb(v3,345):.1f}%', orange)
        text(.025, .065, 'Means / occupancy: 50–100 ns; contact ≤4 Å; two-step edges ≤4.5 Å; sustained events ≥1 ns.', 12, gray)
        text(.025, .033, 'RMSD: ligand heavy atoms vs own start; charged-end distance: Arg368 NE/NH1/NH2 to Asp70 OD1/OD2.', 12, gray)
    else:
        draw_chinese(text, box, pair, hb, contact, charged, joint, npaths, rmsd, start, wt18, n18, v18, v3, wt3, teal, orange, gray)
    for extension in ['png', 'svg']:
        path = out / f'mechanism_evidence_map.{extension}'
        fig.savefig(path, dpi=160)
        if extension == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
    plt.close(fig)


def draw_chinese(text, box, pair, hb, contact, charged, joint, npaths, rmsd, start, wt18, n18, v18, v3, wt3, teal, orange, gray):
    text(.025, .975, '六组 Gpr2 模拟的局部堆积、配体接触与时间', 25, weight='bold')
    text(.025, .905, '六组受体–配体体系；每组一条 100 ns 轨迹；各自起始构象', 14, gray)
    text(.025, .84, 'N84V–ascr18：84 位附近堆积与 Arg189 配体接触', 18, teal)
    box(.03, .515, '84 位与 Gly171 的距离',
        f'N84V–ascr18：{pair(n18):.2f} Å\nWT–ascr18：{pair(wt18):.2f} Å\n50–100 ns 的平均距离', teal)
    box(.355, .515, '84 位至 Phe186 / Arg189',
        f'输入表中的两步路径：{npaths} 条\n两条蛋白近接边均 ≤4.5 Å\n每条边的序列间隔均 >2', teal)
    box(.68, .515, 'Arg189–羧酸端接触',
        f'近接占有率：{contact(n18,189):.0f}%\n任一几何氢键：{hb(n18,189):.1f}%\n核心拟合 RMSD：{rmsd:.2f} Å', teal)
    text(.025, .46, 'V368R：Arg368–Asp70 几何、局部近接与 Arg345 接触', 18, orange)
    box(.03, .13, 'Arg368 带电端–Asp70',
        f'ascr18：{charged(v18,"late_guanidinium_Asp70_mean_A"):.2f} Å；氢键 {charged(v18,"late_any_hbond_occupancy")*100:.1f}%\n'
        f'ascr3：{charged(v3,"late_guanidinium_Asp70_mean_A"):.2f} Å；氢键 {charged(v3,"late_any_hbond_occupancy")*100:.0f}%\n带电端原子：NE / NH1 / NH2', orange)
    box(.355, .13, '368–Lys364–Arg345–配体',
        f'联合近接：ascr3 {joint(v3):.1f}%；WT {joint(wt3):.1f}%\n带电端分离起点：{event_time(charged(v3,"first_1ns_charged_separation_ns"))}\nArg345 近接起点：{event_time(start)}', orange)
    box(.68, .13, 'Arg345–羧酸端接触对照',
        f'WT–ascr3：近接 {contact(wt3,345):.0f}%\nV368R–ascr3：近接 {contact(v3,345):.0f}%\n几何氢键：WT {hb(wt3,345):.1f}%；V368R {hb(v3,345):.1f}%', orange)
    text(.025, .065, '均值及占有率：50–100 ns；接触 ≤4 Å；两步近接边 ≤4.5 Å；持续事件 ≥1 ns。', 12, gray)
    text(.025, .033, 'RMSD：配体重原子相对各自起点；带电端距离：Arg368 NE/NH1/NH2 至 Asp70 OD1/OD2。', 12, gray)


if __name__ == '__main__':
    main()
