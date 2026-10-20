"""Plot the anchor question from completed polar tables, with exact source data."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from .common import ensure_dir

def main():
    p=argparse.ArgumentParser();p.add_argument('--summary-dir',type=Path,required=True);p.add_argument('--data-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    ensure_dir(a.out)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
    ids=['v368r_ascr3_arg368_start','v368r_ascr3_arg189_start','native_wt_ascr3_lid_rebuilt_own_pose_v1','raw_v368r_ascr18','v368r_ascr3_lid_rebuilt_own_pose_v1']
    labels=['Arg368-directed start','Arg189-directed start','WT–ascr3','V368R–ascr18','Original V368R–ascr3']
    colors=['#8254a0','#087f8c','#a4adb5','#a4adb5','#a4adb5'];common='#386ba0'
    fig=plt.figure(figsize=(12.8,7.2));gs=fig.add_gridspec(2,2,width_ratios=[1.25,1],left=.08,right=.97,bottom=.20,top=.86,wspace=.55,hspace=.65)
    plotrows=[]
    for j,target in enumerate([368,189]):
     ax=fig.add_subplot(gs[j,0]);df=pd.read_csv(a.data_root/ids[j]/'polar/polar_timeseries.csv.gz')
     t=df.time_ns.to_numpy();b=np.minimum(t.astype(int),99)
     for feature,color,label in [(f'R{target}_ligand',colors[j],f'Arg{target} (designed)'),('R345_ligand',common,'Arg345')]:
      v=df[feature+'__salt_4A'].to_numpy()
      bins=np.array([v[b==k].mean() for k in range(100)])
      ax.stairs(bins*100,np.arange(101),color=color,lw=2.2,label=label,baseline=None)
      for k,value in enumerate(bins):plotrows.append(dict(case_id=ids[j],feature=feature,start_ns=k,end_ns=k+1,occupancy=value))
     ax.axvspan(50,100,color='#e8edf2',zorder=-2);ax.set(xlim=(0,100),ylim=(-4,105),yticks=[0,50,100],xlabel='Time (ns)',ylabel='Salt-bridge occupancy (%)')
     ax.text(-.16,1.1,'AB'[j],transform=ax.transAxes,fontweight='bold',fontsize=16)
     ax.set_title(labels[j],loc='left',fontweight='bold',pad=14);ax.legend(loc='center right',fontsize=10,frameon=False)
     if j==0:ax.text(74,6,'50–100 ns',ha='center',fontsize=9,color='#58616a')
    ax=fig.add_subplot(gs[:,1]);summary=pd.read_csv(a.summary_dir/'primary_polar_summary.csv');data=[]
    for cid in ids:
     row={}
     for metric in ['salt_4A','short_linear_hbond']:
      row[metric]=float(summary[(summary.case_id==cid)&(summary.feature=='R345_ligand')&(summary.metric==metric)].occupancy.iloc[0])*100
     data.append(row)
    y=np.arange(5);salt=[v['salt_4A'] for v in data];hb=[v['short_linear_hbond'] for v in data]
    ax.barh(y,salt,color=colors,height=.5,alpha=.85)
    ax.scatter(hb,y,s=48,marker='D',color='#25374b',zorder=3)
    for k,value in enumerate(salt):ax.text(value+1.5,k,f'{value:.2f}%',va='center',fontsize=10)
    ax.set(xlim=(0,118),xticks=[0,50,100],yticks=y,yticklabels=labels,xlabel='Occupancy, 50–100 ns (%)');ax.invert_yaxis();ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0,pad=9)
    ax.text(-.55,1.05,'C',transform=ax.transAxes,fontweight='bold',fontsize=16);ax.set_title('Arg345 anchoring is shared',fontweight='bold',pad=20)
    fig.legend(handles=[Line2D([0],[0],color='#9eaab6',lw=8,label='Salt bridge (≤4 Å)'),Line2D([0],[0],marker='D',color='none',markerfacecolor='#25374b',markeredgecolor='#25374b',label='Short linear H bond')],loc='lower center',bbox_to_anchor=(.73,.06),ncol=2,frameon=False,fontsize=10)
    fig.suptitle('Designed anchors are lost; Arg345 anchoring is shared with controls',fontsize=15,fontweight='bold',y=.98)
    fig.text(.08,.015,'A–B: 1 ns bins from all 10 ps samples. C: same 50–100 ns window in all five trajectories. Each case is one trajectory.',fontsize=9,color='#50555b')
    for suffix in ['png','svg','pdf']:fig.savefig(a.out/f'figure1_anchor_persistence.{suffix}',dpi=210,facecolor='white')
    plt.close(fig)
    pd.DataFrame(plotrows).to_csv(a.out/'figure1_timeline_values.csv',index=False)
    pd.DataFrame([dict(case_id=cid,**row) for cid,row in zip(ids,data)]).to_csv(a.out/'figure1_comparison_values.csv',index=False)
    (a.out/'figure1_sources.json').write_text(json.dumps(dict(figure='figure1_anchor_persistence',definition='Frozen 4A salt-bridge and short linear hydrogen-bond geometry',case_ids=ids,bin_ns=1),indent=2)+'\n')

if __name__ == '__main__':
    main()
