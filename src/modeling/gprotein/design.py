"""Generate terminal sequence grafts from Gpa1 and donor FASTA records."""
import argparse
import csv
from pathlib import Path

STUDY = Path(__file__).resolve().parents[3] / "studies" / "gprotein-coupling"

def read_fasta(path):
    lines=Path(path).read_text().splitlines()
    if sum(line.startswith(">") for line in lines)!=1:
        raise ValueError(f"one sequence required: {path}")
    return "".join(line.strip() for line in lines if not line.startswith(">"))

def designs(study=STUDY):
    wt=read_fasta(study/"sequences/yeast_Gpa1_P08539.fasta")
    if len(wt)!=472:
        raise ValueError("Gpa1 backbone must contain 472 residues")
    yield "Gpa1-WT",wt,"ScGpa1",0,472
    for donor,file in [("AoGpa2","Ao_Gpa2_TGJ64397.fasta"),("AfGasA","Af_GasA_A0A436ZRX2.fasta")]:
        seq=read_fasta(study/"sequences"/file)
        for n in [5,10,15,25,37]:
            yield f"Gpa1-{donor}-C{n}",wt[:-n]+seq[-n:],donor,n,472-n

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args();out=args.output.resolve()
    if out.is_relative_to(STUDY) or (out.exists() and (not out.is_dir() or any(out.iterdir()))):
        ap.error("output must be a new/empty directory outside study references")
    rows=list(designs())
    out.mkdir(parents=True,exist_ok=True)
    with (out/"designs.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["design_id","donor","graft_length","backbone_end_1based","length"])
        for name,seq,donor,n,end in rows:
            w.writerow([name,donor,n,end,len(seq)])
            (out/(name+".fasta")).write_text(f">{name}\n{seq}\n")
if __name__=="__main__":
    main()
