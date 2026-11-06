# robustness of grade A to threshold choices
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from samples import *
import itertools
S=pd.read_csv(f"{OUT}/probe_criteria.csv",index_col=0)
gp=pd.read_csv(f"{OUT}/gene_probe_map.csv")
sig=S.paad_mean-S.blood_mean
rows=[]
for bm,fp,pf,pm in itertools.product([0.05,0.10,0.15],[0.05,0.10,0.20],[0.5,1.0],[2,4,6]):
    c2=(S.blood_median<=bm)&(S.blood_p95<=2*bm)
    c3=S[f"sens_f{fp}"]>=0.5
    npos=sum(((S[f"mean_{c}"]-S.blood_mean)>=pf*sig).astype(int) for c in PAN)
    c4=(sig>0)&(npos<=pm)
    ok=S.c1_repro&c2&c3&c4
    genes=sorted(gp[gp.probe.isin(S.index[ok.fillna(False)])].gene.unique())
    rows.append(dict(blood_median_max=bm,dilution_fraction=fp,pan_signal_frac=pf,max_other_cancers=pm,n_gradeA=len(genes),gradeA_genes=";".join(genes)))
T=pd.DataFrame(rows); T.to_csv(f"{OUT}/robustness.csv",index=False)
print(T.n_gradeA.describe()); 
from collections import Counter
c=Counter(g for x in T.gradeA_genes for g in x.split(";") if g); print("times grade A across",len(T),"settings:",c.most_common(15))
base=T[(T.blood_median_max==0.10)&(T.dilution_fraction==0.10)&(T.pan_signal_frac==0.5)&(T.max_other_cancers==4)]
print(base.to_string(index=False))
