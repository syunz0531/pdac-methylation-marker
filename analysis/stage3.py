# marker combination in simulated cfDNA: train thresholds on GSE49149, test on TCGA-PAAD tumours
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from samples import *
rng=np.random.default_rng(7)
S=pd.read_csv(f"{OUT}/probe_criteria.csv",index_col=0)
G=pd.read_csv(f"{OUT}/gene_grades.csv")
C,blood,clin=cohorts()
bt=read_targets("GSE40279")[blood].astype(float)
perm=rng.permutation(len(blood)); bref=[blood[i] for i in perm[:328]]; btest=[blood[i] for i in perm[328:]]
tr=read_targets("GSE49149")[C["GSE49149"][0]].astype(float)
te=read_targets("TCGA_PAAD")[C["TCGA-PAAD"][0]].astype(float)
def panel_eval(probes, fgrid=(0.2,0.1,0.05,0.02,0.01), rep=20):
    probes=[p for p in probes if p in tr.index and p in te.index and p in bt.index]
    thr=bt.loc[probes,bref].quantile(0.99,axis=1)
    out={}
    neg=bt.loc[probes,btest].to_numpy()
    for f in fgrid:
        res={}
        for name,T in [("train_GSE49149",tr),("test_TCGA",te)]:
            t=T.loc[probes].to_numpy()                       # probes x tumours
            ok=~np.isnan(t).any(0); t=t[:,ok]
            hits=[]
            for _ in range(rep):
                b=neg[:,rng.integers(0,neg.shape[1],t.shape[1])]
                mix=f*t+(1-f)*b
                hits.append(((mix-thr.to_numpy()[:,None])>0).sum(0))
            hits=np.concatenate(hits)
            # negatives: held-out healthy blood
            nh=((neg-thr.to_numpy()[:,None])>0).sum(0)
            k=int(np.quantile(nh,0.99))+1                    # positives need > 99th pct of negatives' hit count
            res[name]=(hits>=k).mean()
            res["spec_heldout_blood"]=(nh<k).mean(); res["k_required"]=k
        out[f]=res
    return probes,out
best=lambda genes:[G.set_index("gene").loc[g,"best_probe"] for g in genes]
A=G[G.grade=="A"].gene.tolist()
B_blood=G[(G.grade=="B")&(~G.best_fails.fillna("").str.contains("blood"))].sort_values("best_sens_f10",ascending=False).gene.tolist()
sets={"best single (A)":A[:1],"grade A":A,"A + top3 B":A+B_blood[:3],"A + top5 B":A+B_blood[:5],"A + all blood-suitable B":A+B_blood}
rows=[]
for name,genes in sets.items():
    probes,out=panel_eval(best(genes))
    for f,r in out.items():
        rows.append(dict(panel=name,genes=";".join(genes),n=len(probes),tumor_fraction=f,**r))
T=pd.DataFrame(rows); T.to_csv(f"{OUT}/panel_performance.csv",index=False)
print(T.round(3).to_string(index=False))
