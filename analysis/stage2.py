import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from samples import *
rng=np.random.default_rng(20261106)
gp=pd.read_csv(f"{OUT}/gene_probe_map.csv")
R=pd.read_csv(f"{OUT}/tissue_probe_stats.csv")
C,blood,clin=cohorts()
# ---- thresholds (fixed before evaluating criteria 2-4) ----
DELTA,AUC,FDR=0.10,0.70,0.05
BLOOD_MED,BLOOD_P95=0.10,0.20
F_GRID=[1.0,0.5,0.2,0.1,0.05,0.02,0.01]; F_PASS=0.10; SENS_PASS=0.50
PAN_FRAC,PAN_MAX=0.5,4
# 1 reproducibility
P=R.pivot(index="probe",columns="cohort")
ok=lambda c:(P[("delta",c)]>=DELTA)&(P[("auc",c)]>=AUC)&(P[("fdr",c)]<FDR)
S=pd.DataFrame(index=P.index)
S["repro_TCGA"]=ok("TCGA-PAAD"); S["repro_GSE49149"]=ok("GSE49149")
S["dir_GSE74071"]=P[("delta","GSE74071")]>0
S["measured_both"]=P[("delta","TCGA-PAAD")].notna()&P[("delta","GSE49149")].notna()
S["c1_repro"]=S.repro_TCGA&S.repro_GSE49149
for c in ["TCGA-PAAD","GSE49149","GSE74071"]:
    S[f"auc_{c}"]=P[("auc",c)]; S[f"delta_{c}"]=P[("delta",c)]
# 2 blood
B=pd.read_csv(f"{OUT}/blood_probe_stats.csv",index_col=0)
S=S.join(B)
S["c2_blood"]=(S.blood_median<=BLOOD_MED)&(S.blood_p95<=BLOOD_P95)
# 3 dilution: tumour tissue (GSE49149) mixed into blood; threshold = 99th pct of a separate half of blood
bt=read_targets("GSE40279")[blood].astype(float)
tt=read_targets("GSE49149")[C["GSE49149"][0]].astype(float)
_m=pd.read_csv(f"{D}/GSE49149_samples.csv")
_cel=pd.to_numeric(_m.Sample_characteristics_ch1_2.str.replace("cellularity: ",""),errors="coerce")
HI=[x for x,c in zip(_m["sample"],_cel) if c>=40 and x in tt.columns]
sens_hi={}
perm=rng.permutation(len(blood)); ref=[blood[i] for i in perm[:328]]; mixb=[blood[i] for i in perm[328:]]
sens={f:{} for f in F_GRID}; tau={}
for p in S.index:
    if p not in tt.index or p not in bt.index: continue
    r=bt.loc[p,ref].dropna().to_numpy(); m=bt.loc[p,mixb].dropna().to_numpy(); t=tt.loc[p].dropna().to_numpy()
    if len(r)<50 or len(t)<20: continue
    th=np.quantile(r,0.99); tau[p]=th
    bj=rng.choice(m,size=(20,len(t)))
    for f in F_GRID:
        sens[f][p]=float(((f*t+(1-f)*bj)>th).mean())
    th_=tt.loc[p,HI].dropna().to_numpy(); bh=rng.choice(m,size=(20,len(th_)))
    sens_hi[p]=float(((0.1*th_+0.9*bh)>th).mean())
S["blood_thr99"]=pd.Series(tau)
for f in F_GRID: S[f"sens_f{f}"]=pd.Series(sens[f])
S["sens_hicell_f0.1"]=pd.Series(sens_hi)
S["c3_dilution"]=S[f"sens_f{F_PASS}"]>=SENS_PASS
# 4 specificity across 14 TCGA cancer types (tumour means from per-probe summaries)
paad_mean=read_targets("TCGA_PAAD")[C["TCGA-PAAD"][0]].astype(float).mean(axis=1)
S["paad_mean"]=paad_mean
sig_p=S.paad_mean-S.blood_mean
pos=pd.DataFrame(index=S.index)
for c in PAN:
    sm=pd.read_csv(f"{D}/TCGA_{c}_summary.csv.gz",usecols=["probe","tumor_mean"],index_col=0)["tumor_mean"]
    S[f"mean_{c}"]=sm.reindex(S.index)
    pos[c]=(S[f"mean_{c}"]-S.blood_mean)>=PAN_FRAC*sig_p
S["n_other_cancers"]=pos.sum(axis=1)
S["other_cancers"]=pos.apply(lambda r:";".join(r.index[r]),axis=1)
S["c4_specific"]=(sig_p>0)&(S.n_other_cancers<=PAN_MAX)
S["all4"]=S.c1_repro&S.c2_blood&S.c3_dilution&S.c4_specific
S.index.name="probe"
S.to_csv(f"{OUT}/probe_criteria.csv")
# ---- gene level ----
G=gp.merge(S.reset_index(),on="probe",how="left")
def grade(d):
    if not d.measured_both.fillna(False).any(): return "NE"
    if d.all4.fillna(False).any(): return "A"
    if d.c1_repro.fillna(False).any(): return "B"
    return "C"
rows=[]
for g,d in G.groupby("gene"):
    d=d.copy(); rep=d[d.c1_repro.fillna(False).astype(bool)]
    best=None
    if len(rep):
        score=rep[["c2_blood","c3_dilution","c4_specific"]].fillna(False).astype(int).sum(axis=1)*10+rep["sens_f0.1"].fillna(0)
        best=rep.loc[score.idxmax()]
    rows.append(dict(gene=g,n_probes=len(d),n_measured=int(d.measured_both.fillna(False).sum()),
        n_repro=int(d.c1_repro.fillna(False).sum()),n_repro_blood=int((d.c1_repro&d.c2_blood).fillna(False).sum()),
        n_all4=int(d.all4.fillna(False).sum()),grade=grade(d),
        best_probe=None if best is None else best.probe,
        best_auc_TCGA=None if best is None else round(best["auc_TCGA-PAAD"],3),
        best_auc_GSE49149=None if best is None else round(best["auc_GSE49149"],3),
        best_blood_median=None if best is None else round(best.blood_median,3),
        best_sens_f10=None if best is None else round(best["sens_f0.1"],3),
        best_n_other=None if best is None else int(best.n_other_cancers),
        best_other=None if best is None else best.other_cancers,
        best_fails=None if best is None else ",".join(k for k,v in [("blood",best.c2_blood),("dilution",best.c3_dilution),("specificity",best.c4_specific)] if not v),
        hypo_repro=bool(((G.loc[d.index,"delta_TCGA-PAAD"]<=-DELTA)&(G.loc[d.index,"delta_GSE49149"]<=-DELTA)&(G.loc[d.index,"auc_TCGA-PAAD"]<=1-AUC)&(G.loc[d.index,"auc_GSE49149"]<=1-AUC)).any())))
Gd=pd.DataFrame(rows).sort_values(["grade","best_sens_f10"],ascending=[True,False])
Gd.to_csv(f"{OUT}/gene_grades.csv",index=False)
print(Gd.grade.value_counts().sort_index().to_dict())
print(Gd[Gd.grade!="C"].to_string(index=False))
