import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from samples import *
from scipy.stats import mannwhitneyu
from scipy.stats import false_discovery_control
gp=pd.read_csv(f"{OUT}/gene_probe_map.csv")
C,blood,clin=cohorts()
files={"TCGA-PAAD":"TCGA_PAAD","GSE49149":"GSE49149","GSE74071":"GSE74071"}
rows=[]
for coh,(tc,nc) in C.items():
    df=read_targets(files[coh])
    df=df.loc[df.index.isin(gp.probe)]
    T=df[tc].to_numpy(float); N=df[nc].to_numpy(float)
    for i,p in enumerate(df.index):
        t=T[i][~np.isnan(T[i])]; n=N[i][~np.isnan(N[i])]
        if len(t)<3 or len(n)<3: continue
        U,pv=mannwhitneyu(t,n,alternative="two-sided")
        rows.append(dict(cohort=coh,probe=p,n_t=len(t),n_n=len(n),mean_t=t.mean(),mean_n=n.mean(),
                         med_t=np.median(t),med_n=np.median(n),delta=t.mean()-n.mean(),auc=U/(len(t)*len(n)),p=pv))
R=pd.DataFrame(rows)
R["fdr"]=R.groupby("cohort")["p"].transform(lambda x: false_discovery_control(x.to_numpy(),method="bh"))
R.to_csv(f"{OUT}/tissue_probe_stats.csv",index=False)
b=read_targets("GSE40279"); b=b.loc[b.index.isin(gp.probe), blood].astype(float)
B=pd.DataFrame({"blood_n":b.notna().sum(1),"blood_mean":b.mean(1),"blood_median":b.median(1),
                "blood_p95":b.quantile(.95,axis=1),"blood_p99":b.quantile(.99,axis=1),"blood_max":b.max(1)})
B.index.name="probe"; B.to_csv(f"{OUT}/blood_probe_stats.csv")
print(R.groupby("cohort").size(), len(B))
