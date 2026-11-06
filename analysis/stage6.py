# promoter CpG-island restricted grading (proxy for the MSP assay region)
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
M=pd.read_csv(os.path.join(REF,"450k_manifest_compact.csv.gz"),dtype=str)
gp=pd.read_csv(f"{OUT}/gene_probe_map.csv"); S=pd.read_csv(f"{OUT}/probe_criteria.csv",index_col=0)
LEG={g:[g]+v for g,v in LEGACY.items()}
M=M.set_index("IlmnID")
prom=set(["TSS1500","TSS200","5'UTR","1stExon"])
def is_prom_island(probe,gene):
    if probe not in M.index: return False
    r=M.loc[probe]
    names=str(r.UCSC_RefGene_Name).split(";"); groups=str(r.UCSC_RefGene_Group).split(";")
    ok=any(n in LEG.get(gene,[gene]) and g in prom for n,g in zip(names,groups))
    return ok and r.Relation_to_UCSC_CpG_Island=="Island"
gp["prom_island"]=[is_prom_island(p,g) for p,g in zip(gp.probe,gp.gene)]
X=gp.merge(S.reset_index(),on="probe",how="left")
for c in ["c1_repro","c2_blood","c3_dilution","c4_specific","all4","measured_both"]: X[c]=X[c].fillna(False).astype(bool)
G=pd.read_csv(f"{OUT}/gene_grades.csv").set_index("gene")
rows=[]
for g,d in X.groupby("gene"):
    p=d[d.prom_island&d.measured_both]
    gr="NE" if len(p)==0 else ("A" if p.all4.any() else ("B" if p.c1_repro.any() else "C"))
    rows.append(dict(gene=g,grade_all=G.loc[g,"grade"],n_prom_island=len(p),n_prom_repro=int(p.c1_repro.sum()),grade_promoter=gr,
                     prom_best_auc=round(p[["auc_TCGA-PAAD","auc_GSE49149"]].min(axis=1).max(),3) if len(p) else None))
P=pd.DataFrame(rows); P.to_csv(f"{OUT}/promoter_grades.csv",index=False)
print(pd.crosstab(P.grade_all,P.grade_promoter))
print(P[(P.grade_all!=P.grade_promoter)].to_string(index=False))
print(P[P.gene.isin(["CCND2","HOXA1"])].to_string(index=False))
gp.to_csv(f"{OUT}/gene_probe_map.csv",index=False)
