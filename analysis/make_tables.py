import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from samples import *
T=os.path.join(ROOT,"results","tables"); os.makedirs(T,exist_ok=True)
C,blood,clin=cohorts()
t1=[("TCGA-PAAD","UCSC Xena (TCGA)","USA (multi-center)","PDAC tissue",f"{len(C['TCGA-PAAD'][0])} tumors (ductal type), {len(C['TCGA-PAAD'][1])} normal"),
    ("GSE49149","GEO","Australia","PDAC tissue",f"{len(C['GSE49149'][0])} tumors, {len(C['GSE49149'][1])} adjacent non-tumor"),
    ("GSE74071","GEO","Taiwan","PDAC tissue",f"{len(C['GSE74071'][0])} tumors, {len(C['GSE74071'][1])} adjacent normal"),
    ("GSE40279","GEO","USA","Healthy whole blood",f"{len(blood)} adults")]
for c in PAN:
    x=pd.read_csv(f"{D}/TCGA_{c}_samples.csv"); t1.append((f"TCGA-{c}","UCSC Xena (TCGA)","USA (multi-center)","Other cancer tissue",f"{(x.group=='tumor').sum()} tumors"))
pd.DataFrame(t1,columns=["Dataset","Source","Country","Role","Samples used"]).to_csv(f"{T}/Table1_datasets.csv",index=False)
G=pd.read_csv(f"{OUT}/gene_grades.csv")
lit=pd.read_csv(LIT).set_index("gene_symbol")
G["n_studies"]=G.gene.map(lit.n_studies)
t2=G[["gene","grade","n_studies","n_probes","n_repro","best_probe","best_auc_TCGA","best_auc_GSE49149","best_blood_median","best_sens_f10","best_n_other","best_fails"]].copy()
t2.columns=["Gene","Grade","Studies reporting","Probes","Reproducible probes","Representative probe","AUC TCGA-PAAD","AUC GSE49149","Blood median beta","Sensitivity at 10% (99% spec)","Other cancer types positive (of 14)","Criteria not met"]
t2=t2.sort_values(["Grade","Sensitivity at 10% (99% spec)"],ascending=[True,False])
t2.to_csv(f"{T}/Table2_gene_grades.csv",index=False)
with pd.ExcelWriter(f"{T}/Supplementary_Tables.xlsx") as w:
    pd.read_csv(os.path.join(LITDIR,"pdac_markers_studies.csv")).to_excel(w,"S1 studies",index=False)
    pd.read_csv(os.path.join(LITDIR,"screening_log.csv")).to_excel(w,"S2 PubMed screening",index=False)
    t2.to_excel(w,"S3 gene grades",index=False)
    pd.read_csv(f"{OUT}/probe_criteria.csv").merge(pd.read_csv(f"{OUT}/gene_probe_map.csv"),on="probe").to_excel(w,"S4 probe criteria",index=False)
    pd.read_csv(f"{OUT}/robustness.csv").to_excel(w,"S5 robustness",index=False)
    pd.read_csv(f"{OUT}/panel_performance.csv").to_excel(w,"S6 panels",index=False)
    pd.read_csv(f"{OUT}/study_level.csv").to_excel(w,"S7 study level",index=False)
    pd.read_csv(f"{OUT}/promoter_grades.csv").to_excel(w,"S8 promoter island grades",index=False)
    a=pd.read_csv(f"{OUT}/assay_matched_henriksen.csv"); a[~a.gene.isin(["CHFR","TAC1","VIM","MEST_v1"])].to_excel(w,"S9 assay-region matching",index=False)
print("tables written")
