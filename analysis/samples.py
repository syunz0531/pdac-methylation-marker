import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
def cohorts():
    """returns dict cohort -> (tumor_cols, normal_cols) plus blood cols and pan-cancer tumor cols"""
    C={}
    # TCGA PAAD: ductal-type primary tumours vs solid tissue normal
    clin=pd.read_csv(f"{D}/TCGA_PAAD_clinical.tsv",sep="\t")
    duct=set(clin.loc[clin.histological_type=="Pancreas-Adenocarcinoma Ductal Type","sampleID"])
    s=pd.read_csv(f"{D}/TCGA_PAAD_samples.csv")
    tum=[x for x,g in zip(s["sample"],s["group"]) if g=="tumor" and x.split("-")[3][:2]=="01" and x[:15] in duct]
    nor=[x for x,g in zip(s["sample"],s["group"]) if g=="normal"]
    C["TCGA-PAAD"]=(tum,nor)
    s=pd.read_csv(f"{D}/GSE49149_samples.csv")
    C["GSE49149"]=(s.loc[s.Sample_source_name_ch1=="Tumor","sample"].tolist(), s.loc[s.Sample_source_name_ch1=="Adjacent","sample"].tolist())
    s=pd.read_csv(f"{D}/GSE74071_samples.csv")
    C["GSE74071"]=(s.loc[s.Sample_source_name_ch1=="pancreatic cancer","sample"].tolist(), s.loc[s.Sample_source_name_ch1=="adjacent normal","sample"].tolist())
    blood=pd.read_csv(f"{D}/GSE40279_samples.csv")["sample"].tolist()
    return C, blood, clin
PAN=["CHOL","LIHC","STAD","ESCA","COAD","READ","LUAD","LUSC","BRCA","PRAD","KIRC","BLCA","HNSC","UCEC"]
