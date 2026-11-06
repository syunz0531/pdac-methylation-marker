# study-level comparison: grades of each study's genes vs its sample size / reported performance
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
from scipy.stats import spearmanr, mannwhitneyu
G=pd.read_csv(f"{OUT}/gene_grades.csv").set_index("gene")
L=pd.read_csv(os.path.join(LITDIR,"pdac_markers_long.csv")); St=pd.read_csv(os.path.join(LITDIR,"pdac_markers_studies.csv"))
L=L[L.gene_symbol.isin(G.index)].drop_duplicates(["study_id","gene_symbol"])
L["grade"]=L.gene_symbol.map(G.grade)
agg=L.groupby("study_id").agg(n_genes=("gene_symbol","nunique"),n_A=("grade",lambda x:(x=="A").sum()),
     n_B=("grade",lambda x:(x=="B").sum()),n_C=("grade",lambda x:(x=="C").sum()))
agg["frac_reproducible"]=(agg.n_A+agg.n_B)/agg.n_genes
St=St.set_index("study_id").join(agg)
St["n_pdac_num"]=pd.to_numeric(St.n_pdac,errors="coerce")
St.to_csv(f"{OUT}/study_level.csv")
d=St.dropna(subset=["n_pdac_num","frac_reproducible"])
rho,p=spearmanr(d.n_pdac_num,d.frac_reproducible)
print("studies",len(d),"spearman rho",round(rho,3),"p",round(p,3))
small=d[d.n_pdac_num<50].frac_reproducible; large=d[d.n_pdac_num>=50].frac_reproducible
print("small<50:",len(small),round(small.median(),2),"large>=50:",len(large),round(large.median(),2),"MWU p",round(mannwhitneyu(small,large).pvalue,3))
print(St[["first_author","year","n_pdac","purpose","n_genes","n_A","n_B","n_C","frac_reproducible"]].to_string())
