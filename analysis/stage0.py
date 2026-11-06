# map compiled marker genes to 450K probes (UCSC Xena probe annotation)
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
lit=pd.read_csv(LIT)
genes=[g for g in lit.gene_symbol if g not in ("AK055957","GH05J042948")]   # non-gene identifiers
gp=gene_probes(probemap(),genes)
gp.to_csv(f"{OUT}/gene_probe_map.csv",index=False)
print(len(genes),"genes;",len(gp),"gene-probe pairs;",gp.probe.nunique(),"unique probes")
