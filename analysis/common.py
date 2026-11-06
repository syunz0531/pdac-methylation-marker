import pandas as pd, numpy as np, gzip, os
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D=os.path.join(ROOT,"extraction","output"); DX=os.path.join(ROOT,"extraction","output_extra")
REF=os.path.join(ROOT,"extraction","refs"); LITDIR=os.path.join(ROOT,"literature")
OUT=os.path.join(ROOT,"results","intermediate"); os.makedirs(OUT,exist_ok=True)
LIT=os.path.join(LITDIR,"pdac_markers_genes.csv")
# gene symbol -> names used in the TCGA legacy probe map
LEGACY={"SEPTIN9":["SEPT9"],"ALKAL1":["FAM150A"],"MARCHF11":["MARCH11"],"CAVIN3":["PRKCDBP"],
        "PYCARD":["PYCARD","TMS1"],"MIR663":["MIR663","MIR663A"],"CCDC181":["CCDC181","C1orf114"],
        "LINC01158":["LINC01158","PANTR1","C2orf85"],"MIR129-2":["MIR129-2"]}
def probemap():
    pm=pd.read_csv(f"{D}/probemap_450k.csv.gz"); pm["genes"]=pm["genes"].fillna("")
    return pm
def gene_probes(pm, genes):
    rows=[]
    tok=pm["genes"].str.replace(";",",").str.split(",")
    lookup={}
    for g in genes:
        for n in LEGACY.get(g,[g]): lookup.setdefault(n.upper(),set()).add(g)
    for p,L,ch,st in zip(pm["probe"],tok,pm["chrom"],pm["start"]):
        for x in L:
            x=x.strip().upper()
            if x in lookup:
                for g in lookup[x]: rows.append((g,p,ch,st))
    return pd.DataFrame(rows,columns=["gene","probe","chrom","start"]).drop_duplicates(["gene","probe"])
def read_targets(name, extra=False):
    f=f"{D}/{name}_targets.csv.gz"
    df=pd.read_csv(f,index_col=0)
    fx=f"{DX}/{name}_targets.csv.gz"
    if os.path.exists(fx):
        dx=pd.read_csv(fx,index_col=0); df=pd.concat([df,dx[~dx.index.isin(df.index)]])
    return df
def tcga_type(bc):
    c=int(bc.split("-")[3][:2]); return "tumor" if c<=9 else ("normal" if c<=19 else "other")
