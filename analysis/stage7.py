# exact assay-region matching for Henriksen 2016 MSP amplicons (GRCh38 -> hg19 by chain file)
import sys, os, gzip, bisect, json; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
R=pd.read_csv(os.path.join(LITDIR,"henriksen_regions_hg38.csv"))
need=set(R.chrom)
blocks={}   # chrom -> list of (t_start, t_end, q_chrom, q_start, strand)
with gzip.open(f"{REF}/hg38ToHg19.over.chain.gz","rt") as f:
    cur=None
    for line in f:
        if line.startswith("chain"):
            p=line.split(); tN,tStart=p[2],int(p[5]); qN,qStrand,qStart,qSize=p[7],p[9],int(p[10]),int(p[8])
            cur=(tN,tStart,qN,qStart,qStrand,qSize) if tN in need else None; tpos,qpos=tStart,qStart; continue
        if cur is None or not line.strip(): continue
        p=list(map(int,line.split()))
        size=p[0]
        if cur[4]=="+": blocks.setdefault(cur[0],[]).append((tpos,tpos+size,cur[2],qpos))
        tpos+=size+(p[1] if len(p)>1 else 0); qpos+=size+(p[2] if len(p)>2 else 0)
for c in blocks: blocks[c].sort()
def lift(c,pos1):   # 1-based position
    L=blocks.get(c,[]); x=pos1-1; i=bisect.bisect_right([b[0] for b in L],x)-1
    if i<0: return None
    ts,te,qc,qs=L[i]
    return (qc, qs+(x-ts)+1) if ts<=x<te else None
M=pd.read_csv(f"{REF}/450k_manifest_compact.csv.gz",dtype=str); M["pos"]=M.MAPINFO.astype(float); M["chr"]="chr"+M.CHR
S=pd.read_csv(f"{OUT}/probe_criteria.csv",index_col=0)
ens=json.load(open(f"{REF}/henriksen_regions_grch37.json")) if os.path.exists(f"{REF}/henriksen_regions_grch37.json") else {}
rows=[]
for r in R.itertuples():
    a=lift(r.chrom,r.all_start); b=lift(r.chrom,r.all_end)
    if not a or not b: rows.append(dict(gene=r.gene,lifted=False)); continue
    c,s,e=a[0],min(a[1],b[1]),max(a[1],b[1])
    chk=ens.get(r.gene,{}).get("mappings")
    ens_ok=None if not chk else (chk[0]["mapped"]["start"]==s and chk[0]["mapped"]["end"]==e)
    sub=M[(M.chr==c)]
    inside=sub[(sub.pos>=s)&(sub.pos<=e)].IlmnID.tolist()
    near=sub.assign(d=np.minimum(abs(sub.pos-s),abs(sub.pos-e))).query("pos<@s or pos>@e").nsmallest(3,"d")
    within300=sub[(sub.pos>=s-300)&(sub.pos<=e+300)].IlmnID.tolist()
    def summ(pl):
        pl=[p for p in pl if p in S.index and S.loc[p,"measured_both"]]
        if not pl: return dict(n=0)
        d=S.loc[pl]; best=d.sort_values(["c1_repro","c2_blood","c3_dilution","c4_specific","sens_f0.1"],ascending=False).iloc[0]
        g="A" if d.all4.any() else ("B" if d.c1_repro.any() else "C")
        return dict(n=len(pl),grade=g,best=best.name,auc_min=round(min(best["auc_TCGA-PAAD"],best["auc_GSE49149"]),3),
                    blood_med=round(best.blood_median,3),sens10=round(best["sens_f0.1"],3),n_other=int(best.n_other_cancers))
        
    si=summ(inside); s3=summ(within300)
    rows.append(dict(gene=r.gene,lifted=True,hg19=f"{c}:{s}-{e}",ensembl_agrees=ens_ok,amplicon_bp=e-s+1,
        probes_inside=";".join(inside),n_inside_measured=si.get("n",0),grade_inside=si.get("grade"),
        nearest_outside=";".join(f"{p}({int(d)}bp)" for p,d in zip(near.IlmnID,near.d)),
        n_within300=s3.get("n",0),grade_within300=s3.get("grade"),best_within300=s3.get("best"),
        auc_min=s3.get("auc_min"),blood_med=s3.get("blood_med"),sens10=s3.get("sens10"),n_other=s3.get("n_other")))
T=pd.DataFrame(rows)
G=pd.read_csv(f"{OUT}/gene_grades.csv").set_index("gene"); P=pd.read_csv(f"{OUT}/promoter_grades.csv").set_index("gene")
T["grade_all_probes"]=T.gene.map(G.grade); T["grade_promoter_island"]=T.gene.map(P.grade_promoter)
T.to_csv(f"{OUT}/assay_matched_henriksen.csv",index=False)
pd.set_option("display.width",250)
print(T[["gene","hg19","ensembl_agrees","amplicon_bp","n_inside_measured","grade_inside","n_within300","grade_within300","best_within300","auc_min","blood_med","sens10","n_other","grade_all_probes","grade_promoter_island","nearest_outside"]].to_string(index=False))
