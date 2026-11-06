import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from samples import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,
    "axes.edgecolor":"#52514e","axes.labelcolor":"#0b0b0b","xtick.color":"#52514e","ytick.color":"#52514e",
    "axes.grid":True,"grid.color":"#e6e5e1","grid.linewidth":0.6,"axes.axisbelow":True,"savefig.dpi":300,"figure.facecolor":"white"})
BLUE,ORANGE,AQUA,GRAY,INK,INK2="#2a78d6","#eb6834","#1baf7a","#b8b7b0","#0b0b0b","#52514e"
GC={"A":BLUE,"B":ORANGE,"C":GRAY}
FIG=os.path.join(ROOT,"results","figures"); os.makedirs(FIG,exist_ok=True)
S=pd.read_csv(f"{OUT}/probe_criteria.csv",index_col=0)
G=pd.read_csv(f"{OUT}/gene_grades.csv"); G=G[G.grade!="NE"]
gp=pd.read_csv(f"{OUT}/gene_probe_map.csv")

# Figure 1: funnel of genes through criteria (gene passes a step if >=1 probe passes all criteria up to that step)
M=gp.merge(S.reset_index(),on="probe",how="left").fillna({"c1_repro":False,"c2_blood":False,"c3_dilution":False,"c4_specific":False})
for c in ["c1_repro","c2_blood","c3_dilution","c4_specific"]: M[c]=M[c].astype(bool)
steps=[("Compiled & mappable",M.gene.nunique()),
       ("Reproducible in tissue",M[M.c1_repro].gene.nunique()),
       ("+ low blood background",M[M.c1_repro&M.c2_blood].gene.nunique()),
       ("+ detectable at 10% dilution",M[M.c1_repro&M.c2_blood&M.c3_dilution].gene.nunique()),
       ("+ not pan-cancer",M[M.c1_repro&M.c2_blood&M.c3_dilution&M.c4_specific].gene.nunique())]
fig,ax=plt.subplots(figsize=(6.4,2.8))
ys=np.arange(len(steps))[::-1]
ax.barh(ys,[v for _,v in steps],color=[GRAY,ORANGE,ORANGE,ORANGE,BLUE],height=0.62)
for y,(l,v) in zip(ys,steps): ax.text(v+1,y,str(v),va="center",color=INK,fontsize=9)
ax.set_yticks(ys,[l for l,_ in steps]); ax.set_xlabel("Number of genes"); ax.grid(axis="y",visible=False)
ax.set_xlim(0,max(v for _,v in steps)*1.12)
fig.tight_layout(); fig.savefig(f"{FIG}/fig1_funnel.png"); plt.close(fig)

# Figure 2: probe-level tissue AUC vs blood background
R=S[S.measured_both]
x=R[["auc_TCGA-PAAD","auc_GSE49149"]].min(axis=1); y=R.blood_median
fig,ax=plt.subplots(figsize=(5.6,4.4))
cat=np.where(R.all4,"All four criteria",np.where(R.c1_repro,"Reproducible in tissue","Not reproducible"))
for lab,col,z in [("Not reproducible",GRAY,1),("Reproducible in tissue",ORANGE,2),("All four criteria",BLUE,3)]:
    m=cat==lab; ax.scatter(x[m],y[m],s=10 if z<3 else 30,c=col,alpha=0.6 if z<3 else 1,edgecolors="white",linewidths=0.4,label=f"{lab} (n={m.sum()})",zorder=z)
ax.axvline(0.70,color=INK2,lw=0.8,ls="--"); ax.axhline(0.10,color=INK2,lw=0.8,ls="--")
ax.set_xlabel("Tumor vs normal AUC (lower of TCGA-PAAD and GSE49149)"); ax.set_ylabel("Median methylation in healthy blood (beta)")
ax.legend(frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.14),ncol=3,fontsize=7.5,markerscale=1.5)
for _,r in R[R.all4].iterrows():
    g=gp.loc[gp.probe==r.name,"gene"].iloc[0]; ax.annotate(g,(min(r["auc_TCGA-PAAD"],r["auc_GSE49149"]),r.blood_median),xytext=(4,6),textcoords="offset points",fontsize=8,color=INK)
fig.tight_layout(); fig.savefig(f"{FIG}/fig2_tissue_vs_blood.png"); plt.close(fig)

# Figure 3: dilution curves (best probe per gene) for genes reproducible + blood-suitable, top 8 by sens at 10%
F=[1.0,0.5,0.2,0.1,0.05,0.02,0.01]
_c=G[G.best_probe.notna()&~G.best_fails.fillna("").str.contains("blood")].sort_values("best_sens_f10",ascending=False)
top=pd.concat([_c[_c.grade=="A"],_c[_c.grade=="B"].head(6)])
fig,ax=plt.subplots(figsize=(6.2,3.8))
labels=[]
for _,r in top.iterrows():
    s=[S.loc[r.best_probe,f"sens_f{f}"] for f in F]
    col=BLUE if r.grade=="A" else ORANGE; lw=2 if r.grade=="A" else 1.2
    ax.plot(F,s,marker="o",ms=4,color=col,lw=lw,alpha=1 if r.grade=="A" else 0.75)
    labels.append([r.gene,s[3]])
labels.sort(key=lambda x:x[1]); ypos=[]
for g,yv in labels:
    yy=max(yv,ypos[-1]+0.045) if ypos else yv; ypos.append(yy)
    ax.annotate(g,(0.1,yv),xytext=(0.075,yy),textcoords="data",fontsize=7.5,color=INK,va="center",
                arrowprops=dict(arrowstyle="-",color="#b8b7b0",lw=0.5))
ax.set_xscale("log"); ax.set_xticks(F,[f"{int(f*100)}%" for f in F]); ax.invert_xaxis()
ax.set_xlabel("Tumor tissue DNA fraction in simulated blood sample"); ax.set_ylabel("Sensitivity at 99% specificity")
ax.axhline(0.5,color=INK2,lw=0.8,ls="--"); ax.axvline(0.1,color=INK2,lw=0.8,ls=":")
from matplotlib.lines import Line2D
ax.legend([Line2D([],[],color=BLUE,lw=2),Line2D([],[],color=ORANGE,lw=1.2)],["Grade A","Grade B"],frameon=False,loc="upper right")
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_dilution.png"); plt.close(fig)

# Figure 4: specificity heatmap: signal over blood in each cancer type for the best probe of reproducible, blood-suitable genes
H=G[G.best_probe.notna()&~G.best_fails.fillna("").str.contains("blood")].sort_values(["grade","best_n_other"])
cols=["paad_mean"]+[f"mean_{c}" for c in PAN]
mat=np.array([[S.loc[p,c]-S.loc[p,"blood_mean"] for c in cols] for p in H.best_probe])
fig,ax=plt.subplots(figsize=(7.2,0.18*len(H)+1.4))
im=ax.imshow(np.clip(mat,0,None),aspect="auto",cmap="Blues",vmin=0,vmax=np.nanpercentile(mat,98))
ax.set_xticks(range(len(cols)),["PAAD"]+PAN,rotation=60,ha="right"); ax.set_yticks(range(len(H)),[f"{g} ({gr})" for g,gr in zip(H.gene,H.grade)],fontsize=7)
ax.grid(False); cb=fig.colorbar(im,ax=ax,fraction=0.03,pad=0.02); cb.set_label("Mean tumor beta minus blood mean",fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_specificity.png"); plt.close(fig)

# Figure 5: study-level reproducibility vs sample size
St=pd.read_csv(f"{OUT}/study_level.csv").dropna(subset=["n_pdac_num","frac_reproducible"])
fig,ax=plt.subplots(figsize=(5.2,3.8))
dg=St.purpose.str.startswith("diagnostic")
ax.scatter(St.n_pdac_num[dg],St.frac_reproducible[dg]*100,s=28,c=BLUE,edgecolors="white",lw=0.5,label="Diagnostic")
ax.scatter(St.n_pdac_num[~dg],St.frac_reproducible[~dg]*100,s=28,c=ORANGE,marker="s",edgecolors="white",lw=0.5,label="Prognostic / staging / monitoring")
ax.set_xscale("log"); ax.set_xlabel("Number of PDAC patients in original study (log scale)"); ax.set_ylabel("Genes reproducible in tissue (%)")
ax.set_xticks([5,10,30,100,300],["5","10","30","100","300"])
ax.legend(frameon=False,fontsize=8,loc="upper center",bbox_to_anchor=(0.5,-0.2),ncol=2)
ax.set_ylim(-5,108)
fig.tight_layout(); fig.savefig(f"{FIG}/fig5_study_level.png"); plt.close(fig)
print("figures written")
