# Reassessment of reported blood-based DNA methylation markers for pancreatic cancer

Code, compiled marker list, and results for the study:

> Lee S. *Reassessment of Reported Blood-Based DNA Methylation Markers for Pancreatic Cancer Using Public Data.* (manuscript in preparation)

## Summary

DNA methylation markers in circulating cell-free DNA have been proposed for the early detection of pancreatic ductal adenocarcinoma (PDAC), but most were reported in small single cohorts. This repository compiles 90 genes from 37 blood-based studies and reassesses each in public Illumina 450K data against four criteria:

1. **Reproducibility** – hypermethylated in PDAC tissue in TCGA-PAAD and GSE49149 (delta beta >= 0.10, AUC >= 0.70, FDR < 0.05), with direction checked in GSE74071
2. **Blood suitability** – median beta <= 0.10 and 95th percentile <= 0.20 in healthy whole blood (GSE40279, n = 656)
3. **Detection limit** – sensitivity >= 50% at 99% specificity when tumor tissue DNA is diluted to 10% in blood (in silico)
4. **Specificity** – comparable signal in no more than 4 of 14 other TCGA cancer types

Of 86 evaluable genes, 58 were reproducible in tissue, and two (HOXA1 and CCND2) met all four criteria. Results were checked against promoter CpG island probes and, for 25 genes, the exact published assay regions.

## Repository layout

```
extraction/     Scripts that download and stream-process public data (Python standard library only)
analysis/       Analysis pipeline (stage0-stage7), figures, and tables
literature/     Compiled studies and genes, PubMed search log and screening decisions, published primer coordinates
results/        Figures, tables, supplementary tables, and intermediate probe/gene-level results
```

## Reproducing the analysis

Requirements: Python 3.10 or later. The extraction scripts use only the standard library; the analysis needs the packages in `requirements.txt`.

```
pip install -r requirements.txt

# 1. Download and extract data (several hours; about 10 GB is streamed, about 300 MB is kept)
python extraction/extract.py          # TCGA (15 cancer types) and GEO datasets -> extraction/output/
python extraction/extract_extra.py    # sample-level values for genes added from the PubMed search
python extraction/pubmed_search.py    # PubMed search (results are also provided in literature/)
python extraction/fetch_refs.py       # Illumina 450K manifest and published article files -> extraction/refs/
python extraction/fetch_refs2.py      # coordinate conversion check and liftover chain file

# 2. Run the analysis
cd analysis
python stage0.py   # map genes to array probes
python stage1.py   # tissue and blood probe statistics
python stage2.py   # four criteria and gene grades
python stage3.py   # marker panels
python stage4.py   # study-level comparison
python stage5.py   # threshold robustness
python stage6.py   # promoter CpG island grading
python stage7.py   # exact assay-region matching (Henriksen et al. 2016 primers)
python figures.py
python make_tables.py
```

On Windows, the `.bat` files in `extraction/` run the corresponding scripts by double-clicking. Random processes use fixed seeds, and re-running the pipeline reproduces the files in `results/` exactly.

## Data sources

| Dataset | Source | Use |
|---|---|---|
| TCGA-PAAD | UCSC Xena | PDAC tissue (153 ductal-type tumors, 10 normal) |
| GSE49149 | GEO | PDAC tissue (167 tumors, 29 adjacent non-tumor) |
| GSE74071 | GEO | PDAC tissue (12 tumors, 8 adjacent normal) |
| GSE40279 | GEO | Healthy whole blood (656 adults) |
| TCGA (14 other cancer types) | UCSC Xena | Specificity (5,211 tumors) |

All data are public and were used as provided by the original repositories. Raw data are not redistributed here; the extraction scripts download them from the original sources.

## License

Code is released under the MIT License. The compiled tables summarize published studies, which should be cited directly when used.
