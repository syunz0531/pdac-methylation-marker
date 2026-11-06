# -*- coding: utf-8 -*-
r"""
측정 위치 확인용 자료 내려받기 (결과는 refs\ 폴더):
  1) Illumina 450K 공식 주석(프로브별 유전자 위치, CpG island 관계) - NCBI GEO
  2) 주요 마커 유전자 프로모터 주변의 hg19 염기서열 - UCSC Genome Browser API
  3) 프라이머 서열이 실린 공개 논문 페이지와 보충 파일 원본
"""
import csv, gzip, io, json, os, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "refs"); os.makedirs(OUT, exist_ok=True)
UA = {"User-Agent": "Mozilla/5.0 (pdac-methylation-study)"}

def log(m): print(time.strftime("%H:%M:%S ") + m, flush=True)

def get(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout)

# ---- 1. 450K manifest (compact) ----
MAN = os.path.join(OUT, "450k_manifest_compact.csv.gz")
URLS = ["https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13534/suppl/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz"]
if not os.path.exists(MAN):
    for url in URLS:
        try:
            log("[주석] 내려받는 중: " + url)
            src = io.TextIOWrapper(gzip.GzipFile(fileobj=get(url, 300)), encoding="utf-8", errors="replace")
            header = None; n = 0
            with gzip.open(MAN + ".part", "wt", encoding="utf-8", newline="") as f:
                w = csv.writer(f)
                for row in csv.reader(src):
                    if header is None:
                        if row and row[0] == "IlmnID":
                            header = row
                            idx = [header.index(c) for c in ["IlmnID", "CHR", "MAPINFO", "Strand",
                                   "UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name",
                                   "Relation_to_UCSC_CpG_Island"]]
                            w.writerow([header[i] for i in idx])
                        continue
                    if row and row[0].startswith("[Controls]"):
                        break
                    if len(row) > max(idx) and row[0].startswith("c"):
                        w.writerow([row[i] for i in idx]); n += 1
            os.replace(MAN + ".part", MAN)
            log(f"[주석] 완료: 프로브 {n:,}개")
            break
        except Exception as e:
            log(f"[주석] 실패: {e}")

# ---- 2. promoter sequences for key genes ----
GENES = """BNC1 ADAMTS1 SFRP1 SFRP2 TFPI2 NPTX2 SST BMP3 APC RASSF1 CDKN2A CDKN2B MEST HOXA1 ADAMTS2
PCDH10 SEMA5A SPSB4 CCND2 HOXA9 TWIST1 WT1 RPRM PENK SEPT9 MIR129-2 NEUROG1 ALX4 WNT5A HIC1 GSTP1
MLH1 RUNX3 HOXD8 POU4F1 BRCA1 BRCA2 MUC2 UCHL1 SPARC CD1D GRIN2D PRKCB ZNF781 CLEC11A""".split()
SEQ = os.path.join(OUT, "promoter_sequences_hg19.json")
if os.path.exists(MAN) and not os.path.exists(SEQ):
    prom = {}
    with gzip.open(MAN, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            names = r["UCSC_RefGene_Name"].split(";"); groups = r["UCSC_RefGene_Group"].split(";")
            for nm, gr in zip(names, groups):
                if nm in GENES and gr in ("TSS1500", "TSS200", "5'UTR", "1stExon"):
                    p = int(r["MAPINFO"]); c = "chr" + r["CHR"]
                    lo, hi, cc = prom.get(nm, (p, p, c)); prom[nm] = (min(lo, p), max(hi, p), c)
    seqs = {}
    for g, (lo, hi, c) in prom.items():
        start, end = max(0, lo - 3000), hi + 3000
        url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom={c};start={start};end={end}"
        try:
            d = json.load(get(url))
            seqs[g] = {"chrom": c, "start": start, "end": end, "dna": d["dna"].upper()}
            log(f"[서열] {g} {c}:{start}-{end} ({end-start:,} bp)")
        except Exception as e:
            log(f"[서열] {g} 실패: {e}")
        time.sleep(0.4)
    json.dump(seqs, open(SEQ, "w"))

# ---- 3. open-access article pages / supplements with primer tables ----
DOCS = {
    "Henriksen2016_main.html": "https://link.springer.com/article/10.1186/s13148-016-0286-2",
    "Henriksen2016_ESM1.docx": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13148-016-0286-2/MediaObjects/13148_2016_286_MOESM1_ESM.docx",
    "Eissa2019_main.html": "https://link.springer.com/article/10.1186/s13148-019-0650-0",
    "Shinjo2020_main.html": "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0233782",
    "Kim2024_main.html": "https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2024.1395473/full",
    "Manoochehri2020_main.html": "https://febs.onlinelibrary.wiley.com/doi/full/10.1002/1878-0261.12684",
    "Yi2013_main.html": "https://aacrjournals.org/clincancerres/article/19/23/6544/78255",
    "Ying2021_main.html": "https://academic.oup.com/gpb/article/19/6/890/7229543",
}
for name, url in DOCS.items():
    path = os.path.join(OUT, name)
    if os.path.exists(path):
        continue
    try:
        data = get(url, 60).read()
        open(path, "wb").write(data); log(f"[논문] {name} ({len(data):,} bytes)")
    except Exception as e:
        log(f"[논문] {name} 실패: {e}")
    time.sleep(1.5)
log("=== 끝. 이 창을 닫아도 됩니다 ===")
