# -*- coding: utf-8 -*-
r"""
췌장암 메틸화 마커 재검증용 공개 데이터 추출 스크립트.

공개 데이터(TCGA via UCSC Xena, GEO series matrix)를 내려받으면서 바로바로 처리해서
(디스크에 원본을 저장하지 않음) 다음만 작은 파일로 남긴다.
  1) 검증 대상 유전자 근처 프로브의 전체 샘플 값   -> *_targets.csv.gz
  2) 모든 프로브의 그룹별 요약 통계(평균/SD/분위수) -> *_summary.csv.gz
  3) 샘플 정보(GEO 메타데이터, TCGA 임상정보)       -> *_samples.csv / clinical

이미 끝난 데이터셋은 건너뛰므로, 중간에 끊겨도 다시 실행하면 이어서 진행된다.
결과는 이 파일과 같은 폴더의 output\ 아래에 생긴다.
"""
import csv, gzip, io, math, os, sys, time, traceback
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)
LOG = open(os.path.join(OUT, "log.txt"), "a", encoding="utf-8")

def log(msg):
    line = time.strftime("%H:%M:%S ") + msg
    print(line, flush=True)
    LOG.write(line + "\n"); LOG.flush()

XENA = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

# ---- 검증 대상 유전자 (문헌 23편, 78개) + 옛 이름/별칭 -------------------------
TARGET_GENES = """BNC1 SFRP1 TFPI2 NPTX2 RASSF1 ADAMTS1 APC BMP3 CDKN2A PENK SFRP2 CCND2 CDKN2B
HOXA1 MEST MLH1 MYOD1 PRKCB SEPTIN9 SPARC SST UCHL1 VHL ZNF781 ADAMTS2 ALKAL1 ALX4 BRCA1
CAVIN3 CCDC181 CD1D CDKN1C CLEC11A DAPK1 DCC EPB41L3 ESR1 FER1L4 GPC3 GRIN2D GSTP1 HIC1
HOXA7 LINC01158 LRFN5 LRRC4 MAPT MARCHF11 MGMT MIR129-2 MIR663 MUC2 NEUROG1 NTRK3 PCDH10
PGK1 PGR PLAU PXDN PYCARD RARB RB1 RYR2 SEMA5A SHISA9 SIX3 SLC9A3 SOCS1 SPSB4 SYK TBR1
THBS1 TRIM73 VWC2 WNT5A ZNF154
SEPT9 FAM150A MARCH11 PRKCDBP SDPR TMS1 MIR663A MIR1292 C2orf85 PANTR1 TCONS_00003256
FLJ35776 C19orf70 C14orf162""".split()
TARGET_GENES = set(g.upper() for g in TARGET_GENES)

# ---- 데이터셋 목록 -------------------------------------------------------------
PAN_CANCER = ["CHOL", "LIHC", "STAD", "ESCA", "COAD", "READ", "LUAD", "LUSC",
              "BRCA", "PRAD", "KIRC", "BLCA", "HNSC", "UCEC"]

GEO_SETS = [
    # (accession, 설명, 그룹 규칙 종류)
    ("GSE49149", "PDAC tumor vs adjacent non-tumor (450K)", "tissue"),
    ("GSE74071", "PDAC vs adjacent normal + CAF + juice cells (450K)", "tissue"),
    ("GSE40279", "whole blood, 656 healthy adults (450K)", "blood"),
    ("GSE35069", "sorted leukocyte subsets, 6 donors (450K)", "celltype"),
]

def geo_url(acc):
    stem = acc[:-3] + "nnn"
    return f"{GEO}{stem}/{acc}/matrix/{acc}_series_matrix.txt.gz"

# ---- 공통 유틸 -----------------------------------------------------------------
def open_stream(url, gz=True, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "pdac-methylation-study/1.0"})
            resp = urllib.request.urlopen(req, timeout=120)
            if gz:
                return io.TextIOWrapper(gzip.GzipFile(fileobj=resp), encoding="utf-8", errors="replace")
            return io.TextIOWrapper(resp, encoding="utf-8", errors="replace")
        except Exception as e:
            last = e
            log(f"  접속 실패({i+1}/{tries}): {e} - 30초 후 재시도")
            time.sleep(30)
    raise last

def done_flag(name):
    return os.path.join(OUT, f".done_{name}")

def is_done(name):
    return os.path.exists(done_flag(name))

def mark_done(name):
    open(done_flag(name), "w").close()

def to_float(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except Exception:
        return None

def stats(vals):
    v = sorted(x for x in vals if x is not None)
    n = len(v)
    if n == 0:
        return [0, "", "", "", "", "", ""]
    mean = sum(v) / n
    sd = math.sqrt(sum((x - mean) ** 2 for x in v) / (n - 1)) if n > 1 else 0.0
    q = lambda p: v[min(n - 1, int(round(p * (n - 1))))]
    return [n, round(mean, 4), round(sd, 4), round(q(0.05), 4), round(q(0.5), 4),
            round(q(0.95), 4), round(v[-1], 4)]

STAT_COLS = ["n", "mean", "sd", "p05", "median", "p95", "max"]

# ---- 1. 프로브 -> 유전자 지도 --------------------------------------------------
PROBE_GENE = {}

def load_probemap():
    path = os.path.join(OUT, "probemap_450k.csv.gz")
    if not os.path.exists(path):
        log("[프로브 지도] 내려받는 중")
        url = XENA + "probeMap%2FilluminaMethyl450_hg19_GPL16304_TCGAlegacy"
        src = open_stream(url, gz=False)
        tmp = path + ".part"
        with gzip.open(tmp, "wt", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["probe", "genes", "chrom", "start"])
            for line in src:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) >= 4:
                    w.writerow([p[0], p[1], p[2], p[3]])
        os.replace(tmp, path)
    with gzip.open(path, "rt", encoding="utf-8") as f:
        r = csv.reader(f); next(r)
        for probe, genes, chrom, start in r:
            PROBE_GENE[probe] = genes
    n_t = sum(1 for p in PROBE_GENE if is_target(p))
    log(f"[프로브 지도] 프로브 {len(PROBE_GENE):,}개, 대상 유전자 프로브 {n_t:,}개")

def is_target(probe):
    g = PROBE_GENE.get(probe, "")
    if not g or g == ".":
        return False
    return any(x.strip().upper() in TARGET_GENES for x in g.replace(";", ",").split(","))

# ---- 2. TCGA (UCSC Xena) ------------------------------------------------------
def tcga_group(barcode):
    # TCGA-XX-XXXX-01A : 01~09 종양, 10~19 정상
    try:
        code = int(barcode.split("-")[3][:2])
    except Exception:
        return "other"
    if code <= 9:
        return "tumor"
    if 10 <= code <= 19:
        return "normal"
    return "other"

def do_tcga(project, summary=True):
    name = f"TCGA_{project}"
    if is_done(name):
        log(f"[{name}] 이미 완료 - 건너뜀"); return
    log(f"[{name}] 메틸화 데이터 처리 시작")
    url = XENA + f"TCGA.{project}.sampleMap%2FHumanMethylation450.gz"
    src = open_stream(url)
    header = src.readline().rstrip("\n").split("\t")
    samples = header[1:]
    groups = [tcga_group(s) for s in samples]
    gnames = sorted(set(groups))
    idx = {g: [i for i, x in enumerate(groups) if x == g] for g in gnames}
    log(f"  샘플 {len(samples)}개: " + ", ".join(f"{g} {len(idx[g])}" for g in gnames))

    t_path = os.path.join(OUT, f"{name}_targets.csv.gz")
    s_path = os.path.join(OUT, f"{name}_summary.csv.gz")
    tf = gzip.open(t_path + ".part", "wt", encoding="utf-8", newline="")
    tw = csv.writer(tf); tw.writerow(["probe"] + samples)
    sw = None
    if summary:
        sf = gzip.open(s_path + ".part", "wt", encoding="utf-8", newline="")
        sw = csv.writer(sf)
        sw.writerow(["probe"] + [f"{g}_{c}" for g in gnames for c in STAT_COLS])
    n = 0
    for line in src:
        p = line.rstrip("\n").split("\t")
        probe = p[0]
        if not probe.startswith("cg"):
            continue
        vals = [to_float(x) for x in p[1:]]
        if is_target(probe):
            tw.writerow([probe] + ["" if v is None else round(v, 4) for v in vals])
        if sw:
            row = [probe]
            for g in gnames:
                row += stats([vals[i] for i in idx[g]])
            sw.writerow(row)
        n += 1
        if n % 50000 == 0:
            log(f"  {n:,} 프로브 처리")
    tf.close(); os.replace(t_path + ".part", t_path)
    if sw:
        sf.close(); os.replace(s_path + ".part", s_path)
    with open(os.path.join(OUT, f"{name}_samples.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f); w.writerow(["sample", "group"])
        for s, g in zip(samples, groups):
            w.writerow([s, g])
    mark_done(name)
    log(f"[{name}] 완료 ({n:,} 프로브)")

def do_tcga_clinical(project):
    name = f"TCGA_{project}_clinical"
    if is_done(name):
        return
    log(f"[{name}] 임상정보 내려받는 중")
    url = XENA + f"TCGA.{project}.sampleMap%2F{project}_clinicalMatrix"
    src = open_stream(url, gz=False)
    with open(os.path.join(OUT, f"{name}.tsv"), "w", encoding="utf-8") as f:
        for line in src:
            f.write(line)
    mark_done(name)

# ---- 3. GEO series matrix -----------------------------------------------------
def geo_group(meta_text, rule):
    t = meta_text.lower()
    if rule == "blood":
        return "blood"
    if rule == "celltype":
        for key, g in [("whole blood", "whole_blood"), ("pbmc", "PBMC"), ("granulocyte", "granulocytes"),
                       ("neutrophil", "neutrophils"), ("eosinophil", "eosinophils"),
                       ("cd4", "CD4T"), ("cd8", "CD8T"), ("cd19", "Bcell"), ("cd14", "monocytes"),
                       ("cd56", "NK")]:
            if key in t:
                return g
        return "other"
    for key in ["fibroblast", "caf", "juice", "circulating", "cell line"]:
        if key in t:
            return "other"
    for key in ["non-tumor", "non tumor", "nontumor", "non-tumour", "adjacent", "normal", "benign"]:
        if key in t:
            return "normal"
    for key in ["tumor", "tumour", "cancer", "carcinoma", "pdac", "pdca"]:
        if key in t:
            return "tumor"
    return "other"

def do_geo(acc, desc, rule):
    name = acc
    if is_done(name):
        log(f"[{name}] 이미 완료 - 건너뜀"); return
    log(f"[{name}] {desc} 처리 시작")
    src = open_stream(geo_url(acc))
    meta = {}
    samples = None
    for line in src:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_"):
            p = line.rstrip("\n").split("\t")
            key = p[0][1:]
            vals = [x.strip('"') for x in p[1:]]
            if key in meta:
                k = 2
                while f"{key}_{k}" in meta:
                    k += 1
                key = f"{key}_{k}"
            meta[key] = vals
    header = src.readline().rstrip("\n").split("\t")
    samples = [x.strip('"') for x in header[1:]]
    nS = len(samples)
    text_keys = [k for k in meta if any(s in k for s in ("title", "source_name", "characteristics", "description"))]
    groups = []
    for i in range(nS):
        txt = " | ".join(meta[k][i] for k in text_keys if i < len(meta[k]))
        groups.append(geo_group(txt, rule))
    gnames = sorted(set(groups))
    idx = {g: [i for i, x in enumerate(groups) if x == g] for g in gnames}
    log(f"  샘플 {nS}개: " + ", ".join(f"{g} {len(idx[g])}" for g in gnames))

    with open(os.path.join(OUT, f"{name}_samples.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        keys = [k for k in meta if len(meta[k]) == nS]
        w.writerow(["sample", "auto_group"] + keys)
        for i, s in enumerate(samples):
            w.writerow([s, groups[i]] + [meta[k][i] for k in keys])

    t_path = os.path.join(OUT, f"{name}_targets.csv.gz")
    s_path = os.path.join(OUT, f"{name}_summary.csv.gz")
    tf = gzip.open(t_path + ".part", "wt", encoding="utf-8", newline="")
    tw = csv.writer(tf); tw.writerow(["probe"] + samples)
    sf = gzip.open(s_path + ".part", "wt", encoding="utf-8", newline="")
    sw = csv.writer(sf); sw.writerow(["probe"] + [f"{g}_{c}" for g in gnames for c in STAT_COLS])
    n = 0
    for line in src:
        if line.startswith("!series_matrix_table_end"):
            break
        p = line.rstrip("\n").split("\t")
        probe = p[0].strip('"')
        if not probe.startswith("cg"):
            continue
        vals = [to_float(x.strip('"')) for x in p[1:]]
        if is_target(probe):
            tw.writerow([probe] + ["" if v is None else round(v, 4) for v in vals])
        row = [probe]
        for g in gnames:
            row += stats([vals[i] for i in idx[g] if i < len(vals)])
        sw.writerow(row)
        n += 1
        if n % 50000 == 0:
            log(f"  {n:,} 프로브 처리")
    tf.close(); sf.close()
    if n < 1000:
        log(f"  경고: 프로브가 {n}개뿐입니다. 이 데이터셋은 series matrix에 값이 없을 수 있어요.")
    os.replace(t_path + ".part", t_path); os.replace(s_path + ".part", s_path)
    mark_done(name)
    log(f"[{name}] 완료 ({n:,} 프로브)")

# ---- 실행 순서 -----------------------------------------------------------------
def run(step, *args):
    try:
        step(*args)
    except Exception as e:
        log(f"!! {step.__name__}{args} 실패: {e}")
        LOG.write(traceback.format_exc())

if __name__ == "__main__":
    log("=== 시작 ===")
    load_probemap()
    run(do_tcga_clinical, "PAAD")
    run(do_tcga, "PAAD")
    for acc, desc, rule in GEO_SETS:
        run(do_geo, acc, desc, rule)
    for proj in PAN_CANCER:
        run(do_tcga, proj, True)
    pending = [n for n in (["TCGA_PAAD"] + [a for a, _, _ in GEO_SETS] + [f"TCGA_{p}" for p in PAN_CANCER])
               if not is_done(n)]
    if pending:
        log("=== 끝났지만 실패한 항목이 있어요: " + ", ".join(pending) + " (다시 실행하면 이어서 시도) ===")
    else:
        log("=== 전부 완료! 이 창을 닫아도 됩니다 ===")
