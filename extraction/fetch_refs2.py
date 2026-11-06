# -*- coding: utf-8 -*-
r"""측정 위치 확인 2단계: 원 논문 좌표(GRCh38)를 칩 주석 좌표(hg19)로 변환하고, Shinjo 2020 보충 표를 받는다. 결과는 refs\ 폴더."""
import json, os, time, urllib.request
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,"refs"); os.makedirs(OUT,exist_ok=True)
UA={"User-Agent":"Mozilla/5.0 (pdac-methylation-study)"}
def get(u,t=60): return urllib.request.urlopen(urllib.request.Request(u,headers=UA),timeout=t).read()
REGIONS=[["ALX4", "chr11", 44309869, 44309983], ["APC", "chr5", 112737725, 112737823], ["BMP3", "chr4", 81031015, 81031129], ["BNC1", "chr15", 83284403, 83284505], ["BRCA1", "chr17", 43125327, 43125441], ["CDKN2A", "chr9", 21974925, 21975058], ["CDKN2B", "chr9", 22009019, 22009138], ["CHFR", "chr12", 132887175, 132887286], ["ESR1", "chr6", 151807700, 151807809], ["EYA2", "chr20", 46894860, 46894964], ["GSTP1", "chr11", 67583649, 67583760], ["HIC1", "chr17", 2056054, 2056196], ["MEST", "chr7", 130486225, 130486355], ["MEST_v1", "chr7", 130492052, 130492182], ["MGMT", "chr10", 129467213, 129467342], ["MLH1", "chr3", 36993529, 36993652], ["NEUROG1", "chr5", 135535860, 135535952], ["NPTX2", "chr7", 98617274, 98617368], ["PENK", "chr8", 56446032, 56446165], ["RARB", "chr3", 25428388, 25428515], ["RASSF1", "chr3", 50340782, 50340897], ["SEPTIN9", "chr17", 77373442, 77373560], ["SFRP1", "chr8", 41309435, 41309584], ["SFRP2", "chr4", 153788922, 153789046], ["SST", "chr3", 187670234, 187670299], ["TAC1", "chr7", 97732043, 97732126], ["TFPI2", "chr7", 93890745, 93890831], ["VIM", "chr10", 17229296, 17229439], ["WNT5A", "chr3", 55487159, 55487311]]
res={}
for g,c,s,e in REGIONS:
    u=f"https://rest.ensembl.org/map/human/GRCh38/{c[3:]}:{s}..{e}:1/GRCh37?content-type=application/json"
    try:
        res[g]=json.loads(get(u)); print("[좌표]",g,"ok",flush=True)
    except Exception as ex:
        res[g]={"error":str(ex)}; print("[좌표]",g,"실패",ex,flush=True)
    time.sleep(0.4)
json.dump(res,open(os.path.join(OUT,"henriksen_regions_grch37.json"),"w"))
for i in range(1,4):
    n=f"Shinjo2020_S{i}"; u=f"https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0233782.s00{i}"
    try:
        d=get(u); ext=".xlsx" if d[:2]==b"PK" else (".pdf" if d[:4]==b"%PDF" else ".bin")
        open(os.path.join(OUT,n+ext),"wb").write(d); print("[보충]",n+ext,len(d),flush=True)
    except Exception as ex: print("[보충]",n,"실패",ex,flush=True)
    time.sleep(1.5)
try:
    d=get("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz",120)
    open(os.path.join(OUT,"hg38ToHg19.over.chain.gz"),"wb").write(d); print("[체인] ok",len(d))
except Exception as ex: print("[체인] 실패",ex)
print("=== 끝. 이 창을 닫아도 됩니다 ===")
