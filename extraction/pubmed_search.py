# -*- coding: utf-8 -*-
r"""
PubMed 검색 스크립트: 검색식을 고정해서 실행하고, 결과 논문 전체(PMID, 제목, 학술지, 연도, 초록)를
output\pubmed_results.csv 로 저장한다. 검색일과 검색식, 건수는 output\pubmed_search_log.txt 에 기록된다.
"""
import csv, os, time, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

QUERY = (
    '("pancreatic cancer"[tiab] OR "pancreatic cancers"[tiab] OR "pancreatic adenocarcinoma"[tiab] '
    'OR "pancreatic ductal adenocarcinoma"[tiab] OR PDAC[tiab]) '
    'AND (cfDNA[tiab] OR "cell-free DNA"[tiab] OR "cell free DNA"[tiab] OR "circulating tumor DNA"[tiab] '
    'OR ctDNA[tiab] OR plasma[tiab] OR serum[tiab] OR "liquid biopsy"[tiab] OR blood[tiab]) '
    'AND (methylation[tiab] OR methylated[tiab] OR hypermethylation[tiab] OR methylome[tiab]) '
    'AND (diagnos*[tiab] OR detection[tiab] OR biomarker*[tiab] OR marker*[tiab] OR screening[tiab]) '
    'AND english[la]'
)

def get(url, data=None):
    for i in range(4):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "pdac-methylation-study/1.0"})
            return urllib.request.urlopen(req, timeout=60).read()
        except Exception as e:
            print(f"  접속 실패({i+1}/4): {e}")
            time.sleep(5)
    raise RuntimeError("PubMed 접속 실패")

def text(el):
    return "".join(el.itertext()).strip() if el is not None else ""

def main():
    print("PubMed 검색 중...")
    q = urllib.parse.urlencode({"db": "pubmed", "term": QUERY, "retmax": 5000,
                                "retmode": "xml", "tool": "pdac_methylation_study"})
    root = ET.fromstring(get(EU + "esearch.fcgi?" + q))
    count = int(root.findtext("Count"))
    pmids = [e.text for e in root.findall("./IdList/Id")]
    print(f"검색 결과 {count}건")

    rows = []
    for i in range(0, len(pmids), 200):
        batch = pmids[i:i + 200]
        data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "xml",
                                       "tool": "pdac_methylation_study"}).encode()
        art_root = ET.fromstring(get(EU + "efetch.fcgi", data=data))
        for art in art_root.findall(".//PubmedArticle"):
            mc = art.find("MedlineCitation")
            a = mc.find("Article")
            year = (mc.findtext(".//JournalIssue/PubDate/Year") or
                    (mc.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4])
            authors = a.findall(".//AuthorList/Author")
            first = ""
            if authors:
                first = (authors[0].findtext("LastName") or authors[0].findtext("CollectiveName") or "")
                ini = authors[0].findtext("Initials")
                if ini:
                    first += " " + ini
            doi = ""
            for aid in art.findall(".//ArticleIdList/ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = aid.text or ""
            pubtypes = ";".join(text(p) for p in a.findall(".//PublicationTypeList/PublicationType"))
            abstract = " ".join(text(t) for t in a.findall(".//Abstract/AbstractText"))
            rows.append([mc.findtext("PMID"), year, first, text(a.find("ArticleTitle")),
                         text(a.find("Journal/Title")), doi, pubtypes, abstract])
        print(f"  {min(i + 200, len(pmids))}/{len(pmids)} 건 정보 받음")
        time.sleep(0.5)

    with open(os.path.join(OUT, "pubmed_results.csv"), "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pmid", "year", "first_author", "title", "journal", "doi", "publication_types", "abstract"])
        w.writerows(rows)
    with open(os.path.join(OUT, "pubmed_search_log.txt"), "w", encoding="utf-8") as f:
        f.write(f"search_date: {time.strftime('%Y-%m-%d %H:%M %Z')}\n")
        f.write(f"database: PubMed (NCBI E-utilities)\n")
        f.write(f"query: {QUERY}\n")
        f.write(f"count: {count}\nretrieved: {len(rows)}\n")
    print(f"완료: {len(rows)}건 저장. 이 창을 닫아도 됩니다.")

if __name__ == "__main__":
    main()
