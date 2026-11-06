# -*- coding: utf-8 -*-
r"""
추가 추출: PubMed 검색으로 새로 들어간 유전자 12개의 샘플별 값을
췌장암 조직 3개 + 정상 혈액 데이터에서만 다시 뽑는다. 결과는 output_extra\ 에 저장된다.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import extract as ex

ex.load_probemap()                      # output\probemap_450k.csv.gz 재사용
ex.TARGET_GENES = set("BRCA2 EYA2 HOXA9 HOXD8 LAX1 NR3C1 POU4F1 RCAN3 RPRM RUNX3 TWIST1 WT1".split())
ex.OUT = os.path.join(HERE, "output_extra")
os.makedirs(ex.OUT, exist_ok=True)

ex.log("=== 추가 추출 시작 ===")
ex.run(ex.do_tcga, "PAAD", False)
for acc, desc, rule in ex.GEO_SETS[:3]:
    ex.run(ex.do_geo, acc, desc, rule)
done = [n for n in ["TCGA_PAAD", "GSE49149", "GSE74071", "GSE40279"] if ex.is_done(n)]
ex.log(f"=== 추가 추출 끝: {len(done)}/4 완료 ===")
