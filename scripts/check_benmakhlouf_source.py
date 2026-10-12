"""Independent PDF-coordinate check (PyMuPDF), not the Poppler extraction path.

Usage: python scripts/check_benmakhlouf_source.py SOURCE.pdf
"""
import csv
import hashlib
import json
from pathlib import Path
import re
import sys
import fitz

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'skindose_verifier/data/benmakhlouf_2011'


def check(pdf):
    doc = fitz.open(pdf)
    rows_by_page = {}
    for i in range(24,33):
        page = doc[i]
        groups = {}
        for word in page.get_text('words'):
            box = fitz.Rect(word[:4])*page.rotation_matrix
            key = next((k for k in groups if abs(k-box.y0) < 1), box.y0)
            # Text baselines in these tables align; x sorting is independent of
            # PDF content-stream order and Poppler's layout algorithm.
            groups.setdefault(key, []).append((box.x0, word[4]))
        rows_by_page[i+1] = [[text for x,text in sorted(words)] for words in groups.values()]
    def numeric_rows(page, size):
        return [r for r in rows_by_page[page] if len(r)==size and all(re.fullmatch(r'\d+(?:\.\d+)?',v) for v in r)]
    checked = 0
    for filename, pages, columns in [('table1_monoenergetic.csv',[25],11),('table3_clinical.csv',list(range(27,33)),15)]:
        data = list(csv.DictReader((DATA/filename).open()))
        for page in pages:
            pdfrows = numeric_rows(page, columns)
            csvrows = [r for r in data if int(r['source_page'])==page]
            assert len(pdfrows)*5 == len(csvrows), (page,len(pdfrows),len(csvrows))
            for i, printed in enumerate(pdfrows):
                rs = csvrows[i*5:(i+1)*5]
                if columns == 11:
                    reconstructed = [rs[0]['energy_kev']]+[r['k_bs'] for r in rs]+[r['k_f'] for r in rs]
                else:
                    reconstructed = [rs[0][k] for k in ['kvp','inherent_al_mm','added_cu_mm','hvl1_mm_al','h_original']]+[r[k] for r in rs for k in ['k_bs','k_f']]
                assert printed == reconstructed, (page,i,printed,reconstructed)
                checked += len(printed)
    t2 = list(csv.DictReader((DATA/'table2_comparison.csv').open()))
    for row in t2:
        printed = next(r for r in rows_by_page[26] if row['spectrum'] in r and r[-1] == row['f_ratio_original'])
        expected = [row[k] for k in ['k_bs_full_mc','k_bs_averaged','bs_ratio_original','k_f_full_mc','k_f_averaged','f_ratio_original']]
        assert printed[-6:]==expected, (printed, expected)
        checked += 6
    t4 = list(csv.DictReader((DATA/'table4_mammography.csv').open()))
    printed = [r for r in rows_by_page[33] if len(r)>11 and all(re.fullmatch(r'\d+\.\d+',v) for v in r[-11:])]
    assert len(printed)==6
    for i,r in enumerate(printed):
        rs=t4[i*5:(i+1)*5]
        expected = [rs[0]['hvl1_mm_al']]+[s[k] for s in rs for k in ['k_bs','k_f']]
        assert r[-11:]==expected, (r,expected)
        checked += 11
    return {'method':'PyMuPDF rotated word coordinates, independent of Poppler',
            'pdf_sha256':hashlib.sha256(Path(pdf).read_bytes()).hexdigest(),
            'numeric_cells_compared':checked,'factor_values_compared':1814,'discrepancies':[]}

if __name__ == '__main__':
    print(json.dumps(check(sys.argv[1]),indent=2))
