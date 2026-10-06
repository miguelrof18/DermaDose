"""Reproducible extraction of the supplied report. Requires Poppler pdftotext.

Usage: python scripts/extract_benmakhlouf.py SOURCE.pdf [--check]
Original decimal strings are retained; no OCR, inference or corrections.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'src/skindose_verifier/data/benmakhlouf_2011'
FIELDS = ['5', '10', '20', '25', '35']


def extract(pdf):
    pages = subprocess.check_output(['pdftotext', '-layout', str(pdf), '-']).decode().split('\f')
    outputs = {}
    def save(name, rows):
        out = io.StringIO(newline='')
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
        outputs[name] = out.getvalue()
    mono, clinical, comparison, mammography = [], [], [], []
    for page_number in range(25, 34):
        row_number = 0
        for line in pages[page_number-1].splitlines():
            cells = line.split()
            if page_number == 25 and len(cells) == 11 and all(re.fullmatch(r'\d+(?:\.\d+)?', c) for c in cells):
                row_number += 1
                for j, side in enumerate(FIELDS):
                    mono.append(dict(row_id=f'T1-E{cells[0]}-S{side}', source_table='1', source_page='25',
                        source_row=str(row_number), energy_kev=cells[0], field_side_cm=side,
                        k_bs=cells[j+1], k_f=cells[j+6]))
            if 27 <= page_number <= 32 and len(cells) == 15 and all(re.fullmatch(r'\d+(?:\.\d+)?', c) for c in cells):
                row_number += 1
                for j, side in enumerate(FIELDS):
                    clinical.append(dict(row_id=f'T3-K{cells[0]}-Al{cells[1]}-Cu{cells[2]}-S{side}',
                        source_table='3', source_page=str(page_number), source_row=str(row_number),
                        kvp=cells[0], inherent_al_mm=cells[1], added_cu_mm=cells[2],
                        hvl1_mm_al=cells[3], h_original=cells[4], field_side_cm=side,
                        k_bs=cells[5+2*j], k_f=cells[6+2*j],
                        quality_flag='h_outside_0_1_in_source' if float(cells[4]) > 1 else ''))
            if page_number == 26 and re.match(r'^\s*[a-f]\s+\d+ cm', line):
                c = re.split(r'\s{2,}', line.strip())
                assert len(c) == 10, c
                comparison.append(dict(row_id='T2-'+c[0], source_table='2', source_page='26', spectrum=c[0],
                    field_original=c[1], kvp=c[2], filtration_original=c[3], k_bs_full_mc=c[4],
                    k_bs_averaged=c[5], bs_ratio_original=c[6], k_f_full_mc=c[7],
                    k_f_averaged=c[8], f_ratio_original=c[9]))
            if page_number == 33 and ('RQA-M' in line or 'IEC 1223-3-2' in line):
                c = re.split(r'\s{2,}', line.strip())
                assert len(c) == 16, c
                for j, side in enumerate(FIELDS):
                    mammography.append(dict(row_id=f'T4-{c[0].replace(" ", "_")}-K{c[2]}-S{side}',
                        source_table='4', source_page='33', beam_quality=c[0], target=c[1], kvp=c[2],
                        inherent_filtration_original=c[3], additional_filtration_original=c[4],
                        hvl1_mm_al=c[5], field_side_cm=side, k_bs=c[6+2*j], k_f=c[7+2*j]))
    assert tuple(map(len, (mono, comparison, clinical, mammography))) == (150, 6, 715, 30)
    save('table1_monoenergetic.csv', mono)
    save('table2_comparison.csv', comparison)
    save('table3_clinical.csv', clinical)
    save('table4_mammography.csv', mammography)
    outputs['source_tables_layout.txt'] = '\f'.join(pages[24:33])
    return outputs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('pdf', type=Path); p.add_argument('--check', action='store_true')
    args = p.parse_args()
    manifest = json.loads((DEST/'metadata.json').read_text()) if (DEST/'metadata.json').exists() else None
    digest = hashlib.sha256(args.pdf.read_bytes()).hexdigest()
    if manifest and digest != manifest['source']['pdf_sha256']:
        raise ValueError('Source PDF SHA-256 differs from the versioned primary source')
    outputs = extract(args.pdf)
    for name, content in outputs.items():
        if args.check:
            assert (DEST/name).read_text() == content, name
        else:
            (DEST/name).write_text(content, encoding='utf-8')
    print(json.dumps({'pdf_sha256': digest, 'files': {k: hashlib.sha256(v.encode()).hexdigest() for k,v in outputs.items()}}, indent=2))

if __name__ == '__main__':
    main()
