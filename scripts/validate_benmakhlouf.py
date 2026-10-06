"""Independent Decimal calculations vs provider and complete DermaDose dose chain.

Run from repository: python scripts/validate_benmakhlouf.py
No use of the production interpolation routine for expected values.
"""
from dataclasses import replace
from decimal import Decimal as D
import csv
import json
from pathlib import Path
import tempfile

import numpy as np
from skindose_verifier.config import Equipment
from skindose_verifier.engine import calculate
from skindose_verifier.rdsr import read_rdsr
from skindose_verifier.synthetic import make_synthetic_rdsr
from skindose_verifier.surface_factors import benmakhlouf_2011

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'src/skindose_verifier/data/benmakhlouf_2011'


def independent(rows, query):
    """Successive 1D interpolation using Decimal; a separate implementation."""
    axes=list(query)
    def recurse(candidates, level, column):
        if level==len(axes):
            assert len(candidates)==1
            return D(candidates[0][column])
        axis=axes[level]; x=D(str(query[axis]))
        nodes=sorted({D(r[axis]) for r in candidates})
        low=max(v for v in nodes if v<=x)
        high=min(v for v in nodes if v>=x)
        a=recurse([r for r in candidates if D(r[axis])==low],level+1,column)
        if low==high: return a
        b=recurse([r for r in candidates if D(r[axis])==high],level+1,column)
        return a+(b-a)*(x-low)/(high-low)
    return [recurse(rows,0,k) for k in ('k_bs','k_f')]


def validate():
    clinical=list(csv.DictReader((DATA/'table3_clinical.csv').open()))
    mono=list(csv.DictReader((DATA/'table1_monoenergetic.csv').open()))
    cases=[('clinical',dict(kvp=50,added_cu_mm=0,field_side_cm=5)),
           ('clinical',dict(kvp=80,added_cu_mm=.1,field_side_cm=20)),
           ('clinical',dict(kvp=150,added_cu_mm=.9,field_side_cm=35)),
           ('clinical',dict(kvp=85,added_cu_mm=.15,field_side_cm=15)),
           ('clinical',dict(kvp=73,added_cu_mm=.24,field_side_cm=12)),
           ('monoenergetic',dict(energy_kev=40,field_side_cm=10)),
           ('monoenergetic',dict(energy_kev=42.5,field_side_cm=15))]
    provider=benmakhlouf_2011()
    eq=Equipment.load(ROOT/'examples/equipo_benmakhlouf_investigacion.json')
    with tempfile.TemporaryDirectory() as tmp:
        ev=read_rdsr(make_synthetic_rdsr(Path(tmp)/'test.dcm')).events[0]
    results=[]
    for kind,q in cases:
        rows=[r for r in clinical if D(r['inherent_al_mm'])==D('3.1')] if kind=='clinical' else mono
        bs,f=independent(rows,q)
        actual=provider.evaluate(kind,**q,**({'inherent_al_mm':3.1} if kind=='clinical' else {}))
        expected=[bs,f,bs*f]
        got=[float(actual.k_bs),float(actual.k_f),float(actual.k_bs*actual.k_f)]
        errors=[abs(x-float(y)) for x,y in zip(got,expected)]
        assert max(errors)<1e-14
        entry={'kind':kind,'query':q,'inherent_al_mm':3.1 if kind=='clinical' else None,
               'independent_decimal':dict(zip(['k_bs','k_f','product'],map(str,expected))),
               'dermadose':dict(zip(['k_bs','k_f','product'],got)),
               'absolute_errors':errors,'status':actual.audit['status']}
        if kind=='clinical':
            side=q['field_side_cm']
            event=replace(ev,kvp=q['kvp'],filtration_mm_cu=q['added_cu_mm'],
                          field_width_mm=side*10*1100/600,field_height_mm=side*10*1100/600)
            calculation=calculate([event],eq,'skin_dose',detailed=False)
            mask=calculation.total_gy>0
            # Independent point-source distance calculation; CF=1, T=0.8.
            p=calculation.surface.points_mm[mask]
            r2=p[:,0]**2+p[:,1]**2+(p[:,2]+700)**2
            dose=event.ka_r_gy*550**2/r2*.8*float(bs*f)
            error=float(np.max(np.abs(calculation.total_gy[mask]-dose)))
            assert error<1e-14
            entry['dose_max_absolute_error_gy']=error
            entry['dose_points_checked']=int(mask.sum())
        results.append(entry)
    comparison=[]
    for r in csv.DictReader((DATA/'table2_comparison.csv').open()):
        bmc,fmc=D(r['k_bs_full_mc']),D(r['k_f_full_mc'])
        bav,fav=D(r['k_bs_averaged']),D(r['k_f_averaged'])
        comparison.append({'spectrum':r['spectrum'],
            'bs_relative_difference_percent':float((bav/bmc-1)*100),
            'f_relative_difference_percent':float((fav/fmc-1)*100),
            'product_relative_difference_percent':float((bav*fav/(bmc*fmc)-1)*100)})
    return {'validation_type':'numerical and physical-formalism verification; not clinical validation',
            'cases':results,'source_table2_method_comparison':comparison,
            'max_factor_absolute_error':max(max(r['absolute_errors']) for r in results),
            'max_dose_absolute_error_gy':max(r.get('dose_max_absolute_error_gy',0) for r in results),
            'unexplained_implementation_discrepancies':[]}

if __name__=='__main__':
    print(json.dumps(validate(),indent=2))
