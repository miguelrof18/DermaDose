"""Repeatable model sensitivity, not accuracy or uncertainty assessment."""
import csv,json
from pathlib import Path
from skindose_verifier.spectral import generate_spectrum,integrate_spectrum
root=Path(__file__).resolve().parents[1]
rows=[]
for kvp in (50,80,100,125):
 for cu in (0,.1,.4,.9):
  for side in (5,10,20,35):
   out={'kvp':kvp,'cu_mm':cu,'field_side_cm':side}
   for name,angle,al in [('historical',12,3.1),('angle11',11,3.1),('assembly',11,2.5),('collimator',11,2.6)]:
    # Added Al accompanies Cu only in this explicitly declared comparison scenario.
    e,p,m,a=generate_spectrum(kvp,al+(1 if cu else 0),cu,anode_angle_deg=angle)
    f=integrate_spectrum(e,p,m,side,bin_width_kev=.5)
    out[name+'_C']=f['k_bs']*f['k_f']
    out[name+'_hvl_mm_al']=a['hvl1_mm_al']
   for name in ('angle11','assembly','collimator'):
    out[name+'_change_percent']=100*(out[name+'_C']/out['historical_C']-1)
   rows.append(out)
f=root/'docs/mrc200/model_sensitivity.csv'
with f.open('w',newline='') as stream:
 w=csv.DictWriter(stream,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
summary={name:{'min':min(r[name+'_change_percent'] for r in rows),'max':max(r[name+'_change_percent'] for r in rows)} for name in ('angle11','assembly','collimator')}
(root/'docs/mrc200/model_sensitivity_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary))
