"""Traceable extraction of the supplied QC workbook; conditional spectral research fit."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from .spectral import generate_spectrum,integrate_spectrum,build_dataset

def extract_workbook(path):
    import openpyxl
    w=openpyxl.load_workbook(path,data_only=True)
    s=w['Calidad,tiempo y rendimiento'];rows=[]
    for r in range(31,39):
        value=lambda c:s[f'{c}{r}'].value
        vals=[value(c) for c in ['A','B','D','I','J']]
        if any(not isinstance(x,(int,float)) or not np.isfinite(x) or x<=0 for x in vals):
            raise ValueError(f'Medida inválida en fila {r}')
        nominal,measured,hvl,mas,kerma=vals
        rows.append(dict(row=r,nominal_kvp=nominal,measured_kvp=measured,hvl1_mm_al=hvl,
          yield_gy_per_mas=kerma/mas/1000,source_cells=[f'{c}{r}' for c in ['A','B','D','I','J']]))
    return {'source_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'sheet':s.title,
      'rows':rows,'detector_distance_cm':s['C11'].value,'table_present':s['C22'].value,
      'added_filtration_known':False,'installation_model':'Azurion 7 M20',
      'document_tube':'MRC 200+ GS 0407','document_anode_angle_deg':11,
      'document_model':'Azurion 7 C20 FlexArm',
      'issues':['Filtration cells E20/I20 blank; no automatic zero assumption.',
        'Calidad y rendimiento nominal kVp misaligned with raw log; excluded.',
        'DE015 columns have dose/rate unit ambiguity and PMMA backscatter context; CF not activated.',
        'No measured table transmission supplied; synthetic factor unchanged.',
        'Product document is C20 FlexArm; tube-family agreement does not prove installation identity.']}

def conditional_fit(extracted):
    groups=[]
    for nominal in sorted({r['nominal_kvp'] for r in extracted['rows']}):
        rr=[r for r in extracted['rows'] if r['nominal_kvp']==nominal]
        groups.append({'nominal_kvp':nominal,'measured_kvp':float(np.mean([r['measured_kvp'] for r in rr])),
          'hvl1_mm_al':float(np.mean([r['hvl1_mm_al'] for r in rr])),
          'yield_gy_per_mas':float(np.mean([r['yield_gy_per_mas'] for r in rr])),
          'n':len(rr),'role':'held_out' if nominal==125 else 'fit'})
    def loss(al):
        return np.mean([(generate_spectrum(g['measured_kvp'],al,0,anode_angle_deg=11)[3]['hvl1_mm_al']/g['hvl1_mm_al']-1)**2 for g in groups if g['role']=='fit'])
    opt=minimize_scalar(loss,bounds=(1.,8.),method='bounded',options={'xatol':.001})
    if not opt.success:raise ValueError('Ajuste espectral fallido')
    for g in groups:
        e,p,m,a=generate_spectrum(g['measured_kvp'],opt.x,0,anode_angle_deg=11)
        g['calculated_hvl1_mm_al']=a['hvl1_mm_al'];g['hvl_residual_percent']=100*(a['hvl1_mm_al']/g['hvl1_mm_al']-1)
        g['C_fit_10cm']=integrate_spectrum(e,p,m,10,bin_width_kev=.5)['combined_factor']
        for angle,al,key in [(12,3.1,'C_old_10cm'),(11,3.1,'C_angle11_10cm')]:
            e,p,m,_=generate_spectrum(g['measured_kvp'],al,0,anode_angle_deg=angle)
            g[key]=integrate_spectrum(e,p,m,10,bin_width_kev=.5)['combined_factor']
    return {'effective_al_mm':float(opt.x),'anode_angle_deg':11,'groups':groups,
      'objective':'unweighted mean squared relative HVL1 residual, one weight per nominal quality',
      'fit_assumptions':['Added Cu=0 assumed, not documented','All baseline filtration represented by equivalent Al',
       'Angle 11 degrees from related product document, pending installation confirmation'],
      'clinical_validation':False,'spectral_identifiability_established':False}

def main():
    p=argparse.ArgumentParser();p.add_argument('workbook',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--conditional-fit-assume-cu-zero',action='store_true');p.add_argument('--build-dataset',action='store_true');a=p.parse_args()
    if a.out.exists():p.error('Use a new output directory')
    if a.build_dataset and not a.conditional_fit_assume_cu_zero:p.error('Dataset requires explicit conditional fit opt-in')
    x=extract_workbook(a.workbook);a.out.mkdir(parents=True)
    (a.out/'qc_extracted.json').write_text(json.dumps(x,indent=2,ensure_ascii=False))
    if a.conditional_fit_assume_cu_zero:
        f=conditional_fit(x);(a.out/'conditional_fit.json').write_text(json.dumps(f,indent=2,ensure_ascii=False))
        if a.build_dataset:
            al=f['effective_al_mm'];folder=a.out/'derived_qc_conditional'
            m=build_dataset(folder,kvps=[46,50,60,70,80,90,100,110,120,125,130,140,150],al_families=[al,al+1],copper=[0,.1,.2,.3,.4,.6,.9],fields=[5,10,20,25,35],anode_angle_deg=11)
            m.update(dataset_id='AZURION-QC-CONDITIONAL-EQ13',interpolation_quantity='combined_C',
              qc_conditional_fit=f,qc_source=x,installation_characterization_complete=False,added_filtration_confirmed=False)
            (folder/'metadata.json').write_text(json.dumps(m,indent=2,ensure_ascii=False))
    print(a.out)
if __name__=='__main__':main()
