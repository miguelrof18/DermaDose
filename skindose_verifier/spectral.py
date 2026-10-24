"""Derived spectral factors (Benmakhlouf Eq.13); published datasets stay immutable."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
from importlib.metadata import version
import math
import numpy as np
from .domain import InputError
from .surface_factors import benmakhlouf_2011
from .factor_tables import RectilinearFactorTable


def integrate_spectrum(energy, fluence, muen_air, field_side_cm, *, dataset=None,
                       bin_width_kev=1.0, max_excluded_kerma_fraction=1e-4):
    """Inputs: bin centres keV, differential fluence/keV, mu_en/rho air cm2/g.

    Returns combined C=<B(E)*f(E)>_K and B=<B(E)>_K; f_effective=C/B.
    f_effective is a bookkeeping factor, NOT the separately averaged Eq.12 ratio.
    """
    dataset = dataset or benmakhlouf_2011()
    e,phi,mu=np.asarray(energy,float),np.asarray(fluence,float),np.asarray(muen_air,float)
    if (e.ndim!=1 or e.size<1 or phi.shape!=e.shape or mu.shape!=e.shape or
        not all(np.isfinite(x).all() for x in (e,phi,mu)) or np.any(e<=0) or
        np.any(np.diff(e)<=0) or np.any(phi<0) or np.any(mu<=0)):
        raise InputError('Espectro/coeficientes inválidos.')
    widths=np.broadcast_to(np.asarray(bin_width_kev,float),e.shape)
    if not np.isfinite(widths).all() or np.any(widths<=0) or not 0<=max_excluded_kerma_fraction<1:
        raise InputError('Anchos de bin o tolerancia inválidos.')
    side=float(field_side_cm)
    table=dataset.tables['monoenergetic'][()]
    energies=table.nodes['energy_kev']
    sides=table.nodes['field_side_cm']
    if not math.isfinite(side) or not sides[0]<=side<=sides[-1]:
        raise InputError('Campo fuera del dominio monoenergético; no se extrapola.')
    w=e*phi*mu*widths
    if not np.isfinite(w).all() or w.sum()<=0:raise InputError('Kerma espectral nulo/inválido.')
    mask=(e>=energies[0])&(e<=energies[-1])
    excluded=float(w[~mask].sum()/w.sum())
    if not np.any(mask) or excluded>max_excluded_kerma_fraction:
        raise InputError(f'Fracción de kerma fuera de tabla={excluded:g}; supera tolerancia explícita.')
    # Interpolate the PRODUCT on original monoenergetic nodes, not product of interpolants.
    combined_rows=[{**r,'k_bs':float(r['k_bs'])*float(r['k_f']),'k_f':1.0} for r in table.rows.values()]
    combined=RectilinearFactorTable(combined_rows,list(table.axes))
    b=table.evaluate(energy_kev=e[mask],field_side_cm=side).k_bs
    c=combined.evaluate(energy_kev=e[mask],field_side_cm=side).k_bs
    weights=w[mask]/w[mask].sum()
    B=float(np.dot(weights,b));C=float(np.dot(weights,c))
    return {'k_bs':B,'k_f':C/B,'combined_factor':C,'excluded_kerma_fraction':excluded,
            'retained_kerma_fraction':1-excluded,'renormalized_after_truncation':excluded>0,
            'energy_domain_kev':[float(energies[0]),float(energies[-1])],
            'conversion_definition':'effective_ratio_C_over_B_not_Eq12',
            'mono_data_sha256':dataset.metadata['files_sha256'][dataset.metadata['tables']['monoenergetic']['file']]}


def generate_spectrum(kvp, total_al_mm, cu_mm, *, anode_angle_deg=12., dk_kev=.5):
    try:import spekpy as sp
    except ImportError as exc:raise InputError('Instale SpekPy: python -m pip install spekpy==2.5.4') from exc
    values=(kvp,total_al_mm,cu_mm,anode_angle_deg,dk_kev)
    if not all(math.isfinite(float(v)) for v in values) or not 5<kvp<=150 or min(total_al_mm,cu_mm)<0 or not 0<anode_angle_deg<90 or not 0<dk_kev<=1:
        raise InputError('Parámetros espectrales inválidos; kVp debe estar entre 5 y 150.')
    s=sp.Spek(kvp=float(kvp),th=float(anode_angle_deg),dk=float(dk_kev),targ='W',mu_data_source='nist')
    if total_al_mm:s.filter('Al',float(total_al_mm))
    if cu_mm:s.filter('Cu',float(cu_mm))
    e,phi=s.get_spectrum(edges=False,flu=True,diff=True)
    mu=s.muen_air_data.get_muen_over_rho_air(e)
    hvl1=float(s.get_hvl1(matl='Al'));hvl2=float(s.get_hvl2(matl='Al'))
    audit={'generator':'SpekPy','version':version('spekpy'),'anode_angle_deg':anode_angle_deg,
           'dk_kev':dk_kev,'kvp':kvp,'total_al_mm':total_al_mm,'cu_mm':cu_mm,
           'hvl1_mm_al':hvl1,'hvl2_mm_al':hvl2,'homogeneity':hvl1/hvl2,
           'muen_air_source':s.muen_air_data.muen_over_rho_air_data_source,
           'physics_model':s.state.model_parameters.physics,
           'spectrum_sha256':hashlib.sha256(np.column_stack((e,phi,mu)).astype('<f8').tobytes()).hexdigest()}
    return e,phi,mu,audit


def benchmark_reference(*, anode_angle_deg=12., dk_kev=.5):
    """Compare representative published nodes; report residuals, never fit silently."""
    dataset=benmakhlouf_2011();cases=[]
    for kvp,al,cu,side in [(50,2.5,0.,5.),(80,3.1,.1,20.),(100,4.1,.1,10.),(120,3.1,.3,25.)]:
        e,p,m,a=generate_spectrum(kvp,al,cu,anode_angle_deg=anode_angle_deg,dk_kev=dk_kev)
        derived=integrate_spectrum(e,p,m,side,bin_width_kev=dk_kev)
        ref=dataset.evaluate('clinical',inherent_al_mm=al,kvp=kvp,added_cu_mm=cu,field_side_cm=side)
        C=float(ref.k_bs*ref.k_f)
        refrow=dataset.tables['clinical'][(al,)].rows[(float(kvp),float(cu),float(side))]
        refhvl=float(refrow['hvl1_mm_al'])
        cases.append({**a,'field_side_cm':side,'reference_B':float(ref.k_bs),'reference_f':float(ref.k_f),
          'reference_combined':C,'reference_hvl1_mm_al':refhvl,
          'hvl1_difference_percent':100*(a['hvl1_mm_al']/refhvl-1),'derived':derived,'combined_difference_percent':100*(derived['combined_factor']/C-1),
          'backscatter_difference_percent':100*(derived['k_bs']/float(ref.k_bs)-1)})
    return {'cases':cases,'source':'Benmakhlouf report Table 3 vs Eq13 spectral approximation',
            'note':'Different spectrum generators and Eq13 vs separate averaged factors can produce residuals. No equivalence or clinical validation implied.'}


def build_dataset(destination, *, kvps, al_families, copper, fields, anode_angle_deg=12.,
                  dk_kev=.5, measured_hvl=None, hvl_relative_tolerance=.05):
    destination=Path(destination)
    # Require a new folder so the published data can never be overwritten.
    if destination.exists():raise InputError('Seleccione una nueva carpeta de dataset; no se sobrescribe.')
    if any(not values or len(set(values))!=len(values) for values in (kvps,al_families,copper,fields)):
        raise InputError('Ejes vacíos/duplicados.')
    rows=[];spectra=[];measurements=measured_hvl or []
    if not 0<hvl_relative_tolerance<1:raise InputError('Tolerancia HVL inválida.')
    for al in sorted(al_families):
        for kvp in sorted(kvps):
            for cu in sorted(copper):
                e,p,m,a=generate_spectrum(kvp,al,cu,anode_angle_deg=anode_angle_deg,dk_kev=dk_kev)
                matches=[x for x in measurements if all(math.isclose(float(x[k]),v,abs_tol=1e-8) for k,v in [('kvp',kvp),('total_al_mm',al),('cu_mm',cu)])]
                if len(matches)>1:raise InputError('Mediciones HVL duplicadas.')
                a['hvl_measurement_status']='not_supplied'
                if matches:
                    measured=float(matches[0]['hvl1_mm_al'])
                    if not math.isfinite(measured) or measured<=0:raise InputError('HVL medida inválida.')
                    residual=(a['hvl1_mm_al']/measured-1)
                    a.update(hvl_measurement_status='compared',measured_hvl1_mm_al=measured,hvl1_relative_difference=residual,
                             hvl_within_tolerance=abs(residual)<=hvl_relative_tolerance)
                    if abs(residual)>hvl_relative_tolerance:raise InputError('Espectro incompatible con HVL medida; no se genera dataset.')
                spectra.append(a)
                for side in sorted(fields):
                    f=integrate_spectrum(e,p,m,side,bin_width_kev=dk_kev)
                    rows.append({'row_id':f'CALC-{len(rows)+1}','inherent_al_mm':al,'kvp':kvp,
                       'added_cu_mm':cu,'field_side_cm':side,'k_bs':f['k_bs'],'k_f':f['k_f']})
                    a['integration_audit']={k:f[k] for k in ('excluded_kerma_fraction','retained_kerma_fraction','renormalized_after_truncation','energy_domain_kev','conversion_definition','mono_data_sha256')}
    unmatched=[x for x in measurements if not any(all(math.isclose(float(x[k]),float(a[k]),abs_tol=1e-8) for k in ('kvp','total_al_mm','cu_mm')) for a in spectra)]
    if unmatched:raise InputError('Hay medidas HVL que no corresponden a los ejes solicitados.')
    import io
    stream=io.StringIO(newline='');writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    data=stream.getvalue().encode()
    meta={'schema_version':1,'dataset_id':'DERIVED-SPEKPY-BENMAKHLOUF-EQ13','version':'1.0',
      'medium':'water','formalism':'B_air_times_surface_water_air_ratio',
      'source':'Calculated with SpekPy + Benmakhlouf KS-ASF-201101-IR Table 1 Eq13, not published clinical nodes',
      'quantities':{'k_bs':'kerma weighted B(E)','k_f':'effective C(Q)/B(Q); not Eq12'},
      'geometry':benmakhlouf_2011().metadata['geometry'],'units':{'energy':'keV','Al_Cu':'mm','field_side':'cm'},
      'tables':{'clinical':{'axes':['kvp','added_cu_mm','field_side_cm'],'groups':['inherent_al_mm'],
           'file':'derived_clinical.csv','rows':len(rows),'source_table':'calculated_Eq13_from_Table1'}},
      'files_sha256':{'derived_clinical.csv':hashlib.sha256(data).hexdigest()},
      'spectra':spectra,'clinical_validation':False,'hvl_comparison_complete':all(a['hvl_measurement_status']=='compared' for a in spectra),
      'hvl_relative_tolerance':hvl_relative_tolerance,
      'limitations':['Spectrum at air surface; table/pad spectral hardening not simulated.',
         'Central-axis water phantom normal incidence; no patient-specific scattering.',
         'Excluded spectral kerma limited to 1e-4 and renormalized explicitly.',
         'No interpolation/extrapolation across Al families; no spectrum reconstruction from a single HVL.']}
    destination.mkdir(parents=True)
    try:
        (destination/'derived_clinical.csv').write_bytes(data)
        (destination/'metadata.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    except Exception:
        import shutil
        shutil.rmtree(destination)
        raise
    return meta


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('benchmark');b.add_argument('--out',type=Path,required=True)
    g=sub.add_parser('build');g.add_argument('--out',type=Path,required=True)
    for flag in ('kvps','al','cu','fields'):g.add_argument('--'+flag,required=True,help='Lista separada por comas')
    g.add_argument('--measured-hvl',type=Path);g.add_argument('--hvl-tolerance',type=float,default=.05)
    for parser in (b,g):parser.add_argument('--anode-angle',type=float,required=True);parser.add_argument('--dk',type=float,default=.5)
    a=p.parse_args()
    try:
        if a.command=='benchmark':
            data=benchmark_reference(anode_angle_deg=a.anode_angle,dk_kev=a.dk)
            with a.out.open('x',encoding='utf8') as f:json.dump(data,f,indent=2)
            print(json.dumps([{k:r[k] for k in ('kvp','total_al_mm','cu_mm','combined_difference_percent')} for r in data['cases']],indent=2))
        else:
            axes={k:[float(v) for v in getattr(a,k).split(',')] for k in ('kvps','al','cu','fields')}
            meta=build_dataset(a.out,kvps=axes['kvps'],al_families=axes['al'],copper=axes['cu'],fields=axes['fields'],
                anode_angle_deg=a.anode_angle,dk_kev=a.dk,
                measured_hvl=json.loads(a.measured_hvl.read_text()) if a.measured_hvl else None,
                hvl_relative_tolerance=a.hvl_tolerance)
            print(json.dumps({'dataset':str(a.out),'rows':meta['tables']['clinical']['rows'],'clinical_validation':False,'hvl_comparison_complete':meta['hvl_comparison_complete']}))
    except (InputError,ValueError,OSError) as exc:p.exit(2,str(exc)+'\n')

if __name__=='__main__':main()
