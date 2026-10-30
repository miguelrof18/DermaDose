"""Philips Azurion profile and read-only compatibility inspection.

No patient identifiers are included in the report. The original DICOM is never edited.
"""
from dataclasses import replace
from collections import Counter
import argparse
import json
import math
from pathlib import Path

from .domain import InputError
from .rdsr import read_rdsr

SOURCES = {
    'ifu': 'https://www.documents.philips.com/assets/Instruction%20for%20Use/20201222/c122581e1172411b8c0eac9900949a77.pdf',
    'clarityiq': 'https://www.philips.com/c-dam/b2bhc/master/landing-pages/queen-alexandra/philips-igt-azurion-clarityiq-whitepaper.pdf',
}
SPECIFICATION = {
    'manufacturer': 'Philips', 'model': 'Azurion',
    'minimum_tube_assembly_total_al_equivalent_mm': 2.5,
    'tube_inherent_al_equivalent_mm': 0.8,
    'included_filter_plate_al_mm': 1.7,
    'anode_angle_deg': 11.0,
    'assembly_max_kvp': 125.0,
    'minimum_is_not_measured_total_filtration': True,
    'spectrabeam_added_filters_mm': [{'al': 1.0, 'cu': v} for v in (0.1, 0.4, 0.9)],
    'scope': 'Published family specifications; verify tube, software release and beam limiting device at installation.',
    'sources': SOURCES,
}


def prepare_events(events, equipment):
    """Preserve measured event geometry. Research overrides require explicit opt-in."""
    profile = equipment.raw.get('azurion_profile', {})
    if not profile:
        return tuple(events)
    if profile.get('name') != 'philips_azurion':
        raise InputError('azurion_profile.name debe ser philips_azurion.')
    if profile.get('tube_model') == 'MRC 200 0407':
        from .mrc200 import validate_configuration
        validate_configuration(events, equipment)
    research = profile.get('allow_research_assumptions') is True
    if not research and not all(profile.get(k) is True for k in
                                ('fixed_filtration_characterized', 'geometry_characterized', 'calibration_characterized')):
        raise InputError('Azurion: faltan caracterización de filtración fija, geometría o calibración. Consulte python -m skindose_verifier.azurion.')
    result = []
    first_table = next((e.table_mm for e in events if e.table_mm is not None), None)
    if profile.get('reference_table_from_first_event'):
        if not research or first_table is None:
            raise InputError('Referencia de mesa desde primer evento solo permitida como hipótesis de investigación explícita.')
        equipment.geometry['table_origin_mm'] = list(first_table)
    for event in events:
        updates = {}
        quality = dict(event.quality)
        if profile.get('reference_table_from_first_event'):
            quality['table_reference'] = 'research_reference_first_event_not_clinical_registration'
        for name, fallback in [('sid_mm', 'source_to_detector_mm'), ('sod_mm', 'source_to_isocenter_mm'), ('srp_mm', 'source_to_irp_mm')]:
            value = getattr(event, name)
            if value is None or not math.isfinite(value) or value <= 0:
                if research and profile.get('invalid_distance_policy') == 'configured_research':
                    updates[name] = float(equipment.geometry_for(event.plane)[fallback])
                    quality[name] = 'research_replacement_invalid_or_missing_event_distance'
                else:
                    raise InputError(f'Azurion: evento con {name} ausente/no positivo. No se infiere una distancia clínica.')
        updates['quality'] = quality
        result.append(replace(event, **updates))
    return tuple(result)


def inspect_azurion(path, fixed_al_mm=None):
    report = read_rdsr(path)
    if (report.manufacturer.strip().lower(), report.model.strip().lower()) != ('philips', 'azurion'):
        raise InputError('El RDSR no identifica Philips/Azurion.')
    filters = Counter()
    missing = Counter()
    families = Counter()
    for e in report.events:
        al = sum(p['selected_mm'] for p in e.original.get('filter_stack', [])
                 if (p['material_scheme'], p['material_code']) == ('SRT', 'C-120F9'))
        filters[(e.quality['filtration_mm_cu'], al, e.filtration_mm_cu)] += 1
        for name in ('sid_mm', 'sod_mm', 'srp_mm', 'field_width_mm', 'field_height_mm'):
            if getattr(e, name) is None or getattr(e, name) <= 0:
                missing[name] += 1
        if e.table_mm is None:
            missing['table_mm'] += 1
        if fixed_al_mm is not None and e.filtration_mm_cu is not None:
            families[(round(fixed_al_mm + al, 8), e.filtration_mm_cu)] += 1
    out = {'device': {'manufacturer': report.manufacturer, 'model': report.model},
           'event_count': len(report.events), 'published_specification': SPECIFICATION,
           'filter_groups': [{'status':k[0], 'added_al_mm':k[1], 'added_cu_mm':k[2], 'events':v} for k,v in filters.items()],
           'missing_or_nonpositive_geometry': dict(missing),
           'physical_calculation_ready': False,
           'note': 'Reading compatibility is not physical validation. Fixed filtration, coordinate mapping and calibrations require installation data.'}
    if fixed_al_mm is not None:
        from .surface_factors import benmakhlouf_2011
        dataset = benmakhlouf_2011()
        checks=[]
        for (al,cu),count in families.items():
            try:
                # Only assesses the filtration coordinate; full kVp and local field checks are separate.
                table = dataset.tables['clinical'].get((al,))
                if table is None: raise InputError('Familia Al no tabulada')
                cu_nodes=table.nodes['added_cu_mm']
                if cu < cu_nodes[0] or cu > cu_nodes[-1]:
                    raise InputError(f'Cu fuera de dominio [{cu_nodes[0]}, {cu_nodes[-1]}]')
                checks.append({'total_al_mm':al,'cu_mm':cu,'events':count,'filtration_supported':True})
            except InputError as exc:
                checks.append({'total_al_mm':al,'cu_mm':cu,'events':count,'filtration_supported':False,'reason':str(exc)})
        out['hypothetical_fixed_al_mm']=fixed_al_mm
        out['candidate_table_families']=checks
    return out


def main():
    parser=argparse.ArgumentParser(description='Inspección Azurion sin cálculo ni modificación del RDSR')
    parser.add_argument('rdsr',type=Path)
    parser.add_argument('--fixed-al-mm',type=float,help='Hipótesis explícita para estudiar familias; no modifica el equipo')
    args=parser.parse_args()
    try: print(json.dumps(inspect_azurion(args.rdsr,args.fixed_al_mm),indent=2,ensure_ascii=False))
    except InputError as exc: parser.exit(2,str(exc)+'\n')

if __name__=='__main__':main()
