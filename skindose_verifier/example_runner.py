"""Ejemplos con RDSR, equipo JSON y STL locales; entradas inmutables.

Coloca siempre la malla con la misma cadena (docs/sistemas_de_referencia.md):

    TC (C, LPS) --M_T<-C--> tablero (T) --t_k = t_0 + S (m_k - m_0)--> isocentro (I)

Cada eslabón se toma medido o explícito si existe y, si no, se estima con una
hipótesis declarada que queda en configuracion_efectiva.json y en los resultados:

* M_T<-C: registro_tc.json (DermaDose Surface o marcas) o, si falta, estimación
  automática (postura del RDSR, punto diana del TC, colchoneta y distancia a la
  cabecera). El XCAT usa su anclaje posterior de hombros.
* t_0, m_0: marcador medido del perfil de mesa o, si falta, punto diana en el
  isocentro en el evento de referencia.
* S: perfil de mesa; una lectura sin caracterizar detiene el cálculo si cambia.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

import numpy as np


class ExampleError(ValueError):
    pass


def read_json(path):
    data = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if not isinstance(data, dict):
        raise ExampleError(f'El JSON debe contener un objeto: {path}')
    return data


def choose(paths, label, requested=None, interactive=False, input_fn=input):
    paths = sorted((p.resolve() for p in paths), key=lambda p: p.name.casefold())
    if requested is not None:
        target = Path(requested).resolve()
        if target not in paths:
            raise ExampleError(f'{label}: archivo no disponible: {requested}')
        return target
    if len(paths) == 1:
        return paths[0]
    if not paths:
        raise ExampleError(f'No se encuentra {label} en el ejemplo.')
    if not interactive:
        raise ExampleError(f'Hay varios archivos de {label}: ' + ', '.join(p.name for p in paths))
    print(f'Seleccione {label}:')
    for i, path in enumerate(paths, 1):
        print(f'  {i}: {path.name}')
    answer = input_fn('Número: ').strip()
    if not answer.isdigit() or not 1 <= int(answer) <= len(paths):
        raise ExampleError('Selección inválida.')
    return paths[int(answer)-1]


def prepare_example(folder, *, config=None, rdsr=None, mesh=None,
                    interactive=False, input_fn=input,
                    data_directory=None, table_calibration=None, ct_registration=None, pad_mm=None):
    """Select the inputs of a case folder; placement is resolved in ``place_patient``."""
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise ExampleError(f'No existe la carpeta: {folder}')
    resources = Path(data_directory) if data_directory else Path(__file__).resolve().parent/'data'
    entries = [p for p in folder.iterdir() if p.is_file()]
    def local(value):
        return None if value is None else folder/Path(value)
    configs = []
    for path in entries:
        if path.suffix.lower() == '.json':
            try:
                data = read_json(path)
                if all(k in data for k in ('version','manufacturer','model','geometry','factors')):
                    configs.append(path)
            except (ValueError, OSError):
                if config and path.resolve() == local(config).resolve():
                    raise
    if config is not None:
        explicit = local(config).resolve()
        if explicit not in configs and explicit.is_file():
            configs.append(explicit)
    config_path = choose(configs, 'configuración de equipo', local(config), interactive, input_fn)
    data = deepcopy(read_json(config_path))
    if ('reference_table_from_first_event' in data.get('azurion_profile', {})
            and 'table_calibration' not in data.get('geometry', {}) and table_calibration is None):
        # Azurion profiles before 0.7.0 placed the table from the first event.
        raise ExampleError(f'{config_path.name} es un perfil Azurion anterior a DermaDose 0.7.0 (referencia de mesa '
                           'desde el primer evento). Use equipo_azurion.json o config/azurion_benmakhlouf_spekpy.json '
                           'del repositorio y retire la copia antigua de la carpeta del caso.')
    dicom = choose([p for p in entries if p.suffix.lower()=='.dcm'], 'RDSR (.dcm)', local(rdsr), interactive, input_fn)
    candidates = [p for p in entries if p.suffix.lower()=='.stl']
    if interactive and mesh is None:
        answer = input_fn('Superficie: 1 = maniquí XCAT, 2 = malla propia de la carpeta: ').strip()
        if answer not in ('1', '2'):
            raise ExampleError('Seleccione 1 o 2.')
        mesh = 'xcat' if answer == '1' else 'patient'
    use_xcat = mesh == 'xcat'
    if use_xcat:
        if ct_registration is not None:
            raise ExampleError('El registro TC pertenece a la malla propia.')
        selected = resources/'xcat_wholebody_body.stl'
    else:
        selected = choose(candidates, 'malla STL del paciente', local(mesh) if mesh not in (None, 'patient') else None, interactive, input_fn)
    selected = selected.resolve()
    if not selected.is_file():
        raise ExampleError(f'No existe la malla: {selected}')
    previous = data.get('patient', {})
    # The registration defines units, axes and origin of the surface (CT LPS, mm).
    data['patient'] = {'type': 'mesh', 'mesh_path': str(selected), 'phantom_name': 'XCAT50' if use_xcat else selected.stem,
                       'resolution_mm': previous.get('resolution_mm', 5), 'thickness_mm': previous.get('thickness_mm', 220),
                       **{k: previous[k] for k in ('mesh_max_edge_mm', 'mesh_refinement_roi_mm', 'mesh_max_faces') if k in previous}}
    data.pop('phantom', None)  # legacy blocks must not replace the selected STL
    geometry = data['geometry']
    if table_calibration is not None:
        geometry['table_calibration'] = str(Path(table_calibration).resolve())
    elif 'table_calibration' not in geometry and (folder/'mesa_azurion_hemodinamica_1.json').is_file():
        geometry['table_calibration'] = str(folder/'mesa_azurion_hemodinamica_1.json')
    registration_file = None
    if not use_xcat:
        if ct_registration is not None:
            registration_file = Path(ct_registration).resolve()
        elif (folder/'registro_tc.json').is_file():
            registration_file = folder/'registro_tc.json'
    digest = sha256(selected.read_bytes()).hexdigest()
    sf = data.get('surface_factors', {})
    dataset = sf.get('dataset_directory')
    if dataset:
        directory = Path(dataset)
        absolute = directory.is_absolute()
        if not absolute:
            directory = config_path.parent/directory
        if not (directory/'metadata.json').is_file():
            parts = str(dataset).replace('\\','/').split('/')
            known = {'derived_mrc200_nominal','derived_azurion_hvl','derived_azurion_spekpy','derived_azurion_qc_conditional'}
            if not absolute and len(parts)>=3 and parts[-3:-1]==['skindose_verifier','data'] and parts[-1] in known:
                directory = resources/parts[-1]
            if not (directory/'metadata.json').is_file():
                raise ExampleError(f'Falta metadata.json de las tablas espectrales: {directory}')
        sf['dataset_directory'] = str(directory.resolve())
    audit = dict(source_config=str(config_path), source_config_sha256=sha256(config_path.read_bytes()).hexdigest(),
                 rdsr=str(dicom), mesh=str(selected), mesh_sha256=digest,
                 mesh_choice='xcat' if use_xcat else 'patient_stl',
                 registration_file=None if registration_file is None else str(registration_file),
                 pad_mm=None if pad_mm is None else float(pad_mm),
                 dataset_directory=sf.get('dataset_directory'))
    return data, audit


def _target_regions(dicom_path):
    """Target Region (DCM 123014) of the RDSR; used only to choose a target point."""
    import pydicom
    from .rdsr import descendants, concept
    values = set()
    for node in descendants(pydicom.dcmread(dicom_path, stop_before_pixels=True)):
        if concept(node) == ("DCM", "123014"):
            for code in getattr(node, "ConceptCodeSequence", []):
                values.add(str(getattr(code, "CodeMeaning", "")).strip())
    return sorted(v for v in values if v)


def _mesh_vertices_lps(mesh_path, metadata):
    from trimesh import load_mesh
    from .placement import mesh_vertices_lps
    mesh = load_mesh(mesh_path, process=False)
    if metadata.get('length_unit') != 'mm':
        raise ExampleError('La malla debe documentar length_unit = mm.')
    return mesh_vertices_lps(np.asarray(mesh.vertices), metadata.get('coordinate_system'))


def place_patient(data, audit, events, *, interactive=False, input_fn=input, pad_mm=None,
                  shoulder_distance_mm=None, reference_event='kerma_weighted', allow_estimates=None):
    """Resolve M_T<-C and the room reference (t_0, m_0) into ``data``; return the placement audit."""
    from .placement import (CHEST_BELOW_SHOULDERS_MM, estimate_ct_registration, estimate_room_reference,
                            pose_rotation, select_target)
    from .reference_frames import _settings, configure_table, resolve_registration_pose, room_reference_state
    base = Path(audit['source_config']).parent
    geometry = data['geometry']
    if 'table_calibration' not in geometry:
        raise ExampleError('La colocación TC → tablero → isocentro requiere el perfil de la mesa '
                           '(geometry.table_calibration o --mesa).')
    table_geometry = deepcopy(geometry)
    configure_table(table_geometry, base)
    table = table_geometry['table_calibration']
    geometry['table_calibration'] = deepcopy(table)
    if allow_estimates is None:
        allow_estimates = data.get('azurion_profile', {}).get('allow_research_assumptions') is True
    mesh_path = Path(audit['mesh'])
    metadata = read_json(mesh_path.with_suffix('.json'))
    readings = table.get('radimetrics', {}).get('other_readings_mm', {})
    estimates = []

    def require_estimates(what):
        if not allow_estimates:
            raise ExampleError(f'Falta {what} y la configuración no admite estimaciones de investigación '
                               '(azurion_profile.allow_research_assumptions o --permitir-estimaciones).')

    def mattress():
        nonlocal pad_mm
        source = 'explicit'
        if pad_mm is None and interactive:
            default = readings.get('anchura_del_colchon')
            hint = f' [Intro = {0.6*default:g} mm estimados]' if default else ''
            answer = input_fn(f'Espesor de colchoneta bajo carga en mm (0 sin colchoneta){hint}: ').strip()
            pad_mm = answer.replace(',', '.') if answer else None
        if pad_mm is None and readings.get('anchura_del_colchon') is not None:
            pad_mm, source = 0.6*float(readings['anchura_del_colchon']), 'estimated_60_percent_of_unloaded_mattress'
        if pad_mm is None or not math.isfinite(float(pad_mm)) or float(pad_mm) < 0:
            raise ExampleError('Indique el espesor de colchoneta bajo carga con --espesor-colchoneta-mm.')
        return float(pad_mm), source

    def shoulders():
        if shoulder_distance_mm is not None:
            return float(shoulder_distance_mm), 'explicit'
        if readings.get('shoulder_to_head') is not None:
            return float(readings['shoulder_to_head']), 'table_profile_typical_value'
        raise ExampleError('Indique la distancia cabecera–hombros con --distancia-hombros-mm.')

    cache = {}

    def vertices_lps():
        if 'vertices' not in cache:  # large STL: load only when an estimate needs it
            cache['vertices'] = _mesh_vertices_lps(mesh_path, metadata)
        return cache['vertices']

    regions = _target_regions(audit['rdsr'])
    rotation, pose = pose_rotation(events)
    # 1. M_T<-C ------------------------------------------------------------------
    if audit['mesh_choice'] == 'patient_stl' and audit.get('registration_file'):
        registration = _settings(audit['registration_file'], base, 'ct_registration')
        registration = resolve_registration_pose(registration, events, interactive=interactive, input_fn=input_fn)
        placement_status = 'explicit'
    elif audit['mesh_choice'] == 'xcat':
        require_estimates('el registro del XCAT en el tablero')
        pad, pad_source = mattress()
        distance, distance_source = shoulders()
        anchor = np.asarray(metadata['posterior_shoulder_midline_mm'], dtype=float)
        registration = {'mesh_sha256': audit['mesh_sha256'], 'anchor_ct_lps_mm': anchor.tolist(),
                        'anchor_table_mm': [0., -pad, -distance], 'rotation_table_from_ct': rotation.tolist(),
                        'anchor_definition': 'posterior_shoulder_midline', 'status': 'estimated_xcat_shoulder_anchor',
                        'estimate': {'method': 'xcat_posterior_shoulder_anchor', 'mattress_mm': pad,
                                     'mattress_source': pad_source, 'shoulder_distance_mm': distance,
                                     'shoulder_distance_source': distance_source, 'intervention_pose': pose,
                                     'hypotheses': ['Fantoma XCAT de referencia, no la anatomía del paciente.',
                                                    'Línea media centrada en el tablero.']}}
        placement_status = registration['status']
        estimates.append('M_T<-C (XCAT)')
    else:
        require_estimates('registro_tc.json')
        pad, pad_source = mattress()
        distance, distance_source = shoulders()
        target, info = select_target(vertices_lps(), metadata, regions)
        registration = estimate_ct_registration(
            vertices_lps(), target, rotation, mattress_mm=pad,
            target_distance_from_head_mm=distance + CHEST_BELOW_SHOULDERS_MM,
            mesh_sha256=audit['mesh_sha256'], target_info=info, pose_info=pose)
        registration['estimate'].update(mattress_source=pad_source, shoulder_distance_source=distance_source)
        placement_status = registration['status']
        estimates.append('M_T<-C')
    data['patient']['ct_registration'] = registration
    # 2. t_0, m_0 ----------------------------------------------------------------
    state = room_reference_state(table)
    if state == 'missing':
        require_estimates('la referencia tablero–isocentro (marcador)')
        target, info = select_target(vertices_lps(), metadata, regions, registration)
        reference = estimate_room_reference(registration, target, events, reference=reference_event, target_info=info)
        geometry['table_calibration']['reference'] = reference
        estimates.append('t0, m0')
        for warning in reference['estimate']['warnings']:
            print('AVISO: ' + warning)
    return dict(ct_to_table_status=placement_status,
                room_reference_status=room_reference_state(geometry['table_calibration']),
                estimated_links=estimates, target_regions_rdsr=regions, intervention_pose=pose,
                registration=registration, room_reference=geometry['table_calibration'].get('reference'))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('folder', type=Path, nargs='?')
    parser.add_argument('--examples-root', type=Path)
    for name in ('config','rdsr','mesh'):
        parser.add_argument('--'+name)
    parser.add_argument('--espesor-colchoneta-mm', type=float, help='Espesor efectivo de la colchoneta bajo carga')
    parser.add_argument('--distancia-hombros-mm', type=float, help='Distancia de la cabecera del tablero al plano de hombros')
    parser.add_argument('--evento-referencia', default='kerma_weighted',
                        help="Evento para estimar t0/m0 si no hay marcador: kerma_weighted (posición de trabajo, "
                             "por defecto), first, max_kerma o su número")
    parser.add_argument('--permitir-estimaciones', action='store_true',
                        help='Admitir estimaciones de colocación aunque el perfil no lo declare')
    parser.add_argument('--mesa', type=Path, help='Perfil JSON de la mesa de la sala')
    parser.add_argument('--registro-tc', type=Path, help='Registro JSON del TC de este caso en el tablero')
    for flag in ('interactive','prepare-only','open-viewer'):
        parser.add_argument('--'+flag, action='store_true')
    args = parser.parse_args(argv)
    try:
        folder = args.folder
        if folder is None:
            if not args.examples_root or not args.interactive:
                raise ExampleError('Indique la carpeta del ejemplo.')
            number = input('Escribe el número del ejemplo (por ejemplo 01): ').strip()
            if not number.isdigit():
                raise ExampleError('El número del ejemplo debe contener solo dígitos.')
            folder = args.examples_root/('ejemplo'+number.zfill(2))
        data, audit = prepare_example(folder, config=args.config, rdsr=args.rdsr, mesh=args.mesh,
                                      interactive=args.interactive, table_calibration=args.mesa,
                                      ct_registration=args.registro_tc, pad_mm=args.espesor_colchoneta_mm)
        from .config import Equipment
        from .rdsr import read_rdsr
        from .azurion import prepare_events
        from .engine import normalize
        destination = folder.resolve()/'resultados'
        destination.mkdir(parents=True, exist_ok=True)
        # Remove only the obsolete report generated by earlier launcher versions.
        (destination/'colocacion_malla.html').unlink(missing_ok=True)
        events = read_rdsr(Path(audit['rdsr'])).events
        try:
            placement = place_patient(data, audit, events, interactive=args.interactive,
                                      pad_mm=args.espesor_colchoneta_mm, shoulder_distance_mm=args.distancia_hombros_mm,
                                      reference_event=args.evento_referencia,
                                      allow_estimates=True if args.permitir_estimaciones else None)
        except (ValueError, OSError) as exc:
            audit.update(registration_status='insufficient_information', registration_note=str(exc))
            (destination/'seleccion_entrada.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
            raise
        audit['placement'] = placement
        audit['registration_status'] = 'estimated' if placement['estimated_links'] else 'explicit'
        audit['registration_note'] = ('Cadena TC → tablero → isocentro medida/explícita.' if not placement['estimated_links']
                                      else 'Eslabones estimados: ' + ', '.join(placement['estimated_links'])
                                      + '. Resultado de investigación; la colocación no está calibrada.')
        base = Path(audit['source_config']).parent
        equipment = Equipment.from_dict(data, base_directory=base)
        evidence = equipment.registered_frames.evidence
        geometry_check = evidence.get('ct_slice_geometry')
        if geometry_check:
            print(f"Geometría TC verificada: {geometry_check['slice_count']} cortes, residual máximo "
                  f"{geometry_check['max_position_residual_mm']:.3f} mm.")
        if evidence.get('frame_of_reference_status') == 'uid_ambiguous_slice_grid_consistent':
            print('Aviso: UID DICOM distintos; rejilla geométrica coherente. Identidad del paciente no verificada.')
        pose = placement['intervention_pose']
        print(f"Postura RDSR: {pose.get('position') or 'no determinada'} ({pose.get('status', 'sin datos')}).")
        print(f"TC → tablero: {placement['ct_to_table_status']}. Tablero → isocentro: {placement['room_reference_status']}.")
        prepared = prepare_events(events, equipment)
        normalized = tuple(normalize(event, equipment) for event in prepared)
        for event in normalized:
            equipment.validate_event(event, 'skin_dose')
            equipment.registered_frames.table_origin_i(event)  # S must cover every table move
        audit['table_coordinate_conversions'] = [
            dict(event_index=index, **event.original['table_coordinate_conversion'])
            for index, event in enumerate(normalized, 1)
            if 'table_coordinate_conversion' in event.original]
        effective = destination/'configuracion_efectiva.json'
        # Persist the resolved table/registration inputs so paths in the original
        # configuration are independent of the output directory.
        effective.write_text(json.dumps(equipment.raw, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        (destination/'seleccion_entrada.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        for label, key in [('RDSR','rdsr'),('Configuración','source_config'),('Malla','mesh')]:
            print(label+': '+audit[key])
        if args.prepare_only:
            return 0
        from .cli import run
        result = run(Path(audit['rdsr']), effective, 'skin_dose', destination)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        if args.open_viewer:
            import webbrowser
            webbrowser.open((destination/'mapa_3d_interactivo.html').as_uri())
        return 0
    except (ValueError, OSError, KeyError) as exc:
        print(f'ERROR CONTROLADO: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
