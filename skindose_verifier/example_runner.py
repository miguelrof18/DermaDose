"""Ejemplos con RDSR, equipo JSON y STL locales; entradas inmutables."""
from __future__ import annotations
import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys


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
                    data_directory=None):
    folder = Path(folder).resolve()
    if not folder.is_dir():
        raise ExampleError(f'No existe la carpeta: {folder}')
    resources = Path(data_directory) if data_directory else Path(__file__).resolve().parent/'data'
    entries = [p for p in folder.iterdir() if p.is_file()]
    def local(value):
        return None if value is None else folder/Path(value)
    if mesh is None and interactive:
        print('¿Qué superficie quieres utilizar?')
        print('  1: Maniquí XCAT incluido en el programa')
        print('  2: TC del paciente (malla STL de la carpeta del ejemplo)')
        answer = input_fn('Opción [1/2]: ').strip()
        if answer not in ('1', '2'):
            raise ExampleError('Seleccione 1 para XCAT o 2 para la malla del paciente.')
        mesh = 'xcat' if answer == '1' else 'patient'
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
    config_path = choose(configs, 'configuración de equipo', local(config), interactive, input_fn)
    data = deepcopy(read_json(config_path))
    dicom = choose([p for p in entries if p.suffix.lower()=='.dcm'], 'RDSR (.dcm)', local(rdsr), interactive, input_fn)
    candidates = [p for p in entries if p.suffix.lower()=='.stl']
    selected = resources/'xcat_wholebody_body.stl' if mesh=='xcat' else choose(candidates, 'malla STL del paciente', local(mesh) if mesh not in (None, 'patient') else None, interactive, input_fn)
    selected = selected.resolve()
    if not selected.is_file():
        raise ExampleError(f'No existe la malla: {selected}')
    # Provisional placement; no patient file or registration confirmation.
    defaults = dict(mesh_unit='mm', mesh_axes_matrix=[[1,0,0],[0,1,0],[0,0,1]],
                    mesh_origin_mm=[0,0,0], origin_on_table_mm=[0,0,0],
                    resolution_mm=5, mesh_max_edge_mm=None, thickness_mm=220)
    data['patient'] = {**defaults, **data.get('patient', {})}
    placement_source = 'equipment_patient_assumed' if 'patient' in read_json(config_path) else 'provisional_defaults_assumed'
    data['patient'].update(type='mesh', mesh_path=str(selected), phantom_name='XCAT50' if mesh=='xcat' else selected.stem)
    data.pop('phantom', None)  # legacy blocks must not replace the selected STL
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
            known = {'derived_mrc200_nominal','derived_azurion_spekpy','derived_azurion_qc_conditional'}
            if not absolute and len(parts)>=3 and parts[-3:-1]==['skindose_verifier','data'] and parts[-1] in known:
                directory = resources/parts[-1]
            if not (directory/'metadata.json').is_file():
                raise ExampleError(f'Falta metadata.json de las tablas espectrales: {directory}')
        sf['dataset_directory'] = str(directory.resolve())
    audit = dict(source_config=str(config_path), source_config_sha256=sha256(config_path.read_bytes()).hexdigest(),
                 rdsr=str(dicom), mesh=str(selected), mesh_sha256=digest,
                 placement_source=placement_source, registration_status='not_established',
                 mesh_choice='xcat' if mesh=='xcat' else 'patient_stl',
                 dataset_directory=sf.get('dataset_directory'),
                 registration_note='Seleccionar un STL no establece el registro paciente-equipo.')
    return data, audit


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path, nargs='?')
    parser.add_argument('--examples-root', type=Path)
    for name in ('config','rdsr','mesh'):
        parser.add_argument('--'+name)
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
                                      interactive=args.interactive)
        from .config import Equipment
        Equipment.from_dict(data, base_directory=Path(audit['source_config']).parent)
        destination = folder.resolve()/'resultados'
        destination.mkdir(parents=True, exist_ok=True)
        effective = destination/'configuracion_efectiva.json'
        effective.write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        (destination/'seleccion_entrada.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
        for label, key in [('RDSR','rdsr'),('Configuración','source_config'),('Malla','mesh')]:
            print(label+': '+audit[key])
        from .placement_summary import console_summary, html_summary
        print(console_summary(data['patient']))
        summary = destination/'colocacion_malla.html'
        summary.write_text(html_summary(data['patient'], audit), encoding='utf-8')
        if args.open_viewer:
            import webbrowser
            webbrowser.open(summary.as_uri())
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
