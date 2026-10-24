"""Versioned water-surface factor providers; data, interpolation and dose separate."""
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
import csv
import hashlib
import io
import json
import math

import numpy as np

from .domain import InputError
from .factor_tables import RectilinearFactorTable, FactorPair


class WaterSurfaceDataset:
    """Read the published CSV schema, also usable for future independent datasets.

    Alternate datasets must declare the same physical formalism, their own version,
    provenance, hashes and complete grids. They never overwrite published files.
    """
    def __init__(self, directory):
        self.directory = directory
        self.metadata = json.loads(directory.joinpath('metadata.json').read_text(encoding='utf-8'))
        m = self.metadata
        if (m.get('schema_version') != 1 or m.get('medium') != 'water' or
                m.get('formalism') != 'B_air_times_surface_water_air_ratio'):
            raise InputError('Dataset incompatible: se requiere B_air y cociente agua/aire superficial.')
        for key in ('dataset_id', 'version', 'source', 'quantities', 'geometry', 'units', 'tables', 'files_sha256'):
            if not m.get(key):
                raise InputError(f'Dataset: falta metadato {key}.')
        self.tables = {}
        for name in ('clinical', 'monoenergetic'):
            if name not in m['tables']:
                continue
            spec = m['tables'][name]
            expected_axes = ['kvp', 'added_cu_mm', 'field_side_cm'] if name == 'clinical' else ['energy_kev', 'field_side_cm']
            expected_groups = ['inherent_al_mm'] if name == 'clinical' else []
            if spec.get('axes') != expected_axes or spec.get('groups') != expected_groups:
                raise InputError(f'Dataset: esquema de ejes incompatible en {name}.')
            content = directory.joinpath(spec['file']).read_bytes()
            if hashlib.sha256(content).hexdigest() != m['files_sha256'].get(spec['file']):
                raise InputError(f'Dataset alterado: SHA-256 incorrecto en {spec["file"]}.')
            rows = list(csv.DictReader(io.StringIO(content.decode('utf-8'))))
            if len(rows) != spec['rows']:
                raise InputError('Dataset: número de filas incoherente.')
            groups = {}
            for row in rows:
                key = tuple(float(row[g]) for g in expected_groups)
                groups.setdefault(key, []).append(row)
            self.tables[name] = {k: RectilinearFactorTable(v, expected_axes) for k,v in groups.items()}

    def evaluate(self, kind, *, inherent_al_mm=None, **coordinates):
        key = () if kind == 'monoenergetic' else (inherent_al_mm,)
        if kind not in self.tables or key not in self.tables[kind]:
            raise InputError(f'{kind}: familia Al={inherent_al_mm} no tabulada; no se interpola entre familias Al.')
        result = self.tables[kind][key].evaluate(**coordinates)
        issues = sorted({r['quality_flag'] for r in result.audit['source_nodes'] if r.get('quality_flag')})
        spec = self.metadata['tables'][kind]
        return FactorPair(result.k_bs, result.k_f, {**result.audit,
            'dataset_id':self.metadata['dataset_id'], 'dataset_version':self.metadata['version'],
            'source':self.metadata['source'], 'table':spec['source_table'],
            'data_sha256':self.metadata['files_sha256'][spec['file']],
            'inherent_al_mm':inherent_al_mm, 'source_quality_flags':issues,
            'medium':'water', 'formalism':self.metadata['formalism'], 'clinical_validation':False})


@lru_cache(maxsize=1)
def benmakhlouf_2011():
    return WaterSurfaceDataset(files('skindose_verifier').joinpath('data/benmakhlouf_2011'))


class EventSurfaceFactors:
    """Explicitly activated paired provider; incompatible legacy factors rejected."""
    def __init__(self, spec, dataset):
        self.spec = dict(spec)
        self.dataset = dataset
        allowed = {'provider', 'dataset_directory', 'inherent_al_mm', 'target', 'field_shape_policy',
                   'out_of_domain', 'incident_kerma_basis', 'quality_basis', 'provenance',
                   'allow_water_phantom_approximation', 'aluminium_policy'}
        if set(spec)-allowed:
            raise InputError(f'Opciones de surface_factors desconocidas: {sorted(set(spec)-allowed)}.')
        if spec.get('target') != 'W':
            raise InputError('Tabla clínica: requiere target=W explícito.')
        try:
            al = float(spec['inherent_al_mm'])
        except (KeyError, TypeError, ValueError) as exc:
            raise InputError('Declarar inherent_al_mm caracterizado; no se deduce del RDSR.') from exc
        if (al,) not in dataset.tables.get('clinical', {}):
            raise InputError('Familia de filtración Al no tabulada; no se extrapola.')
        self.al = al
        self.aluminium_policy = spec.get('aluminium_policy', 'reject_additional')
        if self.aluminium_policy not in ('reject_additional', 'fixed_plus_event'):
            raise InputError('aluminium_policy inválida; usar reject_additional o fixed_plus_event.')
        if spec.get('out_of_domain', 'error') != 'error':
            raise InputError('La única política implementada fuera del dominio es error.')
        if spec.get('field_shape_policy', 'square_only') not in ('square_only', 'equal_area'):
            raise InputError('field_shape_policy debe ser square_only o equal_area.')
        if spec.get('incident_kerma_basis') != 'free_in_air_without_backscatter':
            raise InputError('Declarar kerma incidente free_in_air_without_backscatter para evitar doble conteo.')
        if spec.get('quality_basis') != 'event_kvp_filtration_proxy_at_surface' or not spec.get('provenance'):
            raise InputError('Declarar quality_basis y provenance: la calidad del evento es un proxy del espectro incidente superficial.')
        if spec.get('allow_water_phantom_approximation') is not True:
            raise InputError('Declarar allow_water_phantom_approximation=true para trasladar el fantoma al mapa.')

    def evaluate_event(self, event, field_area_cm2):
        added_al = [part for part in event.original.get('filter_stack', [])
                    if (part['material_scheme'], part['material_code']) == ('SRT', 'C-120F9')
                    and part['selected_mm'] > 0]
        synthetic_spectrum = event.quality.get('synthetic_profile_spectrum') == 'additional_al_not_modeled'
        selected_al = self.al
        if self.aluminium_policy == 'fixed_plus_event':
            if event.quality.get('filtration_mm_cu') == 'missing_filter_information' or event.filtration_mm_cu is None:
                raise InputError(f'Evento {event.uid}: filtración desconocida; no se supone cero.')
            selected_al += sum(float(part['selected_mm']) for part in added_al)
            families = [k[0] for k in self.dataset.tables.get('clinical', {})]
            selected_al = next((a for a in families if math.isclose(a, selected_al, abs_tol=1e-8, rel_tol=0)), selected_al)
            if (selected_al,) not in self.dataset.tables.get('clinical', {}):
                raise InputError(f'Evento {event.uid}: Al fijo + añadido = {selected_al:g} mm no tabulado; no se extrapola.')
        if added_al and not synthetic_spectrum and self.aluminium_policy == 'reject_additional':
            raise InputError(f'Evento {event.uid}: el RDSR declara filtración adicional de Al; '
                             'no se puede equiparar a la familia de Al inherente de Benmakhlouf.')
        square = math.isclose(event.field_width_mm, event.field_height_mm, rel_tol=1e-9, abs_tol=1e-9)
        policy = self.spec.get('field_shape_policy', 'square_only')
        if not square and policy == 'square_only':
            raise InputError(f'Evento {event.uid}: campo no cuadrado; declarar equal_area si se acepta esa aproximación.')
        try:
            result = self.dataset.evaluate('clinical', inherent_al_mm=selected_al, kvp=event.kvp,
                added_cu_mm=event.filtration_mm_cu, field_side_cm=np.sqrt(field_area_cm2))
        except InputError as exc:
            raise InputError(f'Evento {event.uid}: {exc}') from exc
        return FactorPair(result.k_bs, result.k_f, {**result.audit,
            'event_uid':event.uid, 'configuration':self.spec,
            'aluminium_mapping': {'policy': self.aluminium_policy, 'fixed_al_mm': self.al,
                                  'event_added_al_mm': sum(float(p['selected_mm']) for p in added_al),
                                  'selected_table_al_mm': selected_al},
            'field_mapping':'square_side_at_local_axial_depth' if square else 'equal_area_square_at_local_axial_depth',
            'approximations':([
                'RDSR additional aluminium not modeled by the configured synthetic equipment profile; spectral quality not reconstructed.'
            ] if added_al and synthetic_spectrum else []) + [
                'Water phantom central-axis factors applied locally to patient surface; no anatomical transport.',
                'Normal-incidence tables reused without an obliquity correction or field-edge scatter model.',
                '100 cm FSD reference reused at modeled distance; inverse-square correction remains separate.',
                'Event kVp/Al/Cu used as spectral proxy; scalar table transmission does not model spectral hardening.'
            ] + ([] if square else ['Rectangular field replaced by equal-area square; not validated by source.'])})


def make_event_provider(spec, *, base_directory=Path('.')):
    if not isinstance(spec, dict):
        raise InputError('surface_factors debe ser un objeto de configuración.')
    if spec.get('provider') == 'benmakhlouf_2011':
        if 'dataset_directory' in spec:
            raise InputError('Los datos publicados no se sustituyen; utilizar provider=tabulated_water para datos propios.')
        dataset = benmakhlouf_2011()
    elif spec.get('provider') == 'tabulated_water':
        if not spec.get('dataset_directory'):
            raise InputError('Dataset propio: falta dataset_directory.')
        directory = Path(spec['dataset_directory'])
        if not directory.is_absolute():
            directory = base_directory/directory
        spec['dataset_directory'] = str(directory.resolve())
        dataset = WaterSurfaceDataset(directory)
    else:
        raise InputError('Proveedor de factores superficiales desconocido.')
    return EventSurfaceFactors(spec, dataset)
