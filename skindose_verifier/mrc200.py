"""Documented MRC200 constraints; catalog minima are not measured beam quality."""
import math
from .domain import InputError

def validate_configuration(events, equipment):
    spec = equipment.raw.get('tube_characterization', {})
    if spec.get('model') != 'MRC 200 0407':
        raise InputError('MRC 200 0407: falta tube_characterization.model compatible.')
    al = spec.get('fixed_effective_al_mm')
    try:
        valid = math.isfinite(float(al)) and float(al) >= 2.5
    except (TypeError, ValueError):
        valid = False
    if not valid:
        raise InputError('MRC 200: Al efectivo debe ser finito y >= mínimo nominal 2.5 mm; no es solo Al inherente.')
    if spec.get('filtration_basis') not in ('nominal_assembly_scenario', 'nominal_assembly_plus_collimator_scenario', 'measured_effective_model'):
        raise InputError('MRC 200: declarar la base de filtración efectiva.')
    if not math.isclose(float(equipment.raw['surface_factors']['inherent_al_mm']), float(al), abs_tol=1e-8):
        raise InputError('MRC 200: filtración del perfil y de las tablas inconsistente.')
    if equipment.raw['surface_factors'].get('expected_anode_angle_deg') != 11:
        raise InputError('MRC200: exigir tablas de 11 grados explícitamente.')
    for event in events:
        if event.kvp is not None and (not math.isfinite(event.kvp) or event.kvp > 125):
            raise InputError(f'MRC 200 0407: evento {event.uid} supera 125 kV del conjunto; verificar identidad/dato, no recortar.')
