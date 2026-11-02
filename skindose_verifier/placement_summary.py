"""Readable configured transforms; anatomical registration is not inferred."""
from html import escape
import math

DIRECTIONS = (('izquierda', 'derecha'), ('cabeza', 'pies'), ('anterior', 'posterior'))


def direction(vector):
    parts = []
    for value, (positive, negative) in zip(vector, DIRECTIONS):
        if abs(value) > 1e-7:
            label = positive if value > 0 else negative
            parts.append(label if math.isclose(abs(value), 1, abs_tol=1e-7)
                         else f'{abs(value):.3f} hacia {label}')
    return ' + '.join(parts) or 'sin componente'


def placement_rows(patient):
    matrix = patient['mesh_axes_matrix']
    rows = [('Unidades del STL', patient['mesh_unit'])]
    for column, axis in enumerate('XYZ'):
        rows.append((f'Eje +{axis} del STL', direction([matrix[i][column] for i in range(3)])))
    for key, title in (('mesh_origin_mm', 'Traslación de la malla (antes de inclinar)'),
                       ('origin_on_table_mm', 'Desplazamiento adicional respecto a la mesa')):
        parts = []
        for value, (positive, negative) in zip(patient[key], DIRECTIONS):
            if abs(value) > 1e-7:
                parts.append(f'{abs(value)/10:g} cm hacia {positive if value > 0 else negative}')
        rows.append((title, '; '.join(parts) or 'sin desplazamiento'))
    rows.append(('Inclinación adicional',
                 f"X: {patient.get('tilt_x_deg', 0):g} grados; Y: {patient.get('tilt_y_deg', 0):g} grados (primero X, después Y)"))
    return rows


NOTE = ('Colocación provisional: se aplica la configuración del equipo, completada con '
        'valores por defecto si faltan. No confirma la posición real del paciente. '
        'Los nombres anatómicos indican los ejes que asume el programa; '
        'no se ha identificado automáticamente la anatomía del STL. '
        'Los movimientos de mesa del RDSR se aplican después, por evento.')


def console_summary(patient):
    lines = ['\nCOLOCACIÓN QUE SE APLICARÁ A LA MALLA']
    lines.extend(f'  {label}: {value}' for label, value in placement_rows(patient))
    lines.extend(['  Convención: +X izquierda; +Y cabeza; +Z anterior.', '  ' + NOTE,
                  '  Se continúa sin solicitar datos de orientación.'])
    return '\n'.join(lines)


def html_summary(patient, audit):
    rows = ''.join(f'<tr><th>{escape(label)}</th><td>{escape(value)}</td></tr>'
                   for label, value in placement_rows(patient))
    return '''<!doctype html><html lang="es"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Colocación de la malla · DermaDose</title>
<style>body{font:17px/1.6 system-ui;background:#eef3f7;color:#173047;margin:0;padding:4vw}
main{max-width:900px;margin:auto;background:white;border-radius:20px;padding:32px}
h1{line-height:1.2}small{color:#536979}table{width:100%;border-collapse:collapse}
th,td{text-align:left;padding:16px;border-bottom:1px solid #dce5ed}th{width:48%}
.axes{display:flex;gap:12px;flex-wrap:wrap;margin:24px 0}.axes span{padding:12px 18px;
border-radius:10px;background:#e5f0ff}.axes span:nth-child(2){background:#e4f5eb}
.axes span:nth-child(3){background:#f4eaff}.note{background:#fff4da;padding:20px;border-radius:12px}
details{margin-top:24px}pre{overflow:auto;font-size:13px}.file{overflow-wrap:anywhere}
</style><main><small>DERMADOSE · CONFIGURACIÓN INICIAL</small>
<h1>Así se colocará la malla</h1><p class="file">''' + escape(audit['mesh']) + '''</p>
<p>Direcciones del sistema de coordenadas que asume el programa:</p>
<div class="axes"><span>+X · Izquierda<br>−X · Derecha</span>
<span>+Y · Cabeza<br>−Y · Pies</span><span>+Z · Anterior<br>−Z · Posterior</span></div>
<table>''' + rows + '''</table><p class="note">''' + escape(NOTE) + '''</p>
<p>Información de la colocación configurada. El cálculo continúa automáticamente.</p>
<details><summary>Ver matriz técnica de orientación del STL</summary><pre>''' + escape(
        '\n'.join(str(row) for row in patient['mesh_axes_matrix'])) + '''</pre></details></main></html>'''
