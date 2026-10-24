# Philips Azurion: soporte de lectura y perfil explícito

El soporte conserva el RDSR y diferencia «No filter» (DCM 111609, tipo 113772) de filtración ausente. Un No filter contradictorio con material/espesor se rechaza. El aluminio identificado sin cobre implica 0 mm Cu; información ausente permanece desconocida.

## Filtración y tablas

`surface_factors.aluminium_policy: fixed_plus_event` opta explícitamente por seleccionar la familia Al con el valor fijo configurado más el Al añadido declarado en el evento. El nombre heredado `inherent_al_mm` pasa, bajo esta política, a representar la hipótesis/caracterización del Al fijo efectivo. El audit registra fijo, añadido y familia seleccionada. La política predeterminada anterior `reject_additional` se conserva para otras configuraciones. No se convierte Al en Cu equivalente; no se interpola entre familias de Al y no se extrapola. No se modifica ninguna tabla publicada.

Los documentos Philips describen una filtración mínima del tubo equivalente a 2,5 mm Al, y SpectraBeam de 0,1/0,4/0,9 mm Cu + 1 mm Al. El mínimo no identifica la filtración total de la instalación ni justifica 3,1 mm Al. El mapeo espectral es aproximado y requiere contrastar calidad del haz y capas hemirreductoras, además de identificación del tubo/versión y configuración del limitador de haz.

Fuentes: IFU Azurion Release 1.2, sección 16.21.5 (pp. impresas 355–356):
https://www.documents.philips.com/assets/Instruction%20for%20Use/20201222/c122581e1172411b8c0eac9900949a77.pdf
ClarityIQ white paper, nota 41:
https://www.philips.com/c-dam/b2bhc/master/landing-pages/queen-alexandra/philips-igt-azurion-clarityiq-whitepaper.pdf

## Configuraciones

`config/equipo_philips_azurion_caracterizacion.json` es una plantilla bloqueada hasta caracterizar filtración fija, geometría y calibración. No basta con marcar los indicadores: antes hay que sustituir parámetros por datos de la instalación.

`config/equipo_philips_azurion_investigacion.json` permite hipótesis explícitas y auditadas: Al fijo 3,1 mm, geometría de demostración, referencia de mesa del primer evento, factores de calibración/mesa sintéticos y aproximación de campo rectangular por cuadrado de igual área. Sustituye distancias ausentes/no positivas SOLO porque `invalid_distance_policy=configured_research` y `allow_research_assumptions=true` están activados. Mantiene distancias válidas del evento. No es perfil físico validado de la instalación. Revisar malla, ejes y registro de paciente antes de interpretar cualquier mapa. Las rutas de malla son relativas al JSON; en ejemplo07 deben adaptarse.

En perfiles Azurion se rechaza fabricante/modelo diferente y se desactiva el fallback automático a factores históricos sintéticos. No se ignoran silenciosamente eventos fuera de dominio.

## Inspección antes de calcular

Desde la raíz del proyecto, en PowerShell:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.azurion "..\ejemplos rdsr\ejemplo07\0_anon.dcm" --fixed-al-mm 3.1
```

No genera dosis ni necesita malla. El parámetro `--fixed-al-mm` es una hipótesis de inspección; informa compatibilidad del eje de filtración, no de kVp/campo/registro completos.

En el RDSR aportado se leen 101 eventos: 30 explícitamente sin filtro, 47 con 1 mm Al + 0,1 mm Cu y 24 con 1 mm Al + 0,4 mm Cu. Hay un evento con kerma no nulo y SID/SOD/SRP=0; no puede inferirse su geometría clínica. Con Al fijo hipotético 3,1, la familia 4,1 mm Al solo cubre Cu de 0 a 0,2 mm en el dataset actual: los 24 eventos de 0,4 mm Cu siguen fuera de dominio. No se afirma que todo el procedimiento sea calculable con Benmakhlouf.

Para completarlo se requiere caracterizar la instalación y, para calidades no cubiertas, un dataset propio compatible con `provider=tabulated_water`, mediciones o cálculos espectrales independientes. No introducir factores inventados para completar la tabla. El informe DICOM del usuario no se incluye en el repositorio ni en el paquete.
