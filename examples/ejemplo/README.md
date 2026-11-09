# Ejemplo

Esta carpeta contiene un perfil de investigación para Philips Azurion. El RDSR y la superficie del paciente se aportan localmente.

1. Instala DermaDose según el [README principal](../../README.md).
2. Crea una carpeta local, por ejemplo `C:\casos\ejemplo`, y añade un único RDSR DICOM del Azurion.
3. Desde la raíz del repositorio, ejecuta con XCAT:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "C:\casos\ejemplo" --config ".\examples\ejemplo\equipo_azurion.json" --mesh xcat
```

Para usar la superficie del paciente, añade a la carpeta local el STL y su JSON homónimo exportados por DermaDose Surface y cambia `--mesh xcat` por `--mesh patient`. Puedes añadir `--prepare-only` para comprobar las entradas antes del cálculo.

Abre `resultados/mapa_3d_interactivo.html` en la carpeta local del caso. La configuración utiliza una referencia provisional de mesa y factores nominales para investigación.
