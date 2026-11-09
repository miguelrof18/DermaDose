# Ejemplo

Esta carpeta contiene un perfil de investigación para Philips Azurion. El RDSR y la superficie del paciente se aportan localmente.

1. Instala DermaDose según el [README principal](../../README.md).
2. Crea una carpeta local, por ejemplo `C:\casos\ejemplo`, y añade un único RDSR DICOM del Azurion, una copia de `equipo_azurion.json`, el STL de la superficie y su JSON homónimo de DermaDose Surface.
3. Desde la raíz del repositorio, ejecuta:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier.example_runner "C:\casos\ejemplo"
```

Puedes añadir `--prepare-only` para comprobar las entradas antes del cálculo.

Abre `resultados/mapa_3d_interactivo.html` en la carpeta local del caso. La configuración utiliza una referencia provisional de mesa y factores nominales para investigación.
