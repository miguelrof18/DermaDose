# Ejemplo XCAT

Esta carpeta contiene el maniquí `xcat_wholebody_body.stl` y su licencia. La demostración usa la misma malla incluida en el paquete y genera un procedimiento ficticio. Para ejecutarla, instala DermaDose según el [README principal](../../README.md) y, desde la raíz del repositorio, ejecuta:

```powershell
.\.venv\Scripts\python.exe -m skindose_verifier demo --model mesh --out resultados/ejemplo_xcat --stage skin_dose
```

Abre `resultados/ejemplo_xcat/mapa_3d_interactivo.html` para ver el mapa y `resultados/ejemplo_xcat/resultados.json` para consultar la dosis pico y la trazabilidad.

Para usar XCAT con tu propio RDSR del Azurion, sigue el [ejemplo de equipo](../ejemplo/README.md) y selecciona `--mesh xcat`.
