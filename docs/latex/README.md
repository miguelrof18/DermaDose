# Memoria DermaDose en LaTeX

La fuente científica canónica es `docs/Memoria_Calculador_Dosis_Piel.md`.
El archivo `docs/latex/main.tex` compila directamente esa versión para evitar mantener
dos memorias divergentes.

La versión del 6 de octubre de 2026 incorpora:

- factores tabulados de Benmakhlouf, Fransson y Andreo (`k_BS` y `k_f`);
- selección e interpolación dentro del dominio publicado, sin extrapolación silenciosa;
- controles contra doble contabilización de factores;
- verificación numérica independiente y ampliación a 87 pruebas;
- limitaciones físicas actualizadas;
- plan experimental ampliado para C_K, HVL, transmisión de mesa/colchoneta,
  comprobación de k_BS·k_f, geometría, perfiles, oblicuidad e incertidumbre.

## Compilación en VS Code

Usa LuaLaTeX. Con LaTeX Workshop selecciona una receta `latexmk -lualatex` o ejecuta:

```bash
latexmk -lualatex -shell-escape main.tex
```

La fuente monoespaciada se fija en Consolas para evitar el problema anterior con
DejaVu Sans Mono en Windows.

Las imágenes se buscan con rutas relativas a `docs/`. Si una copia local no contiene
los ficheros gráficos referenciados por la memoria, LaTeX informará de ello aunque el
texto científico siga siendo la fuente canónica.
