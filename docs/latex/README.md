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

La fuente monoespaciada usa Consolas cuando está instalada y Latin Modern Mono
como alternativa portable. La revisión del 4 de noviembre de 2026 se incluye en la fuente canónica.

Las imágenes se buscan con rutas relativas a `docs/`. Si una copia local no contiene
los ficheros gráficos referenciados por la memoria, LaTeX informará de ello aunque el
texto científico siga siendo la fuente canónica.

## Comprobación del 4 de noviembre de 2026

La compilación se intentó con LuaLaTeX y se detuvo por ausencia de `spanish.ldf`.
En TeX Live se requiere el soporte de español (por ejemplo, `texlive-lang-spanish`).
Además faltan en el repositorio las figuras históricas de `docs/figuras/` y
`docs/evidencias/` referenciadas por la memoria. Deben recuperarse de sus resultados
originales; no se ha generado ni verificado un PDF completo en esta revisión.
