# Memoria en LaTeX

La fuente científica canónica se mantiene en \`docs/Memoria_Calculador_Dosis_Piel.md\`.
\`main.tex\` permite compilar esa misma versión desde LaTeX mediante el paquete
\`markdown\`, evitando mantener dos memorias divergentes.

La memoria fue actualizada el 6 de octubre de 2026 para incorporar la implementación
de los factores tabulados de Benmakhlouf, Fransson y Andreo (KS-ASF-201101-IR):
\`k_BS\`, \`k_f\`, selección/interpolación, trazabilidad, prevención de doble conteo,
límites físicos y las 33 pruebas específicas añadidas (87 pruebas totales).

## Compilar

Desde \`docs/latex\`:

\`\`\`bash
lualatex --shell-escape main.tex
lualatex --shell-escape main.tex
\`\`\`

o con LaTeX Workshop/latexmk:

\`\`\`bash
latexmk -lualatex -shell-escape main.tex
\`\`\`

La opción \`--shell-escape\` puede ser necesaria para el paquete \`markdown\`,
dependiendo de la distribución TeX instalada.
