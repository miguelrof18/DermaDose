# Revisión científica y decisiones de implementación

**Calculador de dosis en piel, fase de I+D v0.5.0 · 28-09-2026**

Esta revisión orienta mejoras de un calculador independiente de **intervencionismo de kV**. Se distingue evidencia de publicaciones, comprobación de nuestro código y propuestas propias. No es una revisión sistemática PRISMA ni demuestra superioridad clínica.

## Búsqueda y alcance de la lectura

Se revisaron los once PDF aportados, con lectura dirigida de métodos, discusión, limitaciones y conclusiones; en las tesis se revisaron los capítulos pertinentes, no se afirma una lectura exhaustiva de todas las páginas. El documento del grupo SEFM de una página describe objetivos, no valida un algoritmo. Se contrastaron fuentes primarias en PubMed/PMC, Medical Physics, JACMP, Physica Medica, Revista de Física Médica, AAPM/EFOMP, IAEA, ICRP/ICRU, BIPM/JCGM, DICOM y documentación de software. Se buscaron también ESTRO y dosimetría superficial en radioterapia. No se obtuvo una guía ESTRO específica que justifique trasladar un algoritmo MV a este caso; no se ha usado un abstract de congreso como calibración. No hubo acceso directo verificable a Google Scholar; la búsqueda se hizo mediante el buscador web y los portales anteriores.

Familias de consultas utilizadas: `skin dose mapping fluoroscopically guided interventions uncertainty mesh validation`, `skin dose future scatter fluoroscopy`, `AAPM TG357 EFOMP`, `patient surface mesh dosimetry`, `surface dose radiotherapy build-up Monte Carlo oblique incidence`, `PySkinDose`, `MCGPU-IR`, `Trimesh ray triangle`, `LatinHypercube`, `JCGM 101`. Se verificaron enlaces/DOI disponibles. Las búsquedas específicas de MCGPU-IR no aportaron un repositorio adicional confirmado; se utilizó el contenido de las tesis aportadas. El estado de lectura de cada referencia figura al final.

## Qué respalda la literatura

Los trabajos de mapeo analítico [R2–R5] apoyan reconstruir geometría por evento, normalizar Ka,r, corregir distancia y mesa y aplicar factores dosimétricos caracterizados. El detalle visual de un fantoma no garantiza su correspondencia con el paciente. Johnson et al. [R3] y Omar et al. [R8] muestran la importancia de la anatomía y del registro paciente–haz. La comparación experimental reciente de Vergara et al. [R13] vuelve a situar la geometría y los factores medidos en el centro del problema; sus resultados en dos procedimientos no se transfieren al XCAT de este proyecto.

El informe AAPM–EFOMP [R1] se toma como marco técnico, no como certificado del programa. La revisión de Malchair et al. [R5] encuentra heterogeneidad entre productos y diferencias entre validación en fantoma y en pacientes. Por tanto, “independiente del fabricante” describe un objetivo: nuestro parser sigue necesitando caracterización y pruebas por equipo.

Los métodos MC [R7,R12] ofrecen vías para modelar transporte, a costa de información y control estadístico adicionales. Han et al. [R12] proponen MC de baja estadística más reducción iterativa de ruido; además advierten del sesgo del máximo de celdas ruidosas. Es evidencia para diseñar una comparación externa, no para suavizar arbitrariamente nuestra dosis determinista ni atribuirle precisión MC.

## Modelo matemático implementado

Para un evento e, un centroide local p y fuente s, se calcula r = p − s, distancia d = |r| y profundidad axial t = r·b, con b unitario. La colimación rectangular exige:

$$0<t\le SID,\qquad |r\cdot u|\le\frac{w t}{2d_{ref}},\qquad |r\cdot v|\le\frac{h t}{2d_{ref}}.$$

Se exige normal exterior dirigida hacia la fuente y, en el mesh, que la primera cara interceptada sea la cara evaluada. Los vectores están en un marco común. Para el ángulo cutáneo:

$$\alpha=\arccos\left(\frac{n\cdot(s-p)}{|s-p|}\right).$$

La dosis aproximada de entrada, bajo las hipótesis del método analítico [R2,R4], es:

$$D_e(p)=M_e(p)\,K_{a,r,e}\left(\frac{d_{IRP}}{d}\right)^2 CF_e\,P_e(p)\,T_e(p)\,BSF_e(p)\,f_{skin,e}(p).$$

M es la máscara de entrada; P es uno por defecto o un perfil **medido** relativo al eje central; T es uno si el rayo no cruza la mesa. La etapa `free_air` omite los factores; `calibrated_air` aplica CF, P y T; `skin_dose` aplica también BSF y f_skin. Una tabla no puede extrapolarse fuera de su calibración. P no se deduce de un valor universal del efecto talón. Si se suministra, el usuario debe normalizarlo al eje central de forma compatible con Ka,r/CF; debe comprobar también su dependencia con kVp/Cu/campo.

La sección del haz normal a b a profundidad t tiene área `w*h*(t/d_ref)^2`; se convierte de mm² a cm² dividiendo por 100. Esta es la corrección frente a v0.4, que usaba d en ese argumento. La huella oblicua sobre piel y el área usada en una tabla BSF no se identifican sin justificar la convención de calibración. Un escalar de área tampoco representa por completo forma de campo, curvatura o retrodispersión lateral.

La suma se hace en coordenadas de paciente: D(p)=ΣD_e(p), y después se busca el máximo. No se suman máximos de campos separados. La composición angular del ejemplo conserva exactamente `Ry(a) @ Rx(b)` para `primary_then_secondary`; ese nombre no sustituye la verificación física del signo y orden en cada equipo.

En un triángulo se usa su centroide como punto representativo. La integral superficial aproximada es Σ A_i D_i. El refinamiento subdivide las aristas largas compartidas por ambas caras, conserva orientación, área y volumen, y reduce el tamaño de las celdas sin cambiar la anatomía. En bordes de colimación discontinuos no se espera convergencia estrictamente monótona ni se aplica extrapolación de Richardson sin demostrar orden asintótico. La integral dosis–área de piel no es el KAP del haz ni energía absorbida.

## Mejoras candidatas y priorización

Costes cualitativos respecto al motor previo; impacto esperado no equivale a exactitud demostrada. “Sí” significa implementado como capacidad, no calibrado para uso clínico.

| Mejora | Problema y fundamento | Evidencia | Dificultad | Coste | Impacto esperado | Validación posible | Decisión |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Refinamiento conformante global/regional | Integración insuficiente en triángulos grandes | [R3,R6,R17]; diseño propio de subdivisión | Media | Aumenta caras y rayos | Menor discretización, medible | Área/volumen + convergencia | Sí |
| Normales de cara robustas y primera entrada | Caras ocultas/posteriores, orientación | [R4,R17] | Media | Ray tracing acotado | Evitar irradiación ficticia | Esfera, cilindro, brazos | Sí, reforzado |
| Área perpendicular del campo | Argumento geométrico de tablas sesgado fuera del eje | Derivación geométrica; [R2] para formalismo | Baja | Despreciable | Corrección física localizada | Solución analítica | Sí |
| BSF angular calibrado | Orientación y calidad no reducibles a un coseno | [R2,R11] | Baja en API; alta en medidas | Interpolación vectorizada | Representar datos medidos | Fantoma angular | Sí, capacidad; sin medidas reales |
| Perfil de haz/heel medido | Haz uniforme | [R2,R11] | Media | Una interpolación | Mejor descripción transversal | Película/cámara en perfiles | Sí, opcional |
| Caché exacta de geometría | Recalcular eventos idénticos | [R9] propone caché; implementación propia | Baja | Memoria para dos proyecciones | Acelerar repeticiones | Igualdad de mapas + benchmark | Sí |
| LHS/MC de parámetros | Incertidumbre geométrica y picos no lineales | [R9,R14,R18] | Media | N recálculos | Distribución condicionada de PSD y mapas | Caso analítico + semillas/N | Sí |
| DAH y clustering de hotspots | Pico aislado no describe extensión | Mapeo [R3,R4]; métricas exploratorias propias | Media | Grafo + suma de áreas | Inspección más completa | Casos con áreas conocidas | Sí; no predictor clínico |
| Supersampling baricéntrico sin cambiar malla | Campo parcial por cara | Cuadratura numérica; propuesta propia | Media | Varias muestras/cara | Área menos sensible al borde | Integral plana exacta | Priorizado refinamiento real; pendiente |
| Refinamiento automático guiado por error | Gastar resolución donde más importa | [R6] motiva resolución; criterio propio | Alta | Iterativo | Menos caras para tolerancia dada | Comparar con refinamiento global | Solo ROI geométrica explícita por ahora |
| Anatomía individual/registro óptico o CT | Fantoma de referencia no coincide | [R3,R8,R13] | Alta | Segmentación/registro | Posible mejora dominante | Superficie medida + fantoma | Pendiente de datos |
| Reparación de autointersecciones/suavizado | Anatomía inválida o triangulación mala | [R17,R21] | Alta | Geometría adicional | Beneficio solo con control de forma | Distancias y topología antes/después | No automático |
| Dispersión fuera del campo y dosis de salida | Área subestimada fuera de huella primaria | [R7,R12] | Alta | Kernels calibrados o MC | Mejor cobertura física | Película y MC | Pendiente |
| GPU/Embree/BVH | Muchas intersecciones | [R6,R7,R17] | Media–alta | Dependencia/hardware | Aceleración a gran escala | Mismos rayos, tolerancias y memoria | No añadido sin necesidad demostrada |
| Profundidad, espesor y heterogeneidad | Superficie sin transporte | [R7,R12,R21,R25] | Alta | Modelo volumétrico | Preguntas físicas nuevas | MC y detectores adecuados | Fuera del núcleo actual |
| Enhanced RDSR y colimador poligonal | Coordenadas/filtros/campo insuficientes | [R15] | Alta | Parser + geometría | Compatibilidad mejor fundamentada | Corpus real y DICOM de referencia | Pendiente |

## Evaluación de las ampliaciones físicas solicitadas

**Adecuadas al núcleo actual:** distancia fuente–superficie, divergencia, orientación local, tamaño de campo, gantry/mesa, múltiples campos y acumulación. Ya están implementadas o reforzadas. Energía del haz entra mediante kVp/Cu y tablas documentadas; kVp no se interpreta como energía monoenergética. `patient_thickness_mm` es un argumento opcional de tablas, no una medición de trayecto ni un cálculo de densidad.

**Adecuadas con caracterización adicional:** perfiles, retrodispersión angular/espectral, dependencia de mesa y campos asimétricos. La interpolación multilineal evita introducir oscilaciones de orden alto sin evidencia; extrapolar sigue bloqueado. Normal por vértice puede servir a la iluminación/interpolación gráfica; la intersección física utiliza la normal de la cara real.

**Requieren otro modelo físico:** espesor/profundidad de piel, densidad, heterogeneidades, build-up electrónico, dosis de salida y scatter fuera del campo. Los trabajos MV como [R25] muestran por qué profundidad y detector importan, pero sus correcciones no se transfieren a kV. No se incorpora MLC, jaws de linac ni VMAT: el código no es un TPS ni lee RTPLAN. Un barrido C-arm uniforme es una hipótesis explícita, no un tratamiento VMAT.

## Ideas propuestas en la literatura para futuros desarrollos

Las propuestas de esta tabla se localizan en el texto de los autores. No se afirma que sigan inéditas o sin implementar en toda la comunidad en 2026.

| Artículo/tesis, año, autores | Propuesta localizada | Problema | Encaje e implementación en nuestro proyecto |
| --- | --- | --- | --- |
| Skin dose mapping…, 2011, Johnson et al. [R3], discusión pp. 5498–5499 | Transferir una dosis calculada en un contorno a un fantoma híbrido para visualizarla; automatizar/recibir datos casi en tiempo real | Exactitud anatómica frente a presentación y disponibilidad de RDSR | Separación cálculo/vista conservada. El registro entre anatomías distintas queda pendiente; no proyectar colores sin conservar correspondencia. |
| Fast skin dose…, 2018, Takata et al. [R7], discusión p. 238 | Mantener fantoma y materiales en memoria GPU entre simulaciones | Transferencias repetidas limitan aceleración | Se aplica el principio de reutilización en CPU. No se atribuye a nuestro motor el rendimiento GPU de ese trabajo. |
| Vendor-independent…, 2021, Krajinović, Kržanović y Ciraj-Bjelac [R4], conclusión p. 156 | Validación clínica y ampliación a procedimientos/montajes; discusión sobre armonización RDSR y posición del paciente | Generalización y registro | Parser conservador, protocolo externo y módulo posicional. Falta corpus clínico multivendedor. |
| PyMCGPU-IR…, 2023, David Fernández Bosman [R9], cap. 12, pp. 114–116 | Visor de piel, escala geométrica, posición sobre mesa, cálculo en tiempo real y caché de eventos | Presentación, anatomía, registro y coste | Visor, propagación posicional y caché exacta implementados. Escalado anatómico y flujo en tiempo real pendientes. |
| PyMCGPU-IR…, 2024, Víctor García Balcaza [R10], p. 162 | Ampliar validación variable por variable, normalizar entradas, mejorar tracking/oclusiones y seguimiento de pantalla suspendida | Fiabilidad de simulación y dosimetría ocupacional | Pruebas por componente y entradas explícitas. Tracking ocupacional/pantallas queda fuera del alcance. |
| Fast and accurate…, 2025, Han, Streitmatter, Kitahara y Lee [R12], conclusión | Integración del método en NCIRF y biblioteca anatómica ampliada, incluidas gestantes | Accesibilidad y personalización | Referencia externa y futura biblioteca; no se añadió denoising a una dosis determinista ni se simuló una gestante. |

## Propuestas propias de investigación, derivadas de lo anterior

**Mapa probabilístico condicionado:** implementado en una primera versión: media, desviación y frecuencia de superar un umbral por cara. Los parámetros son independientes entre sí y su sesgo se comparte entre todos los eventos de cada realización. Son hipótesis del experimento, no incertidumbres medidas. No se descartan silenciosamente muestras fuera de rango. `source_x/y/z_mm` traslada conjuntamente el marco fuente/isocentro/detector; no representa un foco moviéndose dentro de un cabezal fijo. SSD resulta de la geometría; no se añade una segunda incertidumbre SSD independiente que duplique posición.

**Adaptación guiada por error:** comparar centroides/subcentroides y priorizar caras cercanas a bordes y hotspots, preservando aristas conformantes. Una banda entre refinamientos no es por sí sola un intervalo probabilístico. Requiere un criterio de parada para PSD y área y pruebas con campos estrechos.

**Híbrido analítico–MC:** biblioteca de kernels para dispersión y salida, dependiente de calidad, tamaño de campo y geometría. Primero debe reproducir casos externos de referencia. No se pueden ajustar kernels a los mismos casos usados para validar y llamarlo verificación independiente.

**Surrogates, GNN y modelos físicos + aprendizaje:** podrían aproximar una corrección MC costosa; por ahora no hay conjunto de entrenamiento, cobertura de parámetros ni referencia independiente. No se implementan. Evaluar particiones por geometría/equipo, detección fuera de dominio y error en picos antes de optimizar rapidez. No se ha verificado una publicación concreta que justifique una GNN para este calculador: esta es una idea propia, no un método establecido.

**SDF y renderizado diferenciable:** posibles vías para registro anatómico y gradientes respecto a posición. Requieren control de errores de distancia y visibilidad en concavidades. No sustituyen automáticamente primera intersección ni transporte; quedan como propuestas exploratorias.

## Bibliografía verificable y nivel de acceso

- **[R1]** Andersson J et al. (2021). *Estimation of patient skin dose in fluoroscopy: summary of a joint report by AAPM TG357 and EFOMP*. Medical Physics. DOI [10.1002/mp.14910](https://doi.org/10.1002/mp.14910); [AAPM Report 357](https://www.aapm.org/pubs/reports/detail.asp?docid=216), [PMID 33930183](https://pubmed.ncbi.nlm.nih.gov/33930183/). Ficha/resumen y extractos técnicos verificados en web; no se afirma lectura completa de todos los apéndices.
- **[R2]** Gañán Mora A, Sánchez Casanueva RM, Fernández Soto JM (2024). *Metodología para la validación de programas de cálculo de dosis en piel*. Rev Fis Med 25(2):25–34. DOI [10.37004/sefm/2024.25.2.002](https://doi.org/10.37004/sefm/2024.25.2.002). PDF aportado, método/discusión revisados. Su suplemento de RDSR no se ha obtenido.
- **[R3]** Johnson PB, Borrego D, Balter S, Johnson K, Siragusa D, Bolch WE (2011). *Skin dose mapping for fluoroscopically guided interventions*. Medical Physics 38:5490–5499. DOI [10.1118/1.3633935](https://doi.org/10.1118/1.3633935). PDF aportado.
- **[R4]** Krajinović M, Kržanović N, Ciraj-Bjelac O (2021). *Vendor-independent skin dose mapping application for interventional radiology and cardiology*. JACMP 22(2):145–157. DOI [10.1002/acm2.13167](https://doi.org/10.1002/acm2.13167). PDF aportado.
- **[R5]** Malchair F, Dabin J, Deleu M, Sans Merce M, Ciraj Bjelac O, Gallagher A, Maccia C (2020). *Review of skin dose calculation software in interventional cardiology*. Physica Medica 80:75–83. DOI [10.1016/j.ejmp.2020.09.023](https://doi.org/10.1016/j.ejmp.2020.09.023). PDF aportado; revisión, no validación propia.
- **[R6]** Rana V, Rudin S, Bednarek DR (2012). *Use of a graphics processing unit (GPU) to facilitate real-time 3D graphic presentation of the patient skin-dose distribution*. Proc SPIE 8313:831343. DOI [10.1117/12.911344](https://doi.org/10.1117/12.911344). PDF aportado.
- **[R7]** Takata T et al. (2018). *Fast skin dose estimation system for interventional radiology*. Journal of Radiation Research 59(2):233–239. DOI [10.1093/jrr/rrx062](https://doi.org/10.1093/jrr/rrx062). PDF aportado.
- **[R8]** Omar A, Bujila R, Fransson A, Andreo P, Poludniowski G (2016). *A framework for organ dose estimation in x-ray angiography and interventional radiology based on dose-related data in DICOM structured reports*. Physics in Medicine & Biology 61:3063. DOI [10.1088/0031-9155/61/8/3063](https://doi.org/10.1088/0031-9155/61/8/3063). Manuscrito aceptado aportado; dosimetría de órganos, utilizado para registro/arquitectura, no equivalencia con PSD.
- **[R9]** Fernández Bosman D (2023). *PyMCGPU-IR: a new tool for patient dose monitoring in interventional radiology procedures*. Tesis UPC, PDF aportado `TDFB1de1.pdf`, cap. 12. Título/año/capítulo comprobados en el propio documento; no se inventa DOI.
- **[R10]** García Balcaza V (2024). *PyMCGPU-IR, a computational Monte Carlo application for occupational dosimetry in interventional radiology procedures*. Tesis UPC, PDF aportado `TVGB1de1.pdf`, discusión/conclusiones, p. 162. Referencia de tesis, no artículo adicional inventado.
- **[R11]** Omar A (2020). *Radiation dose and X-ray beam modelling in diagnostic and interventional radiology using Monte Carlo methods*. Karolinska Institutet. ISBN 978-91-8016-009-4. Tesis aportada; cap. 6 y modelos espectrales/off-axis.
- **[R12]** Han H, Streitmatter SW, Kitahara CM, Lee C (2025). *Fast and accurate peak skin dose estimation method for interventional fluoroscopy patients*. Medical Physics 52(4):2551–2559. DOI [10.1002/mp.17667](https://doi.org/10.1002/mp.17667); [PMC11972035](https://pmc.ncbi.nlm.nih.gov/articles/PMC11972035/). Texto obtenido en PMC; un acceso posterior devolvió comprobación de navegador. Se utilizan método y propuestas futuras visibles, sin extrapolar sus errores a nuestro código.
- **[R13]** Vergara D, Makkia RS, Li Z, Hoerner M, Fadhel MN (2025). *Experimental validation of a comprehensive fluoroscopy peak skin dose model using four different computational phantoms*. Medical Physics 52(6):3635–3648. DOI [10.1002/mp.17737](https://doi.org/10.1002/mp.17737); [PMID 40091498](https://pubmed.ncbi.nlm.nih.gov/40091498/). Resumen primario consultado; no se atribuyen “future work” no leídos.
- **[R14]** JCGM 101:2008. *Propagation of distributions using a Monte Carlo method*. [BIPM, DOI 10.59161/JCGM101-2008](https://www.bipm.org/en/doi/10.59161/jcgm101-2008). Marco de propagación; nuestra implementación LHS no se declara certificada por JCGM.
- **[R15]** DICOM PS3.17, edición web consultada 2026d, [UUUU.4.5.4: Geometry and Radiation Output Encoding in Enhanced RDSR](https://dicom.nema.org/medical/dicom/current/output/chtml/part17/sect_uuuu.4.5.4.html). Describe marcos, transformaciones y colimador. Enhanced RDSR no implementado aquí.
- **[R16]** Auer B, repositorio [Mesh-based-Human-Phantom-for-Simulation](https://github.com/BenAuer2021/Mesh-based-Human-Phantom-for-Simulation). Activo XCAT, licencia/commit/hash en `ASSET_ORIGEN.md`.
- **[R17]** Trimesh, documentación primaria de [ray tracing](https://trimesh.org/trimesh.ray.ray_triangle.html) y [remeshing](https://trimesh.org/trimesh.remesh.html). Se consultaron también funciones instaladas. El refinamiento conformante local es código propio con invariantes verificados; no se asume que todas las versiones tengan idéntico comportamiento.
- **[R18]** SciPy, [LatinHypercube](https://docs.scipy.org/doc/scipy-1.15.1/reference/generated/scipy.stats.qmc.LatinHypercube.html). Se utiliza `seed` compatible con el rango declarado; versiones ejecutadas en informes.
- **[R19]** IAEA (2007), [TRS 457, Dosimetry in Diagnostic Radiology: An International Code of Practice](https://www.iaea.org/publications/7638/dosimetry-in-diagnostic-radiology-an-international-code-of-practice). Ficha oficial consultada; contexto metrológico, no tabla de calibración importada.
- **[R20]** ICRU (2005), [Report 74, Patient Dosimetry for X Rays used in Medical Imaging](https://www.icru.org/report/patient-dosimetry-for-x-rays-used-in-medical-imaging-report-74/). Resumen oficial consultado sobre magnitudes y métodos; no se afirma lectura del informe completo.
- **[R21]** ICRP, [Publication 145, Adult Mesh-type Reference Computational Phantoms](https://www.icrp.org/publication.asp?id=ICRP+Publication+145). Ficha oficial: un fantoma volumétrico con regiones delgadas no equivale a nuestra superficie externa XCAT.
- **[R22]** NRC, [EGSnrc](https://nrc-cnrc.github.io/EGSnrc/) y [dose scoring](https://nrc-cnrc.github.io/EGSnrc/doc/pirs898/classEGS__DoseScoring.html). Documentación primaria; código no ejecutado en esta fase.
- **[R23]** [GATE / Geant4](https://opengate.readthedocs.io/en/latest/radiotherapy_and_dosimetry_applications.html), [TOPAS scoring](https://github.com/topasmc/topas-docs/blob/master/parameters/scoring/output.rst). Documentación primaria; herramientas externas propuestas, no integradas ni ejecutadas.
- **[R24]** OECD-NEA, [PENELOPE-2024](https://www.oecd-nea.org/jcms/pl_108003/penelope-2024-a-code-system-for-monte-carlo-simulation-of-electron-and-photon-transport). Transporte acoplado electrón–fotón; ficha técnica oficial consultada, no simulaciones propias.
- **[R25]** Xiang HF et al. (2007). *Build-up and surface dose measurements on phantoms using micro-MOSFET in 6 and 10 MV x-ray beams and comparisons with Monte Carlo calculations*. Medical Physics 34(4):1266–1273. DOI [10.1118/1.2710951](https://doi.org/10.1118/1.2710951); [PMID 17500458](https://pubmed.ncbi.nlm.nih.gov/17500458/). Resumen consultado; evidencia MV, fuera del alcance de nuestro formalismo kV.
- **[R26]** [PySkinDose, repositorio oficial](https://github.com/rvbCMTS/PySkinDose). Se mantiene su mesh alternativo con procedencia/licencia. No se ha comparado contra la ejecución de su motor completo.

Los PDF aportados no se redistribuyen dentro del ZIP. Las citas no autorizan atribuir los resultados publicados a la implementación actual.
