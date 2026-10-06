---
title: "Calculador independiente de dosis en piel"
subtitle: "Memoria técnica y científica del prototipo skindose-verifier 0.5.3"
date: "28 de septiembre de 2026"
lang: es-ES
documentclass: report
classoption: [11pt, a4paper, oneside]
geometry: [margin=25mm]
mainfont: Latin Modern Roman
monofont: DejaVu Sans Mono
colorlinks: true
linkcolor: black
urlcolor: teal
toc: true
toc-depth: 2
numbersections: true
header-includes:
  - \usepackage{float}
  - \usepackage{microtype}
  - \usepackage{fvextra}
  - \usepackage{seqsplit}
  - \sloppy
  - \DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\{\}}
  - \setlength{\emergencystretch}{3em}
  - \setcounter{secnumdepth}{2}
  - \AtBeginDocument{\renewcommand{\tablename}{Tabla}\renewcommand{\figurename}{Figura}}
---

# Resumen {.unnumbered}

Se presenta una memoria del estado implementado de `skindose-verifier` 0.5.3, un prototipo Python para reconstruir la distribución de dosis de entrada sobre una superficie de paciente a partir de eventos de irradiación contenidos en informes DICOM RDSR. El programa combina geometría de arco y mesa, proyección de un campo rectangular divergente, selección de la primera superficie interceptada y factores dosimétricos configurables. Puede trabajar con un plano, un cilindro elíptico o una malla triangular cerrada. La demostración antropomórfica utiliza exclusivamente el contorno corporal externo de un fantoma XCAT de referencia.

La revisión comprende el código ejecutable, los parámetros, los modelos, la documentación de desarrollo, las pruebas y los resultados conservados. Se ha repetido la suite de 54 pruebas, se han realizado comprobaciones analíticas y se han medido convergencia, escalabilidad y memoria. En el ejemplo XCAT de dos eventos, una malla de 177 926 caras produce una PSD puntual de 0,0347763062 Gy. Al reducir la arista máxima de 10 a 5 mm, la PSD cambia un 0,0165 %, mientras que el área irradiada cambia un 0,395 %. Estas cifras cuantifican sensibilidad numérica en una anatomía y un caso determinados; no representan exactitud dosimétrica clínica.

El principal resultado es una plataforma verificable para investigación geométrica y dosimétrica, con mapas por evento, métricas ponderadas por área, propagación paramétrica y visor tridimensional autónomo. Los factores disponibles en la demostración son sintéticos, no existen medidas propias de validación y no se ha ejecutado una comparación externa con transporte Monte Carlo. Por tanto, no puede asignarse una precisión clínica al programa. La siguiente fase debe priorizar caracterización de equipo, registro anatómico, pruebas con RDSR reales y validación experimental independiente.

**Palabras clave:** dosis en piel; radiología intervencionista; RDSR; fantoma triangular; dosis pico; verificación computacional; incertidumbre.

# Introducción

## Contexto clínico y tecnológico

Los procedimientos guiados por fluoroscopia permiten diagnosticar y tratar enfermedades mediante acceso mínimamente invasivo. La imagen radiológica sirve para dirigir catéteres, posicionar dispositivos y comprobar el resultado de una intervención. En cardiología, neurorradiología y radiología vascular pueden combinarse periodos prolongados de fluoroscopia con adquisiciones de mayor calidad de imagen. La complejidad técnica del procedimiento, el espesor atravesado y la necesidad de determinadas proyecciones condicionan la irradiación. El beneficio clínico de la intervención no elimina la necesidad de conocer y optimizar su exposición radiológica [1, 23].

La irradiación cutánea no es uniforme. Cada posición del arco deposita dosis principalmente en una región de entrada; los desplazamientos del paciente y de la mesa modifican esa región. Cuando varias proyecciones se superponen, una zona puede recibir una dosis considerablemente mayor que el resto de la piel. La distribución acumulada, y no solamente un indicador global de exposición, es la información relevante para localizar regiones de mayor dosis. El riesgo de reacciones tisulares depende además de la magnitud y distribución de la dosis, del tiempo entre exposiciones y de características individuales. El calculador no incorpora un modelo radiobiológico ni predice una lesión [2, 3, 24].

En los equipos de intervencionismo se dispone habitualmente de indicadores como el kerma en aire en el punto de referencia intervencionista, $K_{a,r}$, y el producto kerma-área, $P_{KA}$. Son magnitudes útiles, pero ninguna especifica por sí sola el máximo cutáneo. El punto de referencia puede encontrarse separado de la piel; entre él y la superficie cambian distancia, tamaño de campo, interposición de mesa y retrodispersión. Además, dos procedimientos con idéntico kerma acumulado pueden producir mapas muy distintos si sus campos coinciden o se distribuyen entre varias regiones [1, 2, 5].

La información necesaria para reconstruir cada exposición puede almacenarse en un RDSR. A diferencia de un informe con un único total, su estructura permite describir eventos de irradiación y sus parámetros. Sin embargo, disponer de un objeto DICOM no garantiza que contenga toda la información requerida por un algoritmo: pueden faltar el registro del paciente, las dimensiones completas del campo o convenciones inequívocas de posición. El diseño del lector debe distinguir ausencia, incompatibilidad y dato conocido, en lugar de suplir silenciosamente la información que determina el resultado [4, 8, 15].

## Métodos de estimación y medida

La dosimetría directa permite obtener una referencia experimental. La película radiocrómica ofrece distribución espacial de alta resolución, pero requiere calibración, protocolo de lectura y control de la respuesta a la calidad del haz. Los detectores TLD, OSLD o de vidrio pueden distribuirse sobre la superficie, aunque un conjunto finito de puntos puede no coincidir con el verdadero máximo. Las cámaras de ionización aportan trazabilidad metrológica, con limitaciones de volumen sensible, posicionamiento y resolución espacial. La interpretación debe incluir la perturbación introducida por el detector y su correspondencia con la magnitud calculada [2, 4, 7].

Una alternativa indirecta transforma los indicadores del equipo mediante factores de calibración y correcciones geométricas. La estimación más simple usa supuestos globales acerca de distancia y superposición; una reconstrucción por evento permite localizar la huella de cada campo y sumar dosis espacialmente. Su ventaja es aprovechar datos rutinarios, con coste computacional moderado. Su limitación es que la información omitida en los datos de entrada no se recupera por aumentar la sofisticación de la visualización [3, 5, 26].

Los métodos de transporte Monte Carlo representan las interacciones de los fotones y, según el código y la aproximación elegida, el transporte de electrones secundarios. Pueden tratar materiales, dispersión y dosis en profundidad. Exigen definir espectro, geometría, materiales y normalización, y controlar la estadística en regiones pequeñas. Takata y colaboradores utilizaron simulación acelerada por GPU con geometría derivada de CT; ello ilustra una familia de métodos físicamente distinta del trazado geométrico superficial implementado aquí [7].

La tabla siguiente sitúa el proyecto respecto a esas alternativas. La comparación se refiere al alcance de los métodos, no a una clasificación universal de exactitud.

| Método | Información que aporta | Limitación dominante |
|:--|:--|:--|
| Cámara o detectores puntuales | Dosis o kerma en posiciones muestreadas | Muestreo espacial y respuesta del detector |
| Película calibrada | Mapa sobre la superficie cubierta | Calibración, lectura y registro |
| Estimación global desde indicadores | Orden de magnitud del máximo bajo supuestos | No reconstruye superposición real |
| Reconstrucción superficial por evento | Distribución de entrada y acumulación | Geometría, datos y factores aproximados |
| Transporte Monte Carlo | Distribución dependiente de materiales y transporte | Modelo de fuente, coste y estadística |

## Motivación y problema abordado

El proyecto nace de la necesidad de disponer de un calculador transparente que permita inspeccionar cómo se transforma un evento RDSR en una contribución de dosis y por qué aparece un máximo en una posición concreta. La independencia buscada es computacional y metodológica: el usuario puede revisar código, configuración, factores y resultados, sin depender de una implementación comercial cerrada. Esto no implica independencia metrológica respecto a los datos del equipo. La rama principal utiliza el $K_{a,r}$ informado por el propio equipo; un error de ese indicador puede propagarse a todos los puntos calculados.

Una segunda motivación es separar problemas que suelen confundirse. Un campo mal situado puede deberse al signo de un ángulo, a un origen de mesa incorrecto o a un registro anatómico inadecuado, incluso cuando la ecuación dosimétrica es correcta. Una distribución visualmente atractiva puede proceder de una malla insuficiente. Un resultado reproducible puede mantener un sesgo físico. El programa y esta memoria organizan esas cuestiones en entradas, geometría, física, discretización, representación y evidencia experimental.

La revisión de Malchair et al. identificó 22 herramientas con diferencias importantes de capacidad y evaluación. La validación en fantoma no garantizaba un comportamiento comparable en pacientes ni compatibilidad con todos los fabricantes [5]. El estudio VERIDIC de Dabin et al., sobre diez productos y cuatro sistemas radiológicos, encontró dificultades especialmente en irradiaciones laterales y en la estandarización de los datos [27]. Estos resultados justifican una estrategia de evaluación por componentes y por equipo, sin heredar una precisión de la literatura.

La revisión de Wunderle y Balter de 2026 mantiene vigente esta cautela: los avances en indicadores y RDSR no eliminan las incertidumbres de geometría, hardware y mapeo. Su planteamiento respalda documentar los supuestos del sistema y reforzar los datos de entrada antes de atribuir valor clínico a una automatización [28]. La publicación IAEA de 2025 sobre garantía de calidad y optimización en procedimientos guiados por fluoroscopia aporta un marco actual de QA, aunque esta memoria no demuestra conformidad del software con ese documento [29].

## Objetivos y alcance

El objetivo general es desarrollar y evaluar un prototipo trazable de estimación de dosis de entrada en piel, capaz de reconstruir y acumular eventos de intervencionismo sobre una superficie de referencia. Los objetivos específicos son interpretar un subconjunto explícito de RDSR; reconstruir fuente, campo, mesa y paciente; calcular la superficie visible; aplicar factores documentados; conservar las contribuciones individuales; cuantificar máximos y áreas; explorar sensibilidad e incertidumbre; y producir resultados numéricos y gráficos auditables.

El alcance actual es retrospectivo y de investigación. Incluye configuraciones monoplano y capacidad de parametrización por plano, con pruebas sintéticas biplano. No incluye adquisición en tiempo real desde un equipo, conexión PACS, segmentación CT, registro automático, deformación anatómica, dosimetría de órganos ni transporte en volumen. Tampoco interpreta RTPLAN ni RTDOSE de radioterapia. La palabra «piel» designa aquí una superficie de entrada, no una capa histológica con espesor y composición representados.

El documento de diseño `spec.txt` plantea además caracterización experimental y comparación con sistemas comerciales. Esas fases están previstas y parcialmente soportadas por interfaces de comparación, pero no ejecutadas con datos independientes. El objetivo de disponer de un prototipo computacional se ha alcanzado; el de disponer de un verificador clínico caracterizado permanece abierto.

## Método de revisión y nivel de evidencia

Se revisaron los 22 módulos Python del paquete, ambos archivos de pruebas, las configuraciones, la plantilla HTML, los lanzadores, los scripts de estudio, las licencias de mallas y los documentos de evolución. Se inventariaron los archivos con tamaño y SHA-256. Los resultados previos se contrastaron con nuevas ejecuciones; el código del motor se conservó sin modificaciones durante esta auditoría. No existe un historial Git en la copia disponible: las decisiones históricas se reconstruyen desde código, documentación y artefactos, no desde una secuencia de commits.

La revisión bibliográfica combina los once PDF aportados, búsqueda por títulos y conceptos en PubMed y editoriales, y fuentes oficiales de AAPM, IAEA, ICRP, ICRU, NIST, DICOM y BIPM. Se buscaron trabajos recientes con términos sobre *skin dose mapping*, *peak skin dose*, *interventional fluoroscopy*, validación y años 2025–2026, y se siguieron referencias del formalismo. El anexo documenta acceso y exclusiones. Es una revisión técnica amplia y dirigida al modelo; no se presenta como revisión sistemática con protocolo registrado, cribado doble o garantía de recuperar todas las publicaciones.

En adelante, «resultado propio» indica una salida obtenida del código auditado; «verificación» indica concordancia con una ecuación o invariante; «validación externa» requiere referencia física independiente. Cuando una justificación no está registrada como decisión histórica, se presenta como interpretación técnica. La bibliografía fundamenta métodos y limitaciones, pero no convierte resultados ajenos en validación del proyecto.

# Fundamentos físicos y teóricos

## Magnitudes y unidades

La fluencia de partículas $\Phi=dN/da$ cuenta partículas que atraviesan una esfera de sección $da$ por unidad de área. La fluencia energética $\Psi=dR/da$ pondera cada partícula por su energía. Para fotones monoenergéticos, $\Psi=E\Phi$; para un espectro, se integra la fluencia diferencial ponderada por energía. Esta distinción importa porque dos haces con igual número de fotones pueden transferir distinta energía a la materia [22].

El kerma es la energía cinética inicial transferida por partículas no cargadas a partículas cargadas por unidad de masa. La dosis absorbida es la energía media impartida a la materia por unidad de masa. Ambas se expresan en gray, $1\,\mathrm{Gy}=1\,\mathrm{J\,kg^{-1}}$, pero no son la misma magnitud. La energía puede transferirse en un lugar y depositarse en otro debido al movimiento de electrones secundarios. Para un haz monoenergético:

$$K=\Psi\left(\frac{\mu_{\mathrm{tr}}}{\rho}\right),\qquad
K_{\mathrm{col}}=\Psi\left(\frac{\mu_{\mathrm{en}}}{\rho}\right).
\tag{1}$$

Los coeficientes másicos de transferencia y absorción energética tienen unidades de área por masa. El segundo descuenta la fracción de energía cinética que se pierde radiativamente. Bajo equilibrio de partículas cargadas, la dosis puede relacionarse con el kerma de colisión. En una interfaz aire-tejido, esa equivalencia no debe imponerse sin considerar condiciones de medida y transporte [19, 20, 22].

El producto kerma-área integra el kerma sobre una sección transversal del haz. Para un campo uniforme, $P_{KA}=K_a A$. Puede conservarse aproximadamente al cambiar la distancia en aire, porque el kerma disminuye como el inverso del cuadrado y el área aumenta con el cuadrado. Su unidad habitual es Gy·cm², mientras que el programa almacena Gy·m². La conversión es $1\,\mathrm{Gy\,cm^2}=10^{-4}\,\mathrm{Gy\,m^2}$. No equivale a la integral de dosis sobre la piel, cuya superficie puede ser oblicua y cuya dosis incluye otras correcciones [2, 19].

El núcleo trabaja con distancias en mm, ángulos en grados, áreas de muestreo en mm² y dosis en Gy. Las áreas se convierten a cm² dividiendo por 100. Las figuras del visor muestran mGy, mediante multiplicación por 1000. Estas conversiones son operaciones explícitas; el lector rechaza unidades desconocidas. `domain.py` define los objetos físicos y `rdsr.py` aplica las conversiones por campo.

## Interacción de los fotones con la materia

En el intervalo diagnóstico, la absorción fotoeléctrica y la dispersión Compton condicionan la transferencia energética y la atenuación; la dispersión coherente cambia la dirección con transferencia energética despreciable. La importancia relativa depende de energía y composición. Un espectro de rayos X contiene una distribución de energías, no fotones de energía igual al kVp. El kVp limita la energía máxima del espectro y, junto con filtración, ánodo y geometría, caracteriza su calidad [11, 22].

La ley exponencial de atenuación se obtiene considerando que la probabilidad de retirada del haz primario en un espesor $dx$ es proporcional a $\mu dx$. Para un medio homogéneo y un haz monoenergético estrecho:

$$\frac{d\Phi}{dx}=-\mu\Phi,\qquad
\Phi(x)=\Phi(0)e^{-\mu x}.
\tag{2}$$

En un haz polienergético, cada componente se atenúa de forma distinta. El espectro transmitido se endurece y una única constante efectiva puede no conservar su validez al cambiar el espesor. La ley exponencial describe la componente primaria; no incorpora por sí sola la radiación dispersa que llega al punto de evaluación [11, 22].

Esta física justifica la necesidad de factores dependientes de calidad y geometría, pero el programa no evalúa la ecuación (2), no consulta una base de coeficientes y no calcula el espesor atravesado en tejido. `patient_thickness_mm` es un argumento escalar opcional para las tablas de factores. No es una integral de trayecto ni una densidad obtenida del fantoma. La distinción es decisiva: disponer de una malla cerrada no significa simular atenuación dentro del cuerpo.

## Aproximación superficial y conversión aire a tejido

El formalismo implementado aproxima la dosis de entrada a partir del kerma incidente y de factores que representan efectos no transportados explícitamente. Para un espectro local, una conversión aire-medio puede expresarse conceptualmente como:

$$f_{m,a}=
\frac{\int \Psi(E)(\mu_{\mathrm{en}}/\rho)_m(E)\,dE}
{\int \Psi(E)(\mu_{\mathrm{en}}/\rho)_a(E)\,dE}.
\tag{3}$$

Aquí $m$ representa el medio de interés y $a$ el aire. La distribución espectral utilizada debe corresponder a la definición dosimétrica del factor, incluyendo o separando de forma coherente la retrodispersión. El trabajo de Benmakhlouf et al. proporciona un fundamento para calcular factores de retrodispersión y cocientes de absorción mediante espectros y Monte Carlo; sus tablas no están incorporadas automáticamente a este proyecto [25].

El código multiplica por `f_skin`, un factor suministrado por configuración. Permite etiquetar el resultado como dosis en piel o en agua mediante `dose_medium`, pero cambiar esa etiqueta no cambia materiales ni recalcula coeficientes. El usuario debe proporcionar un factor compatible con el medio elegido. En la demostración `f_skin=1,05` es una constante sintética; no procede de una caracterización del equipo ni de una integración espectral ejecutada por el programa.

Se escoge una superficie porque el problema inicial es localizar y acumular dosis de entrada con datos RDSR, sin necesitar una anatomía volumétrica. Frente a una discretización voxelizada, reduce las incógnitas a elementos cutáneos y facilita localizar huellas y solapamientos. El coste físico es excluir dosis en profundidad, heterogeneidades y salida. Esta justificación es compatible con la arquitectura y la documentación del proyecto; no demuestra que una superficie sea suficiente para todos los usos clínicos [2–4].

## Fuente puntual y ley del inverso del cuadrado

Si la emisión en un pequeño ángulo sólido se propaga sin atenuación apreciable, la misma energía se reparte sobre áreas proporcionales a $d^2$. Para el evento $e$, el kerma en un punto de superficie $p$ se obtiene trasladando el kerma de referencia desde una distancia $d_{r,e}$:

$$K_{a,e}(p)=K_{a,r,e}\left(\frac{d_{r,e}}{\|p-s_e\|}\right)^2.
\tag{4}$$

La posición $s_e$ es el foco. Se usa la distancia euclídea foco-punto, no la distancia axial ni una distancia común para toda la cara posterior. El cociente es adimensional y el resultado mantiene Gy. En `engine._dose`, esta expresión aparece como `ka_r_gy * (irp_distance_mm / distance)**2`.

Las hipótesis son una fuente puntual, propagación rectilínea y ausencia de atenuación explícita en aire. No se modela el tamaño finito del foco ni su penumbra. La aproximación proporciona una dependencia geométrica clara y verificable, y evita exigir información de fuente que el RDSR no aporta. Una extensión a foco finito requeriría integrar contribuciones de varias posiciones y validar bordes de campo, con mayor coste por evento.

Para pequeños errores de distancia, $\delta D/D\simeq-2\delta d/d$ si el resto permanece fijo. A 600 mm, un error axial de 5 mm implica aproximadamente un 1,67 % de variación por distancia. Es una derivación local, no una incertidumbre medida del proyecto: cerca de bordes o de campos superpuestos, el cambio de máscara puede dominar y producir una respuesta no lineal.

## Geometría divergente del campo

Sea $b_e$ la dirección unitaria del eje del haz y $u_e,v_e$ sus ejes transversales. Para $\Delta=p-s_e$, la profundidad axial es $t=\Delta\cdot b_e$. Si el campo tiene dimensiones $w_e,h_e$ a distancia de referencia $d_{f,e}$, el punto pertenece al cono piramidal rectangular cuando:

$$0<t\leq SID_e,\qquad
|\Delta\cdot u_e|\leq\frac{w_e t}{2d_{f,e}},\qquad
|\Delta\cdot v_e|\leq\frac{h_e t}{2d_{f,e}}.
\tag{5}$$

La referencia del campo se declara como detector o isocentro. El área utilizada como argumento de las tablas dosimétricas es la sección perpendicular al eje a profundidad $t$:

$$A_\perp(p)=w_e h_e\left(\frac{t}{d_{f,e}}\right)^2/100
\quad [\mathrm{cm^2}].\tag{6}$$

En la versión actual se diferencia correctamente $t$ de $d=\|\Delta\|$: la primera determina divergencia y área perpendicular, mientras que la segunda interviene en el inverso del cuadrado. Usar $d$ para la ecuación (6) aumentaría artificialmente el campo tabulado fuera del eje. La documentación de evolución registra esa corrección respecto a v0.4 y existe una prueba analítica específica.

El área perpendicular no es la huella sobre una superficie oblicua. Para un pequeño elemento plano, la proyección incluye un factor geométrico de coseno, pero eso no autoriza multiplicar la dosis por ese coseno. La fluencia energética se define sobre un área normal a la dirección de transporte; la orientación cutánea cambia entrada y distribución espacial, y las correcciones físicas angulares deben derivarse de medidas o de un modelo dosimétrico. El código calcula el ángulo y lo ofrece a las tablas, sin imponer una corrección angular universal [2, 11].

## Mesa y colchoneta

La mesa se representa como un rectángulo horizontal acotado en el sistema de sala, desplazado con sus coordenadas. Se calcula si el segmento foco-punto cruza su plano antes de llegar al paciente. Solo entonces se aplica transmisión. La incidencia sobre la mesa se obtiene para cada rayo divergente; no se identifica con un único ángulo del gantry.

La transmisión puede configurarse como `table_total` o como producto de `table_normal` y `table_angular`. El programa prohíbe suministrar simultáneamente ambas parametrizaciones para evitar duplicación. Sin embargo, no puede detectar que una tabla experimental de transmisión ya incluya retrodispersión o conversión de medio. La coherencia física de la definición sigue siendo responsabilidad de la caracterización [2, 4].

Un factor tabulado efectivo permite representar mesa y colchoneta sin modelar sus capas; resulta más sencillo que transportar fotones por materiales detallados. Su coste es que no describe raíles, accesorios, inhomogeneidades laterales ni cambios espectrales de forma explícita. La aproximación $T(\theta)=\exp[-\mu x/|\cos\theta|]$ podría proponerse para una lámina homogénea monoenergética, pero no está implementada y no debe atribuirse al programa. En la configuración sintética, $T=0,8$ para los rayos que cruzan la mesa, sin dependencia angular real.

## Retrodispersión y región superficial

La radiación que interacciona en el paciente puede regresar hacia la superficie de entrada. Un factor de retrodispersión resume el cociente entre una magnitud con fantoma y la correspondiente magnitud en aire, bajo una convención experimental o calculada concreta. Depende de calidad, tamaño de campo, material, espesor y condiciones geométricas. El formalismo multiplicativo es útil cuando esas dependencias pueden representarse mediante datos calibrados [2, 25].

El programa permite que `bsf` dependa de kVp, Cu, área perpendicular, espesor o ángulo cutáneo. No estima esas dependencias a partir de composición o transporte. Un BSF constante de 1,25 solo demuestra cómo se aplica un factor, no reproduce automáticamente los datos de una publicación. Tampoco distribuye dosis dispersa fuera de la máscara primaria; allí la contribución del evento permanece exactamente cero.

En la región superficial, equilibrio electrónico, profundidad efectiva de evaluación y respuesta del detector condicionan la relación entre kerma y dosis. El término *build-up* describe el establecimiento de la población de partículas cargadas y de la dosis con profundidad; no debe importarse una curva de un haz de megavoltaje a un problema de kilovoltaje. El programa no contiene una dimensión de profundidad, transporte electrónico ni una capa de piel. Su valor no puede identificarse automáticamente con dosis a 0,07 mm, con una media en 1 mm de tejido o con el resultado de un dosímetro concreto [2, 11, 20].

Esta limitación no se corrige aumentando los triángulos. El refinamiento mejora el muestreo sobre la misma frontera geométrica, pero no añade materiales ni espesor. Para resolver dosis superficial a profundidad definida sería necesario modificar la formulación física, establecer el medio y usar una referencia experimental o de transporte compatible.

## Ecuación de dosis implementada

La contribución de un evento a un elemento superficial queda expresada por:

$$D_e(p)=M_e(p)\,K_{a,r,e}
\left(\frac{d_{r,e}}{\|p-s_e\|}\right)^2
CF_e\,P_e(p)\,T_e(p)\,BSF_e(p)\,f_{m,a,e}(p).
\tag{7}$$

La máscara $M_e$ vale uno para puntos orientados hacia el foco, dentro del campo y visibles como primera entrada; vale cero en el resto. $CF$ corrige el indicador del equipo; $P$ es un perfil relativo opcional; $T$ es transmisión de mesa o uno cuando no se cruza; $BSF$ representa retrodispersión y $f$ la conversión de medio. Todos los factores son adimensionales. Esta estructura conecta el formalismo de la literatura con `engine._dose` [2–4].

El programa distingue tres etapas. `free_air` aplica exclusivamente máscara e inverso del cuadrado. `calibrated_air` añade calibración, perfil y mesa. `skin_dose` añade retrodispersión y conversión. La etiqueta de salida cambia con la magnitud. Ejecutar `free_air` no produce dosis absorbida, aunque el resumen conserve por razones de estructura el nombre de campo `psd_point_gy`; este nombre debe leerse junto con `quantity` y `stage`.

En la demostración, $CF=1$, $P=1$, $T=0,8$, $BSF=1,25$ y $f=1,05$. Cuando el rayo cruza la mesa, el producto de factores vale 1,05. La cancelación numérica entre $0,8$ y $1,25$ es una elección sintética, no una propiedad física universal. Con esas constantes, el cálculo verifica geometría y multiplicación; no demuestra respuesta energética del equipo.

## Dependencias tabuladas y perfil del haz

Cada factor exige identificador y procedencia, y puede ser escalar o una tabla sobre ejes crecientes. `RegularGridInterpolator` aplica interpolación multilineal con error fuera del dominio. Entre dos nodos de un eje, $f(x)=(1-\lambda)f_0+\lambda f_1$, donde $\lambda=(x-x_0)/(x_1-x_0)$. En varias dimensiones se combinan los pesos de las esquinas de la celda. Se evita introducir oscilaciones de interpolaciones de orden alto sin datos que las justifiquen [18].

El bloqueo de extrapolación protege tablas, pero no limita el dominio de un factor escalar. Una constante positiva se acepta para cualquier kVp o espesor válidos sintácticamente. Es una carencia de control de aplicabilidad: futuras configuraciones deberían declarar rangos incluso cuando usen valores constantes.

El perfil opcional `beam_profile` utiliza coordenadas normalizadas del campo, calidad y otros ejes admitidos. Puede representar un perfil medido y parte del efecto talón. El software no verifica que su valor central sea uno ni su coherencia con el plano y la normalización de $K_{a,r}$. Un perfil mal normalizado cambia la dosis global; la procedencia textual no sustituye una prueba metrológica.

## Acumulación y métricas

El mapa acumulado es $D(p)=\sum_e D_e(p)$ en coordenadas fijas del paciente. La PSD puntual discreta es $\max_i D(p_i)$. No se suman los máximos de eventos separados, porque pueden aparecer en posiciones diferentes. Para dosis no negativas, la suma de máximos individuales constituye una cota superior del máximo acumulado; la igualdad exige coincidencia espacial adecuada.

Para caras de área $A_i$, se aproxima la integral superficial por $\sum_i A_iD_i$. El histograma dosis-área acumulado evalúa $A(D\geq d)=\sum_i A_i\,\mathbf{1}(D_i\geq d)$. La media es ponderada por área, no por número de triángulos. $D_{2\%}$, $D_{5\%}$ y $D_{10\%}$ son valores superados por al menos esas fracciones del área declarada, con inversión escalonada de la distribución. Se calculan por separado sobre todo el modelo y sobre la superficie irradiada.

La media local denominada de 1 cm² busca vecindades aproximadamente circulares. En una malla usa caminos entre centroides de caras adyacentes, evitando saltar entre brazo y torso próximos en espacio. Su radio objetivo es $\sqrt{100/\pi}$ mm y admite áreas de 85 a 125 mm². Por tanto, no es una integración exacta de 100 mm² ni necesariamente la huella de un dosímetro. Puede devolver ausencia de resultado. En convergencia de mallas se omite además deliberadamente mediante `detailed=False`; un `null` en ese informe no prueba incapacidad de resolución.

Los hotspots agrupan caras conectadas con dosis al menos del 90 % del máximo acumulado. Son métricas descriptivas propias de la herramienta, no umbrales de toxicidad. La integral dosis-área cutánea tiene unidades Gy·cm², pero no es KAP ni energía absorbida: esta última exigiría integrar dosis con masa, ausente en una superficie sin espesor.

## Incertidumbre y alternativas de mayor rigor

El resumen básico combina cuadráticamente las incertidumbres relativas disponibles, suponiendo independencia, y declara componentes no cuantificadas. Es incompleto incluso con todos los factores informados, porque no propaga toda la geometría ni la discretización. Cuando no hay ninguna componente conocida devuelve cero como suma vacía; la bandera `complete_budget=False` impide interpretarlo como incertidumbre total nula, pero la presentación numérica puede inducir a error.

El módulo paramétrico toma muestras normales, uniformes o lognormales mediante LHS o Monte Carlo ordinario. Cada realización recalcula el mapa. Los parámetros son independientes entre sí; el sesgo de cada parámetro se comparte entre todos los eventos. Se calculan media y varianza por punto mediante actualización incremental, máximos por realización y cuantiles. Los principios de propagación de distribuciones son coherentes con JCGM 101, sin que ello certifique el presupuesto ni las distribuciones elegidas [14].

Este Monte Carlo de parámetros no transporta fotones. La alternativa de transporte permitiría tratar heterogeneidad y dispersión, pero necesitaría datos adicionales y estadística suficiente. Han et al. proponen combinar simulación Monte Carlo de menor estadística con reducción iterativa de ruido para estimar el máximo; ese filtrado responde al ruido estadístico del transporte y no justifica suavizar el mapa determinista de este calculador [12]. Un enfoque híbrido futuro podría tabular correcciones de transporte, siempre con evaluación fuera de los casos usados para ajustar las tablas.

# Desarrollo y procedimiento computacional

## Arquitectura y flujo de información

El paquete separa objetos de entrada, configuración, superficie, proyección geométrica, dosis, métricas y exportación. Python coordina el flujo; NumPy realiza operaciones vectorizadas, SciPy interpola y maneja grafos, pydicom interpreta SR, Trimesh y Rtree resuelven geometría, Matplotlib genera figuras y Plotly produce el visor autónomo. No hay un servidor web de cálculo ni una GUI para editar calibraciones: la interfaz principal es una CLI y archivos JSON.

El flujo comienza al leer RDSR y configuración. Se comprueba fabricante/modelo, se normalizan unidades y se completan únicamente distancias ausentes desde una configuración explícita. Se prepara una superficie local fija. Cada evento se transforma al sistema de sala, se proyecta el campo, se resuelve visibilidad y se evalúan factores. Los mapas se acumulan sobre los mismos índices de paciente; las exportaciones conservan tanto la suma como cada contribución.

![Flujo implementado y separación entre cálculo, visualización y evaluación externa. Elaboración propia a partir del código.](figuras/arquitectura.png){width=90%}

La separación permite ensayar geometría sin correcciones dosimétricas, reutilizar la superficie en estudios paramétricos y comprobar cada factor sin depender del visor. Evita mezclar transformaciones visuales con cálculos físicos. Frente a un monolito gráfico, mejora trazabilidad, aunque obliga a mantener contratos consistentes entre módulos y a informar con claridad las unidades.

## Lectura y normalización del RDSR

El lector exige SOP Class `1.2.840.10008.5.1.4.1.1.88.67` y recorre `ContentSequence`. Reconoce conceptos por pareja esquema-código, no por texto libre de `CodeMeaning`. Localiza contenedores de evento `113706`, requiere UID y plano y rechaza UIDs repetidos. Los campos numéricos deben incluir unidades UCUM conocidas. Las lecturas repetidas de un mismo concepto se aceptan solo si sus valores coinciden dentro de tolerancia; no se promedian valores discrepantes.

Se conservan las lecturas originales y marcas de calidad para rastrear datos presentes, ausentes o completados. El material filtrante admitido es un único filtro de cobre identificado; un espesor de aluminio no se interpreta como Cu. Esta política evita resultados aparentemente válidos con una calidad de haz mal identificada, a costa de rechazar RDSR clínicos que requerirían una representación más rica. La filtración inherente del equipo puede estar incorporada implícitamente en una calibración, pero no está reconstruida espectralmente.

El lector distingue coordenadas de mesa isocéntricas y heredadas. La configuración puede seleccionar el conjunto adecuado si ambos aparecen. No infiere mAs a partir de corriente y tiempo porque, en emisión pulsada, sus promedios y convenciones pueden no representar la integral requerida. Los totales acumulados de fluoroscopia y adquisición se comparan con la suma de eventos con tolerancia relativa del 2 % o absoluta de $10^{-7}$ en unidades internas.

La selección por descendientes es un lector de subconjunto, no un validador completo de plantillas DICOM. La jerarquía de totales y planos necesita ensayos específicos con RDSR reales; el recorrido general puede encontrar conceptos repetidos con distinto alcance. El generador sintético y el lector comparten códigos, por lo que su concordancia no demuestra conformidad completa con TID 10001/10003 ni portabilidad multivendedor. Enhanced RDSR dispone de una representación geométrica más explícita, pero no está implementado [15].

## Coordenadas de equipo y paciente

La convención interna toma $x$ hacia la izquierda del paciente, $y$ craneal y $z$ anterior. El haz inicial apunta en $+z$ con el foco bajo la mesa. El eje del haz se obtiene de la matriz de rotación. Para la opción denominada `primary_then_secondary`, el código usa exactamente $R_y(a)R_x(b)$; el nombre no sustituye la interpretación matemática del producto ni la verificación de signos en el equipo real.

El foco se sitúa en $s=c_{iso}-SOD\,b$ y el detector en $s+SID\,b$. La traslación de mesa es $M_{mesa}(q-q_0)$, donde $q$ son coordenadas del evento y $q_0$ el origen caracterizado. Los puntos del paciente se rotan según inclinaciones configuradas y se trasladan con la mesa y con el desplazamiento paciente-mesa. Las normales se rotan sin traslación. En la rama triangular, el foco se transforma de vuelta al marco local para intersectar una malla que permanece fija.

Esta decisión evita reconstruir un árbol espacial por cada posición de mesa. La anatomía permanece en su sistema local, y cambiar la geometría requiere transformar fuente y puntos, no modificar permanentemente el STL. La alternativa de mover todos los vértices y reconstruir el índice espacial incrementaría trabajo. La contrapartida es exigir coherencia exacta entre marcos; una prueba con rotación e inversión global de normales verifica parte de ese contrato, pero no los signos de un fabricante concreto.

## Superficies y procedencia anatómica

El plano se discretiza mediante centros de celdas rectangulares y normales verticales. La demostración tiene 360 × 360 mm, paso objetivo de 5 mm y 72 × 72 puntos. El cilindro elíptico usa $x=a\cos\theta$, $z=b\sin\theta$ y una coordenada longitudinal uniforme. Las normales se obtienen del gradiente de la elipse; las áreas aproximan $\sqrt{a^2\sin^2\theta+b^2\cos^2\theta}\,d\theta\,dy$. Es un cilindro abierto: se rechazan focos fuera de sus extremos en el sistema inclinado del paciente.

La importación triangular acepta STL, OBJ y PLY. Exige una sola superficie conectada, cerrada y con orientación coherente. Declara unidades mm, cm o m, porque STL no las codifica. La matriz de ejes debe ser ortonormal y el origen finito. Una orientación global invertida se corrige si el volumen es negativo; se rechazan caras degeneradas. Estos controles no detectan todas las autointersecciones, limitación que se declara en `mesh_quality`.

El modelo principal es `WholeBody_Body.stl`, del repositorio de Ben Auer solicitado durante el desarrollo [16]. El archivo contiene 12 148 caras y representa el contorno externo del fantoma masculino mesh50 XCAT. Su SHA-256 es `2edd0de03344a1081f2a6567bbcca87b965fbbadd822756e05830cd570cdcb3b`. La revisión de origen registrada es `1f20b67fd7fa700d27a3dc88ee9ece436658eee0`. Se conserva licencia y procedencia. No se importaron órganos internos.

La demostración aplica $x=X$, $y=Z-500$ mm y $z=-Y+40$ mm para ubicar aproximadamente el tórax respecto al campo. Esta colocación es un supuesto sintético, no un registro anatómico. Las dimensiones generales del contorno son aproximadamente 658 × 280 × 1764 mm antes de reordenar ejes. El modelo alternativo procedente de PySkinDose contiene 26 756 triángulos y requiere escala de cm a mm. No se ha ejecutado el motor de PySkinDose para compararlo con el presente calculador.

La anatomía de referencia permite ensayar oclusión y proyecciones sobre un cuerpo plausible. Frente a un plano, evita algunas inconsistencias en vistas laterales. No incorpora la forma, postura, peso, apoyo sobre la mesa o respiración de un paciente concreto. Los fantomas de referencia de tipo mesh de ICRP son otra familia de modelos, con regiones anatómicas específicas; no son equivalentes a este STL externo [21].

## Discretización y refinamiento real

Cada cara se evalúa en su centroide y su dosis se considera representativa de toda el área. Esta regla es económica y permite asociar directamente un escalar a un triángulo del visor. Si un borde de colimación atraviesa una cara, el algoritmo asigna toda su área según el centroide. Se introduce así error de cuadratura y aliasing de huella; un campo muy estrecho podría no contener ningún centroide aun cruzando la superficie.

El parámetro `resolution_mm` no modifica un STL. La resolución triangular se cambia mediante `mesh_max_edge_mm`. El refinamiento biseca aristas largas compartidas y genera entre dos y cuatro subcaras según el patrón. Al compartir puntos medios evita uniones no conformes. Puede limitarse a una caja regional; las caras cuyos límites intersectan esa caja impulsan refinamiento, y los vecinos comparten aristas subdivididas.

Se limita el número de caras antes de asignar la nueva malla y el proceso dispone de un máximo de 24 iteraciones. Se comprueba que área, volumen, cierre y orientación se conservan, con tolerancias numéricas. El método añade puntos sobre la geometría original; no suaviza ni reconstruye un contorno anatómico más exacto. El coste crece aproximadamente con el inverso del cuadrado del tamaño de elemento cuando la superficie queda homogéneamente refinada.

Frente a evaluar varios puntos baricéntricos por cara sin cambiar topología, el refinamiento ofrece también una representación gráfica y un grafo más finos. Sin embargo, aumenta memoria e intersecciones. El supersampling podría mejorar integrales sin multiplicar toda la estructura del visor; queda pendiente estudiar cómo representar en una misma cara una dosis no uniforme y conservar máximos. El refinamiento actual es geométrico, no adaptativo según un estimador de error.

## Primera entrada y oclusión

Para cada centroide se comprueban campo, profundidad y orientación $n\cdot(s-p)>10^{-9}$. En una superficie convexa esa selección resuelve gran parte de la entrada; en anatomías con concavidades o extremidades, una cara orientada hacia el foco puede estar oculta por otra. La rama mesh lanza rayos solo a los candidatos y conserva un punto cuando el primer triángulo intersectado coincide con su propia cara.

El intersector utiliza Trimesh con un índice R-tree de cajas de triángulos. Los rayos normalizados se procesan en bloques de 64. El contador de candidatos registra pares rayo-triángulo entregados al núcleo, incluyendo los que finalmente se rechazan. No es el número de triángulos iluminados ni un conteo de fotones. Tampoco es una implementación propia de BVH o un trazador GPU [17].

Conceptualmente, una intersección resuelve $s+\lambda d=v_0+u(v_1-v_0)+v(v_2-v_0)$, con $\lambda>0$, $u,v\geq0$ y $u+v\leq1$. Entre soluciones válidas se selecciona la de menor recorrido positivo. La implementación delegada filtra candidatos y calcula coordenadas baricéntricas; los tests verifican primera entrada, simetría y sombra brazo-torso. No se atribuye a esos tests una demostración exhaustiva de robustez en todas las degeneraciones geométricas.

La fuente se comprueba respecto al volumen cuando cae en la caja de la malla. Un foco interior genera error. El modelo no calcula la contribución que atraviesa el cuerpo y emerge: la eliminación de caras posteriores es correcta para su aproximación de dosis de entrada, pero excluye dosis de salida y dispersión que sí existirían físicamente.

## Eventos, caché y barridos

La evaluación puede resumirse en el siguiente pseudocódigo, que conserva el orden lógico real:

```text
leer eventos y configuracion; preparar superficie local
para cada evento:
    completar distancias declaradas; comprobar datos y factores
    dividir barrido si la configuracion autoriza tasa uniforme
    para cada muestra geometrica:
        buscar proyeccion exacta en cache local
        si falta: proyectar campo, primera entrada y mesa
        calcular inverso del cuadrado y factores de la etapa
        sumar contribucion en los indices fijos del paciente
    conservar mapa e informe del evento
sumar mapas; localizar maximo; obtener areas y metricas
exportar datos originales y visualizaciones
```

La caché guarda como máximo dos proyecciones. Su clave incluye plano, ángulos, posición de mesa, SID, SOD, distancia al punto de referencia y dimensiones del campo. Es local a una llamada de cálculo, donde superficie y equipo permanecen fijos. No usa tolerancia angular ni agrupa geometrías parecidas: solo reutiliza igualdad exacta. Los factores y el kerma se evalúan nuevamente en cada evento.

Esta elección conserva resultados y limita memoria adicional. Una caché mayor aceleraría secuencias que alternen más geometrías a costa de almacenar varias matrices por superficie. Agrupar eventos próximos mediante tolerancias podría reducir coste, pero cambiaría el modelo y exigiría un estudio de error en bordes y máximos. En esta versión no hay paralelización explícita por hilos, procesos o GPU.

Los eventos con ángulos inicial y final diferentes se bloquean por defecto. Al activar `allow_uniform_sweep`, se reparten kerma, KAP y mAs entre `sweep_samples`, evaluando ángulos en los puntos medios de intervalos iguales. Ello asume simultáneamente evolución angular y aporte de dosis compatibles con ese reparto. No se usa un perfil temporal medido, ni se interpola movimiento de mesa dentro del evento. Conservar kerma total no demuestra que el mapa de un barrido real esté correctamente distribuido.

## Informes, visor y experiencia de uso

La salida contiene `resultados.json`, `mapa_dosis.npz`, tabla de eventos, imágenes 2D y 3D, histograma dosis-área y HTML interactivo. El NPZ conserva puntos, normales, áreas, caras, vértices y mapas individuales. El JSON incluye magnitud, advertencias, factores, tiempos, hashes de RDSR/configuración/malla y versiones de dependencias. Esa información permite atribuir un resultado a entradas concretas.

El visor utiliza Plotly con su biblioteca incluida en el HTML. Permite rotación, zoom, perspectivas, selección de evento, representación esquemática de fuente-isocentro-detector, cortes por tres ejes, inspección de caras y hotspots. Los planos recortan polígonos únicamente para mostrarlos. No generan cortes de una distribución volumétrica, no recalculan dosis y no modifican las métricas acumuladas.

La escala cromática se fija al máximo acumulado para comparar eventos. El suavizado opcional calcula colores de vértice a partir de dosis de caras ponderadas por área; no altera las dosis originales. El inspector sigue mostrando datos de centroide. La discretización física debe juzgarse con esos datos, no con la continuidad aparente de color. El arco de gantry es esquemático, no un modelo mecánico del equipo.

El lanzador `EJECUTAR_INVESTIGACION.bat` ejecuta pruebas, demo XCAT refinada, convergencia, benchmark, propagación y escalabilidad mediante `reproduce_research.py`. No procesa automáticamente el ejemplo intervencionista de 24 eventos. Ese caso tiene su propio `PROBAR_EJEMPLO.bat`. Sustituir el lanzador de investigación por una llamada `run` haría un cálculo del ejemplo, pero perdería el estudio completo: la distinción corrige una posible confusión de uso.

Los informes omiten nombre e identificador del paciente, pero conservan UIDs de evento y otros metadatos técnicos. Por tanto, «anonimizado» en comentarios del código no equivale a una desidentificación DICOM completa certificada. Antes de compartir casos reales debe evaluarse el contenido exportado y la trazabilidad que se conserva.

## Trazabilidad de ecuaciones y decisiones

| Elemento | Implementación principal | Evidencia disponible |
|:--|:--|:--|
| Unidades y eventos | `rdsr.numeric`, `parse_event` | Conversiones, ausencias, UID duplicado |
| Fuente y marco | `geometry.beam_from_event`, `table_shift` | Signos sintéticos y traslaciones |
| Campo y área axial | `geometry.project` | Plano y solución off-axis |
| Primera entrada | `raytrace.first_hits` | Esfera, cilindro y autosombra |
| Ecuación (7) | `engine._dose` | Factores escalares y tablas sintéticas |
| Suma espacial | `engine.calculate` | Igualdad suma de mapas y orden |
| Refinamiento | `mesh_processing.refine_mesh` | Área, volumen, cierre y convergencia |
| DAH y percentiles | `metrics.surface_metrics` | Áreas desiguales y monotonicidad |
| Incertidumbre | `uncertainty.propagate_parameters` | Semilla, escala analítica y varianza nula |
| Comparación externa | `validation` | Interfaces probadas; medidas reales ausentes |

# Resultados, validación y discusión

## Evaluación del algoritmo y rendimiento computacional

### Protocolo reproducible

Se ejecutó la versión 0.5.3 bajo Linux y Python 3.12.14, con NumPy 2.3.5, SciPy 1.17.0, Trimesh 5.1.0, Rtree 1.4.1, pydicom 3.0.2, Plotly 6.9.0 y Matplotlib 3.10.8. El entorno virtualizado informa procesador AMD EPYC 9V74 y nueve CPU lógicas visibles. No se reserva uso exclusivo del procesador. Se fijaron `OPENBLAS_NUM_THREADS=1` y `OMP_NUM_THREADS=1` para las mediciones; los tiempos no constituyen una garantía para Windows o para el ordenador del usuario.

La suite de 54 pruebas terminó sin fallos y la compilación sintáctica con `compileall` devolvió código cero. Los archivos `pytest.xml`, `pytest.txt` y `compilacion.json` conservan la evidencia. Python no necesita una compilación nativa para ejecutar este paquete; esa comprobación detecta errores de sintaxis, mientras que las pruebas ejercitan rutas de cálculo. Ninguna demuestra ausencia de todos los defectos posibles.

La convergencia se midió en una ejecución por resolución, con preparación incluida y métricas detalladas omitidas en malla. La escalabilidad por eventos se midió con superficie preparada, tres repeticiones y mediana, separando geometrías idénticas y ángulos distintos entre -20° y +20°. Es un benchmark controlado, no una muestra estadística de procedimientos clínicos. No deben compararse directamente sus tiempos con una ejecución completa que incluye importación, refinamiento y exportación.

### Convergencia espacial

| Arista máx. (mm) | Caras | PSD (Gy) | Área (cm²) | Tiempo (s) |
| :-- | :-- | :-- | :-- | :-- |
| 40 | 16 390 | 0,03472008 | 142,448 | 0,218 |
| 20 | 46 886 | 0,03476189 | 141,417 | 0,555 |
| 10 | 177 926 | 0,03477631 | 139,855 | 2,212 |
| 5 | 709 122 | 0,03478205 | 140,408 | 11,544 |

La PSD cambia 0,01652 % entre 10 y 5 mm, la integral dosis-área 0,34646 % y el área irradiada 0,39537 %. El punto de máximo se desplaza 2,333 mm. La estabilización del valor del máximo no garantiza igual estabilidad en su localización ni en la extensión de isodosis. Las áreas no convergen monótonamente, comportamiento compatible con clasificar caras completas por un centroide frente a bordes discontinuos.

![Convergencia propia sobre la misma anatomía XCAT. Se comparan PSD, área y tiempo; la malla más fina actúa como referencia numérica, no como verdad física.](figuras/convergencia.png){width=100%}

La malla original produce 0,0345982715 Gy, frente a 0,0347763062 Gy a 10 mm. El cambio relativo es aproximadamente 0,515 %. La primera cifra no queda invalidada en todo uso, pero demuestra que el tamaño irregular de caras puede afectar al máximo. La malla de 5 mm requiere 709 122 caras, casi cuatro veces las de 10 mm, para una mejora marginal del máximo en este caso. Por ello, 10 mm es una opción de trabajo razonable para esta demostración; no se establece una resolución clínica universal.

La prueba de refinamiento conserva el poliedro. Una esfera aproximada por pocos triángulos puede tener distancias distintas de la esfera matemática aunque se subdividan sus caras. La verificación separa ese problema mediante esferas generadas a distintas resoluciones anatómicas, y compara la intersección con la raíz analítica del rayo-esfera. La convergencia del muestreo y la convergencia de la forma son estudios diferentes.

### Tiempo por etapas y cuellos de botella

| Etapa | Caras | Mediana (s) | Repeticiones |
| :-- | :-- | :-- | :-- |
| Preparación de superficie | 177 926 | 1,358 | 3 |
| Cálculo con métricas detalladas | 177 926 | 0,285 | 3 |
| Cálculo con métricas mínimas | 177 926 | 0,257 | 3 |
| Exportación completa, malla original | 12 148 | 1,622 | 1 (tiempo único) |

En la malla refinada, el perfil de CPU sobre una superficie e índices ya usados atribuye 0,275 s de 0,326 s a proyección y 0,241 s a primera intersección. Estos tiempos son acumulativos y se solapan jerárquicamente: no deben sumarse. Las consultas al R-tree explican gran parte del tiempo de ray tracing. La multiplicación de factores representa una fracción pequeña. El coste dominante del arranque es distinto: preparar y refinar la malla puede superar al cálculo de dos eventos.

La exportación completa se midió aparte sobre la malla original de 12 148 caras y tardó 1,622 s; su HTML ocupó aproximadamente 6,77 MB. Esta medición no corresponde al mesh de 177 926 caras. Una exportación refinada iniciada durante la primera fase quedó interrumpida y no se usa como tiempo válido. Los resultados de cálculo refinado sí quedaron completados y se conservan en JSON. Esta separación evita declarar como rendimiento de cálculo el tiempo necesario para dibujar y serializar miles de caras.

### Número de eventos y efecto de la caché

| Eventos | Geometría | Sin caché (s) | Con caché (s) | Cociente |
| :-- | :-- | :-- | :-- | :-- |
| 1 | Idéntica | 0,110 | 0,112 | 0,98 |
| 4 | Idéntica | 0,385 | 0,106 | 3,64 |
| 12 | Idéntica | 1,163 | 0,126 | 9,24 |
| 24 | Idéntica | 2,354 | 0,156 | 15,07 |
| 1 | Variable | 0,376 | 0,333 | 1,13 |
| 4 | Variable | 1,018 | 0,974 | 1,04 |
| 12 | Variable | 2,592 | 2,581 | 1,00 |
| 24 | Variable | 4,879 | 4,871 | 1,00 |

Las geometrías idénticas reutilizan proyección, pero todavía requieren evaluar factores, generar mapas y acumular. Las geometrías distintas no obtienen ese ahorro. Además, un ángulo oblicuo puede intersectar muchas más cajas candidatas que una proyección frontal, aunque ilumine un número parecido de caras. Por ello, tiempo por evento no es constante respecto al ángulo.

![Escalabilidad medida con 177 926 caras, tres repeticiones y superficie preparada. Las líneas comparan caché activada y desactivada; los dos paneles son cargas de trabajo distintas.](figuras/eventos.png){width=100%}

El benchmark histórico de doce campos idénticos daba una aceleración de 8,94 veces. Es un resultado de aquella ejecución, no una promesa general. La prueba actual conserva las duraciones individuales y permite recalcular el cociente para cada carga. La igualdad exacta de mapas con y sin caché está cubierta por una prueba automatizada; la caché cambia trabajo geométrico, no el modelo físico.

### Memoria y complejidad

| Malla / arista | Caras | RSS máximo (MB) | Pares candidatos |
| :-- | :-- | :-- | :-- |
| Original | 12 148 | 94,98 | 1 087 |
| 20 mm | 46 886 | 125,37 | 21 294 |
| 10 mm | 177 926 | 226,07 | 216 970 |
| 5 mm | 709 122 | 640,18 | 2 718 483 |

Se informa RSS máximo del proceso completo, que incluye intérprete y bibliotecas y excluye exportación. No equivale a memoria de las matrices de dosis ni a asignaciones Python rastreadas por `tracemalloc`. Las primeras mediciones lanzadas desde un padre que ya había preparado la malla fina devolvieron el mismo valor elevado, 742,38 MB, para todos los casos; se descartaron como no discriminantes. La repetición desde un lanzador ligero produjo el escalado de la tabla. Ambos registros se conservan para explicar la corrección metodológica.

Con $N$ caras y $E$ eventos, almacenar todos los mapas cuesta al menos $8NE$ bytes para valores float64, más copias transitorias al apilar, geometría, árbol y resultados. La suma actual crea `np.stack(maps)`, por lo que el máximo de memoria puede superar la colección original. Una ejecución con muchos eventos necesita estudiar almacenamiento incremental o en disco. La caché de dos proyecciones aporta memoria de orden $N$; ampliar su capacidad no es gratuito.

La preselección vectorizada tiene coste aproximadamente $O(EN)$. El trazado se limita a $R$ rayos candidatos, pero depende de $C$ pares rayo-triángulo seleccionados por el índice. No puede garantizarse un coste logarítmico global solo por utilizar R-tree. En los casos de 10 y 5 mm se registraron respectivamente 216 970 y 2 718 483 pares, mostrando que la estructura de candidatos y la orientación importan además del número de caras.

### Posibilidades de optimización

La prioridad de rendimiento depende del uso. Para pocos eventos y anatomía fija, reutilizar superficie y árbol evita trabajo de preparación. Para muchos eventos repetidos, agrupar geometrías exactamente iguales o ampliar de forma acotada la caché puede ser eficaz. Para campos variables sobre mallas grandes, el perfil actual justifica evaluar un intersector con BVH o Embree antes de optimizar la multiplicación de factores. Debe compararse el mismo conjunto de rayos y controlar discrepancias en aristas y primera entrada.

La propagación LHS admite realizaciones independientes y es candidata a multiprocessing. Sin embargo, copiar una malla grande por proceso puede multiplicar RSS; hay que ensayar memoria compartida, número de trabajadores y orden reproducible de semillas. La GPU puede acelerar operaciones geométricas o transporte, pero añade transferencias, precisión y dependencias de hardware. Los resultados GPU publicados por Rana y Takata corresponden a sistemas y cargas distintas y no predicen el rendimiento de este programa [6, 7].

El visor también necesita optimización propia: clipping en JavaScript y actualización de todas las caras pueden limitar la interacción aun cuando el cálculo termine rápido. Una malla gráfica reducida solo sería aceptable manteniendo correspondencia con la malla de cálculo y evitando suavizar máximos. El tiempo de primer cálculo, el tiempo de exportación y la fluidez del navegador deben medirse como variables distintas.

## Verificación física del modelo implementado

### Casos analíticos y consistencia

La comprobación punto a punto del inverso del cuadrado en el plano obtuvo error relativo máximo $2,22\times10^{-16}$ frente a la fórmula independiente. Al duplicar la distancia axial de 600 a 1200 mm, el cociente en el centroide próximo al eje fue 0,25000651036, idéntico a la expresión con su desplazamiento transversal. La pequeña diferencia respecto a 0,25 se debe a que ese centroide no está exactamente sobre el eje, no a una violación de la ley.

Para ángulos -40°, -20°, 0°, +20° y +40°, con campo suficientemente amplio, la comparación de dosis en los puntos iluminados dio errores relativos máximos menores que $7\times10^{-16}$. Esto verifica distancia y evaluación algebraica bajo los supuestos del modelo. No valida la retrodispersión angular, la respuesta real de mesa ni la convención angular de un equipo clínico.

La suite incluye máscara rectangular analítica, transmisión únicamente si se cruza la mesa, factores por modo, interpolación y rechazo de extrapolación, acumulación y ausencia de intersección, barridos con conservación de kerma, coordenadas transformadas, caras ocultas, orientación global y topología de refinamiento. Las pruebas de esfera y cilindro aportan referencias geométricas distintas de los ejemplos antropomórficos. La prueba de autosombra evita atribuir entrada a una cara visible por normal pero oculta por una extremidad.

Se rechazaron mediante pruebas entradas con unidades no reconocidas, kerma ausente, calidad de haz insuficiente, UID duplicado, fabricante distinto, malla abierta y números no finitos. Esta cobertura es relevante para fallo controlado, aunque no constituye una auditoría completa de cada campo JSON ni de todas las jerarquías DICOM posibles.

### Dependencias presentes y ausentes

Manteniendo $K_{a,r}$ y geometría, cambiar el kVp de 60 a 80 o 120 produjo exactamente 0,0174993924 Gy en el caso plano de un evento. Cambiar únicamente `thickness_mm` entre 100, 220 y 400 mm produjo la misma cifra. Es el comportamiento esperado con factores constantes, y demuestra una limitación del ejemplo: no existe una respuesta espectral o de espesor emergente del motor.

Las tablas pueden introducir esas dependencias y las pruebas sintéticas verifican su interpolación. No obstante, una tabla artificial de BSF lineal con ángulo solo demuestra evaluación correcta de una función. No prueba que la función represente retrodispersión física. Para verificar atenuación con espesor atravesado sería necesario un modelo que la implementase o una tabla calibrada con alcance explícito; el motor actual no permite afirmar esa validación.

El nombre de una prueba histórica incluye «energy conservation», pero sus aserciones comprueban producto de factores y máscara, no un balance energético completo. La conservación demostrada en el refinamiento corresponde a área y volumen geométricos; la del barrido corresponde al kerma integrado repartido. Ninguna equivale a conservación de energía en transporte radiativo.

### Resultados de las demostraciones

| Caso sintético | Eventos | PSD (Gy) | Área irradiada (cm²) |
| :-- | :-- | :-- | :-- |
| XCAT original | 2 | 0,03459827 | 143,625 |
| XCAT, arista 10 mm | 2 | 0,03477631 | 139,855 |
| Secuencia intervencionista, 10 mm | 24 | 1,05631146 | 415,049 |

En el caso XCAT refinado, la media local aproximada es 0,0343266488 Gy sobre 89,8518 mm², menor que el máximo puntual de 0,0347763062 Gy. Esa diferencia es coherente con promediar una distribución no uniforme. El área irradiada es 139,855 cm² y aparece una región conectada por encima del 90 % del máximo, de 104,789 cm². La posición del máximo es aproximadamente (-51,257; -30,033; -96,623) mm en coordenadas locales.

![Mapa superficial de la demostración XCAT original regenerado durante la auditoría. Se muestra una geometría de referencia, no una anatomía individual. La resolución numérica de esta imagen es 12 148 caras.](evidencias/xcat_original/mapa_3d.png){width=75%}

El ejemplo intervencionista separado contiene 24 eventos sintéticos, con 18 eventos de fluoroscopia y seis de adquisición. Su kerma acumulado es 1,26 Gy y su PSD calculada aproximadamente 1,05631 Gy. Una PSD menor que el kerma acumulado es posible por distancia, factores y distribución espacial; no demuestra que sea una dosis realista para todo procedimiento de esa clase. Las cifras se han generado para ejercitar el programa, no proceden de un paciente ni de una exposición medida.

La semejanza del ejemplo con una secuencia intervencionista es estructural: varias proyecciones, modos, campos y contribuciones desiguales. Los tiempos de fluoroscopia del guion de ejemplo no constituyen un determinante independiente de dosis en este motor, que usa kerma integrado por evento. No puede inferirse exactitud clínica comparando únicamente el orden de magnitud con procedimientos publicados.

### Propagación paramétrica

| Estadístico de PSD | Resultado (Gy) |
| :-- | :-- |
| Media, 64 muestras LHS | 0,03478935 |
| Desviación estándar muestral | 0,00208814 |
| Cuantil 2,5 % | 0,03143112 |
| Cuantil 50 % | 0,03471816 |
| Cuantil 97,5 % | 0,03952601 |

La propagación de 64 muestras usa desplazamientos normales del paciente con desviación típica de 3 mm por eje, sesgos angulares de 0,5°, escala lognormal de kerma con desviación 0,05 y de BSF con desviación 0,03, ambas con media uno. La semilla es 20260928. Son distribuciones hipotéticas, no estimadas de medidas. Los sesgos son comunes a todos los eventos de cada realización.

La ejecución histórica de 256 muestras con otra semilla dio media 0,0347825785 Gy y desviación 0,0021075133 Gy, próximas a las de 64 muestras. Esa proximidad es una comprobación parcial de estabilidad, no una prueba de convergencia de colas: al cambiar simultáneamente tamaño y semilla no se aísla el efecto de cada uno. Un estudio formal debería repetir varias semillas por tamaño y analizar cuantiles, localización de picos y probabilidades por cara.

Las muestras fuera del dominio o físicamente inválidas detienen la propagación; no se eliminan silenciosamente. `source_x/y/z_mm` desplaza el sistema fuente-isocentro-detector, no un foco móvil respecto a un detector fijo. El mapa de frecuencia de excedencia describe únicamente las distribuciones introducidas y no es una probabilidad de lesión. Excluye error anatómico, discretización, dispersión no modelada, correlaciones e incertidumbre de interpolación.

## Validación externa y comparación bibliográfica

La evidencia propia no contiene medidas de cámara, película, TLD, OSLD ni resultados de transporte independientes. Los CSV de `examples` son plantillas. Las pruebas del módulo de comparación usan referencias artificiales o valores derivados del propio mapa; verifican fórmulas y registro, pero no el fenómeno físico. El estado `validated=False` se mantiene deliberadamente.

Gañán Mora et al. realizaron 27 irradiaciones en dos equipos, con referencia de cámara en un fantoma de agua; informaron incertidumbre experimental del 13 % con $k=2$ y diferencias de -9 % a +8 % para SkinDose 2D [2]. Esta campaña ofrece un modelo reproducible de evaluación por parámetros. No se han procesado sus RDSR suplementarios en el proyecto ni se puede asignar ese intervalo al calculador.

Krajinović et al. validaron SkinCare con película y procedimientos simplificados, y discutieron la heterogeneidad de RDSR y la localización del paciente [4]. Malchair et al. observaron que resultados aceptables en fantoma podían deteriorarse en pacientes [5]. Dabin et al. encontraron que ningún producto estudiado cumplía ±40 % en todas las combinaciones de sistemas compatibles y configuraciones; ese dato motiva ensayos laterales y multiequipo, no adopción automática de ±40 % como tolerancia [27].

Vergara et al. compararon cuatro fantomas computacionales con TLD en dos procedimientos sobre fantoma antropomórfico. El modelo humano dio PSD de 3,163 y 3,728 Gy, frente a 3,161 y 3,713 Gy medidos [13]. Son resultados de sus configuraciones y calibraciones; no prueban que cualquier malla humana produzca ese acuerdo. La lección trasladable es combinar representación anatómica, factores medidos y referencia dosimétrica definida.

Takata et al. informaron diferencias de máximo de hasta 6,1 % respecto a PHITS y aceleración GPU en su sistema de transporte [7]. Han et al. informaron diferencias de PSD dentro del 3 % respecto a simulaciones de alta estadística para los casos satisfactorios de su metodología de reducción de ruido [12]. No existe aquí una reproducción de esos experimentos, y las distintas definiciones de piel, discretización, normalización y estadística impiden usar sus porcentajes como validación indirecta.

La comparación comercial `compare-radimetrics` comprueba medio y devuelve diferencias, pero el $K_{a,r}$ de nuestra rama principal sigue procediendo del RDSR. Dos aplicaciones pueden coincidir por compartir datos y supuestos. La rama `tube_output` solo es independiente si rendimiento y mAs son fiables y caracterizados; en la demo, tanto el rendimiento como el evento son construidos para coincidir.

## Protocolo de validación pendiente

Una primera campaña debe fijar la magnitud: kerma libre en aire, kerma incidente o dosis de entrada en un medio y profundidad especificados. El área o volumen sensible del detector debe relacionarse con el muestreo superficial. Hay que registrar fantoma, fuente, mesa y campo en un marco común y documentar incertidumbres de posición. Sin ese contrato, una diferencia de dosis puede ser una diferencia de definición.

| Fase propuesta | Experimento | Criterio que permite investigar |
|:--|:--|:--|
| Indicador y rendimiento | Cámara trazable, distancia y mAs variados | Linealidad, calibración y rama independiente |
| Mesa | Mismo haz con y sin mesa/colchoneta | Transmisión por calidad, campo y ángulo |
| Plano de agua | Distancias, campos y ángulos por separado | Ecuación dosimétrica y geometría básica |
| Perfil y borde | Película calibrada y registrada | Talón, penumbra y dosis fuera del campo |
| Curvatura | Cilindro y fantoma antropomórfico | Superficie, proyecciones laterales y oclusión |
| Superposición | Campos coincidentes, separados y traslaciones | Acumulación, pico y área |
| Portabilidad | RDSR de cada modelo/versión | Códigos, jerarquía, signos y distancias |
| Repetibilidad | Reposicionamiento y operadores distintos | Distribuciones y correlaciones de entrada |

Cada fase debe separar datos de calibración de datos de validación. Las tolerancias se fijarán antes del análisis para un uso previsto concreto y con presupuesto de incertidumbre. Puede evaluarse $\Delta=D_{calc}-D_{ref}$ y una diferencia normalizada por la incertidumbre combinada; el módulo actual solo divide por la incertidumbre de la medida, de modo que no debe interpretarse como una puntuación estadística completa del acuerdo.

Para mapas, `validate-map` asocia cada punto medido al centroide más próximo mediante KD-tree, exige distancia máxima y evita asociaciones duplicadas. Calcula RMSE, separación de picos entre posiciones muestreadas y Dice binario en esos puntos. No registra imágenes, no integra la huella del detector y no calcula gamma ni Dice de toda la superficie. Es una base para comparar, no un protocolo metrológico resuelto.

Una comparación Monte Carlo debe modelar espectro, mesa, fantoma, materiales y normalización, así como scorer, espesor de piel, cortes y estadística. El máximo de muchas celdas ruidosas puede sesgarse al alza; debe estudiarse su estabilidad. Solo después de comparar geometrías simples sería razonable abordar un fantoma completo y secuencias clínicas. No se propone usar un TPS de megavoltaje como referencia de un arco kV sin demostrar equivalencia del modelo de haz.

## Discrepancias documentales y revisión crítica

Se identifican varias diferencias entre intención, descripción previa e implementación. El indicador `resolution_mm` no refina mallas; la versión actual dispone de un parámetro específico. La escala cromática del HTML es fija, aunque un informe previo menciona escala automática. «Cortes» significa ocultación de superficie, no planos de dosis en órganos. El módulo de incertidumbre no es transporte Monte Carlo y la etiqueta `skin` no crea un material de piel.

La documentación de cambios menciona prueba de escala lognormal, pero la comparación analítica explícita de la suite usa una escala uniforme; la distribución lognormal se ejercita en el estudio, sin equivaler a esa prueba analítica. El nombre de test sobre conservación energética es más amplio que lo comprobado. La revisión corrige la interpretación, manteniendo los resultados verificables.

Hay riesgos de uso no resueltos: factores escalares sin rango declarado, normalización de perfil no comprobada, incertidumbre parcial cero cuando no hay componentes conocidas, autointersecciones no verificadas y ausencia de validación DICOM completa. El software registra varios de ellos, pero un consumidor que lea solo PSD o un número de incertidumbre puede ignorar el contexto. La robustez de presentación debe evolucionar junto al modelo.

# Limitaciones

## Límites físicos

La formulación describe entrada primaria con correcciones multiplicativas. Excluye transporte en tejido, heterogeneidad, dosis de salida y dispersión lateral fuera de la huella. El BSF local no representa explícitamente el acoplamiento espacial entre regiones. La calidad del haz interviene únicamente a través de datos tabulados; las constantes sintéticas no cambian con kVp ni espesor. No existe predicción de profundidad ni de respuesta biológica.

La mesa es un plano acotado con transmisión efectiva. Un rayo que atraviese un raíl o un accesorio no se distingue por material. Tampoco se modelan colimadores asimétricos, filtros de cuña ni penumbra de foco finito. Estas limitaciones pueden ser especialmente importantes en oblicuidades, campos pequeños y regiones cercanas a bordes. Son límites del modelo, no fallos que desaparezcan con más pruebas unitarias.

## Límites geométricos y numéricos

XCAT es una anatomía de referencia colocada de forma aproximada. No hay segmentación individual, registro automático, deformación ni aplanamiento por apoyo. La distancia fuente-piel y la coincidencia de campos pueden ser erróneas aun con cálculo geométrico perfecto. Una malla cerrada y de volumen positivo no garantiza ausencia de autointersecciones.

La dosis por centroide aproxima toda una cara. Las métricas de área y el máximo son sensibles a bordes y a resolución desigual. La media local de 1 cm² tiene área variable y distancia de grafo aproximada. Las tolerancias de primera entrada necesitan más casos degenerados. El cilindro abierto solo es admisible en el dominio geométrico controlado por el programa.

## Límites de implementación y evidencia

El lector RDSR cubre un subconjunto y un filtro Cu; no es compatible por demostración con todos los fabricantes. Como avance práctico derivado de las notas de revisión, se ha añadido un comando de búsqueda automática de RDSR en carpetas DICOM. Esta herramienta no amplía la física del modelo, pero reduce un problema operativo real: distinguir el informe estructurado de dosis de las imágenes XA/RF de un procedimiento. El comando solo lee metadatos, resume los objetos no candidatos y devuelve hashes de UIDs en lugar de identificadores directos. No hay corpus clínico independiente ni comparación con la consola de un equipo real. El cálculo biplano depende de caracterizar cada plano y de registrar el mismo paciente. Los barridos uniformes no sustituyen logs temporales.

La memoria crece con caras y eventos; la exportación HTML duplica representaciones y puede ser pesada. No hay servicio en tiempo real ni pruebas de rendimiento en Windows. La lógica del visor se prueba con un sustituto de DOM/Plotly, no con renderizado WebGL real en esta auditoría. Las capturas estáticas verificadas no garantizan fluidez ni ausencia de problemas gráficos en todos los navegadores.

La exactitud física y clínica sigue sin cuantificar. No hay calibración real, medidas independientes ni validación de dosis a una profundidad definida. La incertidumbre paramétrica depende de hipótesis incompletas. Ningún porcentaje de convergencia, prueba superada o acuerdo con valores sintéticos puede convertirse en un error clínico global.

# Conclusiones

Se ha desarrollado un prototipo modular capaz de leer eventos RDSR, reconstruir geometría, seleccionar primera entrada sobre mallas triangulares y acumular dosis superficial mediante factores configurables. La salida conserva trazabilidad y permite explorar mapas, áreas y contribuciones individuales. Además de conservar los mapas por evento, la versión actual ordena las contribuciones de cada evento en el punto de PSD acumulada, lo que permite interpretar si el máximo procede de un campo dominante o de la superposición de varios campos. El refinamiento real de la malla y la separación entre cálculo y visualización resuelven necesidades computacionales concretas del proyecto.

La nueva evaluación confirma 54 pruebas superadas, consistencia con soluciones geométricas y estabilidad del máximo en la demostración XCAT. El refinamiento de 10 a 5 mm cambia la PSD un 0,0165 %, con mayor variación en áreas y mayor coste. El rendimiento depende de preparación, ángulo, número de candidatos y repetición exacta de geometrías; no de un único tiempo universal por evento.

El objetivo de construir una plataforma de investigación verificable está cumplido. El de estimar dosis clínica con precisión conocida no lo está: faltan datos de equipo, registro anatómico y validación externa. La aportación actual es una implementación auditable del formalismo aproximado y una base de pruebas reproducibles, acompañada de una delimitación explícita de lo que todavía no demuestra.

# Líneas futuras de desarrollo

## Caracterización dosimétrica y control del dominio

La prioridad es sustituir factores sintéticos por medidas trazables para equipo, modo, calidad, campo y ángulo. Deben definirse con precisión las magnitudes de CF, transmisión, BSF y conversión para evitar doble conteo. La configuración debería exigir rangos válidos incluso para constantes, registrar incertidumbre y fecha de calibración, y rechazar perfiles incompatibles con su normalización. Esta propuesta responde directamente a la invariancia artificial observada con kVp y espesor, y al problema de caracterización señalado por la literatura [2, 5, 27].

La salida de incertidumbre parcial debería expresar «no cuantificada» cuando no existe ninguna componente conocida, en lugar de ofrecer un cero aislado. La interfaz debe mostrar junto al máximo etapa, medio, estado de calibración y dominio. Es una mejora de interpretación verificable con pruebas de informes y casos incompletos, sin alterar la ecuación de dosis.

## Anatomía, registro y datos clínicos

La siguiente mejora geométrica es registrar una superficie medida o segmentada en el sistema de mesa. Deben evaluarse distancias de superficie, landmarks y sensibilidad de PSD al error de posición. El escalado antropométrico de una referencia es una opción intermedia, pero no garantiza reproducir apoyo y concavidades. Johnson et al. mostraron que el beneficio del ajuste anatómico depende de la proyección; Omar et al. abordaron explícitamente la reconstrucción de posición desde datos radiológicos [3, 8].

El lector necesita un corpus de RDSR reales anonimizados, pruebas de jerarquía y perfiles por fabricante/modelo/versión. Enhanced RDSR puede facilitar marcos y transformaciones, pero requiere otro soporte explícito [15]. La integración PACS o DICOMweb solo debería realizarse después de resolver identidad, unidades, dominio y trazabilidad. La ampliación a movimiento demanda logs temporales o seguimiento; no debe inferirse movimiento no registrado.

## Resolución y cálculo eficiente

Se propone refinamiento guiado por error: comparar centroides con submuestras en caras cercanas a bordes y máximos, refinar donde PSD o área no alcancen una tolerancia y conservar topología. Su aceptación exige casos de campo estrecho, oblicuidad, autosombra y mapas con varios máximos. El criterio debería incluir área e integral, no solo PSD. La diferencia entre refinamientos no debe llamarse intervalo probabilístico.

El perfil medido justifica evaluar BVH/Embree para primera entrada, acumulación por bloques para reducir memoria y caché configurable para secuencias repetidas. La propagación admite paralelización por realizaciones con control de RSS. La GPU se considerará después de un benchmark equivalente que incluya transferencias y precisión, siguiendo la lección de Takata sobre conservar datos en memoria del dispositivo [7]. No se presupone que añadir GPU sea la mejora dominante para dos eventos.

## Física adicional y validación

Para abordar dispersión fuera del campo o salida, se propone primero comparar con transporte independiente y medir perfiles. Si el sesgo es sistemático y parametrizable, podrían construirse kernels o tablas validadas. Para dosis en profundidad o heterogeneidades sería necesario un modelo volumétrico. No se recomienda añadir factores heurísticos sin separar datos de ajuste y validación.

La validación por componentes de García Balcaza y la preocupación por posición y automatización de Fernández Bosman son trasladables al protocolo, aunque sus sistemas de transporte y dosimetría ocupacional son diferentes [9, 10]. Las líneas futuras de SkinCare subrayan ampliar validación clínica y montajes [4]. Han et al. proponen integración en NCIRF y biblioteca anatómica ampliada [12]; aquí ello inspira una biblioteca versionada de superficies y casos, no una declaración de integración ya realizada.

## Usabilidad y criterios de avance

El visor debe someterse a pruebas WebGL en Windows, con mallas grandes, varias resoluciones de pantalla y medición de latencia al mover cortes. Es conveniente un modo gráfico ligero con correspondencia inequívoca hacia los valores originales. Un selector de caso y configuración reduciría confusión entre demo, investigación y ejemplo, pero debe conservar comandos reproducibles y mensajes de error completos.

El avance se considerará demostrado cuando existan calibraciones documentadas, reproducción de casos analíticos y externos, tolerancias previamente definidas, presupuesto de incertidumbre adecuado al uso, corpus multiequipo y controles de regresión tras cada cambio. Hasta entonces, ampliar la complejidad física o mejorar la apariencia 3D no sustituye la evidencia que falta.

## Comentarios de revisión incorporados

Comentario añadido sobre RDSR reales: una carpeta de procedimiento no debe introducirse completa en el motor de cálculo sin identificar antes el objeto `X-Ray Radiation Dose SR`. Por ello se ha implementado `find-rdsr`, que localiza candidatos RDSR y señala de forma explícita los ficheros que son imágenes angiográficas u otros DICOM no utilizables como entrada principal. Esta mejora ya está incorporada al calculador y validada con pruebas sintéticas.

Comentario añadido sobre compatibilidad de RDSR reales: se ha incorporado `compatibility-report`, que genera un informe previo para RDSR o carpetas DICOM. Este informe comprueba lectura del RDSR, campos críticos por evento, disponibilidad de calidad de haz para la etapa solicitada y coincidencia de fabricante/modelo con la configuración. No equivale a validación física, pero responde a la necesidad de aceptar RDSR reales anonimizados de forma controlada.

Comentario añadido sobre análisis por evento: el cálculo acumulado puede ocultar qué evento origina el máximo. La nueva clave `global_psd_event_contributions` cuantifica la contribución de cada evento en el punto de PSD acumulada. Esta función no modifica la dosis, pero mejora la interpretación del resultado y responde al objetivo de aprovechar mejor la información evento a evento del RDSR.

Comentario añadido sobre trazabilidad de factores: las notas insistían en diferenciar factores sintéticos, medidos, tabulados o calculados, y en conservar su rango de validez e incertidumbre. Por ello se ha añadido `factor_metadata` a la salida. Este bloque documenta para cada corrección su identificador, procedencia, representación, incertidumbre relativa y dominio tabulado cuando existe. Así, las futuras tablas de \(BSF\), transmisión de mesa o \(f_{\text{skin}}\) podrán incorporarse sin confundirse con valores de demostración.

Comentario añadido sobre trazabilidad de las anotaciones: se ha creado `NOTAS_IMPLEMENTACION.md`, que agrupa todas las líneas iniciadas por `% Nota:` y declara si cada propuesta está implementada, preparada para datos calibrados o pendiente por requerir medidas, tablas externas o simulaciones.

Comentario añadido sobre mejoras no implementadas todavía: la calibración mesa--equipo, las tablas de BSF/transmisión/factor aire--piel, el uso de superficies segmentadas de paciente y la incorporación de modelos antropométricos personalizados se mantienen como líneas prioritarias, pero no deben codificarse como constantes o geometrías ficticias sin medidas, procedencia y rango de validez. En esta versión se ha preparado el software para aceptarlas mediante configuración e interpolación acotada, pero la generación de datos calibrados queda pendiente.

# Bibliografía {.unnumbered}

[1] Andersson J, Bednarek DR, Bolch W, et al. Estimation of patient skin dose in fluoroscopy: summary of a joint report by AAPM TG357 and EFOMP. *Medical Physics*. 2021;48(7):e671–e696. DOI: [10.1002/mp.14910](https://doi.org/10.1002/mp.14910). Ficha oficial AAPM y metadatos/resumen consultados; no lectura completa de todos los apéndices.

[2] Gañán Mora A, Sánchez Casanueva RM, Fernández Soto JM. Metodología para la validación de programas de cálculo de dosis en piel. *Revista de Física Médica*. 2024;25(2):25–34. DOI: [10.37004/sefm/2024.25.2.002](https://doi.org/10.37004/sefm/2024.25.2.002). PDF aportado, texto completo disponible.

[3] Johnson PB, Borrego D, Balter S, Johnson K, Siragusa D, Bolch WE. Skin dose mapping for fluoroscopically guided interventions. *Medical Physics*. 2011;38(10):5490–5499. DOI: [10.1118/1.3633935](https://doi.org/10.1118/1.3633935). PDF aportado.

[4] Krajinović M, Kržanović N, Ciraj-Bjelac O. Vendor-independent skin dose mapping application for interventional radiology and cardiology. *Journal of Applied Clinical Medical Physics*. 2021;22(2):145–157. DOI: [10.1002/acm2.13167](https://doi.org/10.1002/acm2.13167). PDF aportado.

[5] Malchair F, Dabin J, Deleu M, Sans Merce M, Ciraj Bjelac O, Gallagher A, Maccia C. Review of skin dose calculation software in interventional cardiology. *Physica Medica*. 2020;80:75–83. DOI: [10.1016/j.ejmp.2020.09.023](https://doi.org/10.1016/j.ejmp.2020.09.023). PDF aportado.

[6] Rana V, Rudin S, Bednarek DR. Use of a graphics processing unit (GPU) to facilitate real-time 3D graphic presentation of the patient skin-dose distribution. *Proceedings of SPIE*. 2012;8313:831343. DOI: [10.1117/12.911344](https://doi.org/10.1117/12.911344). Manuscrito aportado.

[7] Takata T, Kotoku J, Maejima H, et al. Fast skin dose estimation system for interventional radiology. *Journal of Radiation Research*. 2018;59(2):233–239. DOI: [10.1093/jrr/rrx062](https://doi.org/10.1093/jrr/rrx062). PDF aportado; publicación anticipada en 2017.

[8] Omar A, Bujila R, Fransson A, Andreo P, Poludniowski G. A framework for organ dose estimation in x-ray angiography and interventional radiology based on dose-related data in DICOM structured reports. *Physics in Medicine & Biology*. 2016;61(8):3063–3083. DOI: [10.1088/0031-9155/61/8/3063](https://doi.org/10.1088/0031-9155/61/8/3063). Manuscrito aceptado aportado; no confundir fecha de depósito con año de publicación.

[9] Fernández Bosman D. *PyMCGPU-IR: a new tool for patient dose monitoring in interventional radiology procedures*. Tesis doctoral, Universitat Politècnica de Catalunya, 2023. PDF `TDFB1de1.pdf`; especialmente capítulo 12, pp. 114–116.

[10] García Balcaza V. *PyMCGPU-IR, a computational Monte Carlo application for occupational dosimetry in interventional radiology procedures*. Tesis doctoral, Universitat Politècnica de Catalunya, 2024. PDF `TVGB1de1.pdf`; discusión y conclusiones, p. 162.

[11] Omar A. *Radiation dose and X-ray beam modelling in diagnostic and interventional radiology using Monte Carlo methods*. Tesis doctoral, Karolinska Institutet, Estocolmo, 2020. ISBN 978-91-8016-009-4. PDF aportado.

[12] Han H, Streitmatter SW, Kitahara CM, Lee C. Fast and accurate peak skin dose estimation method for interventional fluoroscopy patients. *Medical Physics*. 2025;52(4):2551–2559. DOI: [10.1002/mp.17667](https://doi.org/10.1002/mp.17667). PDF editorial consultado.

[13] Vergara D, Makkia RS, Li Z, Hoerner M, Fadhel MN. Experimental validation of a comprehensive fluoroscopy peak skin dose model using four different computational phantoms. *Medical Physics*. 2025;52(6):3635–3648. DOI: [10.1002/mp.17737](https://doi.org/10.1002/mp.17737). Resumen primario PubMed, PMID 40091498; no lectura completa atribuida.

[14] JCGM. *Evaluation of measurement data: Supplement 1 to the Guide to the expression of uncertainty in measurement. Propagation of distributions using a Monte Carlo method*. JCGM 101:2008. BIPM. DOI: [10.59161/JCGM101-2008](https://doi.org/10.59161/JCGM101-2008).

[15] DICOM Standards Committee. *DICOM PS3.16, TID 10003 Irradiation Event X-Ray Data*; *PS3.17, UUUU.4.5.4 Geometry and Radiation Output Encoding in Enhanced RDSR*. Edición web consultada el 28-09-2026. [TID 10003](https://dicom.nema.org/medical/dicom/current/output/chtml/part16/sect_TID_10003.html); [Enhanced RDSR](https://dicom.nema.org/medical/dicom/current/output/chtml/part17/sect_uuuu.4.5.4.html).

[16] Auer B. *Mesh-based Human Phantom for Simulation*. Repositorio científico, revisión del activo `1f20b67fd7fa700d27a3dc88ee9ece436658eee0`. [Repositorio](https://github.com/BenAuer2021/Mesh-based-Human-Phantom-for-Simulation). Licencia y atribución conservadas en el proyecto.

[17] Trimesh contributors. *trimesh.ray.ray_triangle* y documentación de geometría de Trimesh. [Documentación primaria](https://trimesh.org/trimesh.ray.ray_triangle.html). Versión ejecutada 5.1.0; consulta 28-09-2026.

[18] SciPy contributors. *RegularGridInterpolator* y *scipy.stats.qmc.LatinHypercube*. [Interpolación](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html); [LHS](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.LatinHypercube.html). Documentación web consultada 28-09-2026; versión ejecutada 1.17.0, distinta de la documentación web actual 1.18.0.

[19] IAEA. *Dosimetry in Diagnostic Radiology: An International Code of Practice*. Technical Reports Series 457. Viena: IAEA; 2007. ISBN 92-0-115406-2. [Documento oficial](https://www-pub.iaea.org/MTCD/Publications/PDF/TRS457_web.pdf).

[20] ICRU. *Patient Dosimetry for X Rays Used in Medical Imaging*. ICRU Report 74. 2005. [Ficha oficial](https://www.icru.org/report/patient-dosimetry-for-x-rays-used-in-medical-imaging-report-74/). Se consultó la ficha, no se atribuye lectura íntegra del informe.

[21] ICRP. *Adult Mesh-type Reference Computational Phantoms*. ICRP Publication 145. 2020. [Ficha oficial](https://www.icrp.org/publication.asp?id=ICRP+Publication+145).

[22] Hubbell JH, Seltzer SM. *Tables of X-Ray Mass Attenuation Coefficients and Mass Energy-Absorption Coefficients*. NIST, base de referencia y sección 3 sobre coeficientes de transferencia y absorción. [Tablas](https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients); [fundamento](https://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html). Consulta 28-09-2026.

[23] ICRP. *Radiological Protection in Cardiology*. ICRP Publication 120. 2013. [Ficha oficial](https://www.icrp.org/publication.asp?id=ICRP+Publication+120).

[24] ICRP. *ICRP Statement on Tissue Reactions / Early and Late Effects of Radiation in Normal Tissues and Organs: Threshold Doses for Tissue Reactions in a Radiation Protection Context*. ICRP Publication 118. 2012. [Ficha oficial](https://www.icrp.org/publication.asp?id=ICRP+Publication+118).

[25] Benmakhlouf H, Bouchard H, Fransson A, Andreo P. Backscatter factors and mass energy-absorption coefficient ratios for diagnostic radiology dosimetry. *Physics in Medicine & Biology*. 2011;56(22):7179–7204. DOI: [10.1088/0031-9155/56/22/012](https://doi.org/10.1088/0031-9155/56/22/012). Resumen primario, PMID 22024474.

[26] Jones AK, Pasciak AS. Calculating the peak skin dose resulting from fluoroscopically guided interventions. Part I: Methods. *Journal of Applied Clinical Medical Physics*. 2011;12(4):3670. DOI: [10.1120/jacmp.v12i4.3670](https://doi.org/10.1120/jacmp.v12i4.3670). Metadatos y referencia del método; no reproducción de su validación.

[27] Dabin J, Blidéanu V, Ciraj Bjelac O, Deleu M, De Monte F, Feghali JA, Gallagher A, Knežević Ž, Maccia C, Malchair F, Sans Merce M, Simantirakis G. Accuracy of skin dose mapping in interventional cardiology: Comparison of 10 software products following a common protocol. *Physica Medica*. 2021;82:279–294. DOI: [10.1016/j.ejmp.2021.02.016](https://doi.org/10.1016/j.ejmp.2021.02.016). Resumen PubMed y registro institucional SCK CEN consultados.

[28] Wunderle KA, Balter S. Patient radiation dosimetry for fluoroscopically guided interventions: a review of the history and limitations. *Medical Physics*. 2026;53(1):e70211. DOI: [10.1002/mp.70211](https://doi.org/10.1002/mp.70211). Resumen editorial y PubMed, PMID 41490982; acceso directo PMC intermitente.

[29] IAEA. *Quality Assurance and Optimization for Fluoroscopically Guided Interventional Procedures*. Human Health Series 48. Viena: IAEA; 2025. DOI: [10.61092/iaea.f7bd-coh3](https://doi.org/10.61092/iaea.f7bd-coh3). Ficha oficial consultada; no se declara conformidad del programa.

[30] Grupo de Trabajo de Dosimetría en Radiología Intervencionista. *Objetivos del grupo de trabajo*. Documento aportado `Grupo-de-Trabajo-sobre-Dosimetria-en-RI.pdf`, sin fecha editorial identificada. Se utiliza como contexto de objetivos profesionales, no como validación ni norma.

# Anexos de reproducibilidad {.unnumbered}

## Alcance de la auditoría {.unnumbered}

El paquete de evidencias contiene inventario y huellas del proyecto, entorno de ejecución, resultados de pruebas, convergencia, tiempos por etapa, escalabilidad, memoria corregida y descartada, perfil de CPU, comprobaciones físicas, resultado XCAT y ejemplo intervencionista. `audit_study.py` genera el estudio sin modificar el motor. Los resultados de 256 muestras se identifican como históricos y se conservan aparte, evitando presentarlos como una segunda ejecución de esta memoria.

Las rutas de código de la tabla de trazabilidad son relativas al paquete `skindose_verifier`. El inventario permite identificar el archivo exacto aunque cambie el número de línea. Los hashes de código se registran en `entorno.json`; el hash del STL principal se muestra en el apartado anatómico. Los binarios gráficos y mapas previos se inventariaron y contrastaron mediante las salidas numéricas pertinentes, sin asumir que cada copia generada constituye un desarrollo distinto.

## Reproducción del estudio {.unnumbered}

Desde una carpeta que contenga `skin_dose_project`, `memoria` y el ejemplo en `deliverables/Ejemplo_Intervencion_Sintetica`, instalar el paquete y ejecutar el script del anexo. La variable `PYTHONPATH` específica de la ejecución original se sustituye por una instalación normal del paquete:

```bash
python -m pip install -e "skin_dose_project[test]"
python memoria/audit_study.py
```

Para reproducir versiones en el mismo Python y plataforma, consultar `requirements-research-lock.txt`. La medición de memoria válida se lanza desde un proceso ligero; el bucle Linux utilizado es:

```bash
for edge in 0 20 10 5; do
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python \
    skin_dose_project/scripts/measure_memory.py \
    --config skin_dose_project/results/research_xcat/equipo_sintetico.json \
    --rdsr skin_dose_project/results/research_xcat/rdsr_sintetico.dcm \
    --edge "$edge"
done
```

El script usa `/proc/cpuinfo` para describir el entorno Linux. Para Windows debe adaptarse esa captura y utilizarse una medida de memoria específica de la plataforma; el script original de RSS informa ausencia cuando no dispone de `resource`. Los BAT no se ejecutaron en Windows en esta auditoría.

## Estrategia y límites bibliográficos {.unnumbered}

La búsqueda se organizó por formalismo dosimétrico, validación, anatomía, aceleración y estándares. Se incluyeron estudios directamente relacionados con fluoroscopia kV, documentos de organismos y documentación primaria de las bibliotecas utilizadas. Se descartaron como evidencia del modelo resultados de radioterapia MV, contaminación beta y fuentes comerciales sin validación verificable. Se distinguió fecha editorial de fecha de actualización web; un copyright de 2026 no se interpretó como publicación nueva.

Los once documentos aportados se incorporaron como corpus local. En tesis extensas se revisaron los capítulos pertinentes de métodos, limitaciones y futuro; no se afirma haber reproducido todos sus experimentos. El trabajo sobre dosis de órganos [8] se utiliza para arquitectura y registro, y la tesis ocupacional [10] para validación y captura geométrica, no para atribuir capacidades orgánicas u ocupacionales al calculador.

No se obtuvo el suplemento RDSR de [2], no se accedió íntegramente a todos los informes de pago y no se verificaron todas las cláusulas IEC. Por ello, no se declara cumplimiento normativo ni se inventa una edición o tolerancia reglamentaria. Los enlaces y DOI permiten ampliar la revisión. Los PDF de terceros no se redistribuyen en el paquete de evidencias.

## Comprobación final de coherencia {.unnumbered}

Se contrastaron ecuación central, orden de factores, marcos, máscara, medio, unidades, tamaños de malla y métricas con sus implementaciones. Se distinguieron resultados nuevos de históricos y de literatura. Se revisaron afirmaciones sobre energía, espesor, GPU, resolución, incertidumbre y visor para evitar describir funciones no presentes. Los defectos o carencias encontrados se documentaron sin modificar el motor durante la evaluación, de modo que los números corresponden a una única versión auditada.

La conclusión de validez se limita a lo observado: el prototipo resuelve consistentemente los casos verificados de su formalismo aproximado. Su precisión frente a dosis física o clínica permanece por determinar mediante el protocolo externo descrito.
