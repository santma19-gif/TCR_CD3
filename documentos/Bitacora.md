# 21 de Octubre del 2026

Por el momento no se va a iniciar nada que tenga que ver con el repositorio en GitHub, esto es más como una manera de ir adelantando de forma paralela el proyecto. 

--- 

Se crea la carpeta general de `Proyecto_TCR_CD3` en VS code, ademas se incluyen las subcarpetas para tener ordenado los datos:

- `analisis`
- `datos`
- `documentos`
- `interfaz`
- `pruebas`

Además se creo el ambiente en conda llamada `tcrcd3`, con Python 3.11, Streamlit, pandas y Biopyhton. 

Se creo una primera interfaz en `interfaz/app.py` funcionando: Una página con título, un botón para subir el FASTA, un filtro por gen y una tabla con variantes de ejemplo. 

``` python
import streamlit as st
import pandas as pd

st.title("Variantes del complejo TCR-CD3")

archivo = st.file_uploader("Sube tu archivo FASTA")

datos = pd.DataFrame({
    "Gen": ["CD3E", "CD3D", "CD247"],
    "Cambio": ["p.Arg174Ter", "p.Gly25Ser", "p.Gln70Leu"],
    "Consecuencia": ["nonsense", "missense", "missense"],
})

gen = st.selectbox("Filtrar por gen", ["Todos", "CD3D", "CD3E", "CD3G", "CD247"])
if gen != "Todos":
    datos = datos[datos["Gen"] == gen]

st.dataframe(datos)
```

Por el momento para correr la interfaz se hace: 
```bash 
conda activate tcrcd3 
streamlit run interfaz/app.py 
```
Para detener la interfaz se presiona Control + C en la terminal (cerrar la pestaña del navegador no la detiene).

# 22 de Octubre del 2026

Se empieza a construir el dataset de referencia: la información confiable (tomada de UniProt) que la herramienta va a consultar. Se decidió descargar los datos a mano desde la página de UniProt y no con código, para entender qué significa cada dato. Por ahora solo se hizo CD3ε (CD3E), para revisar que el procedimiento funcione antes de repetirlo con las demás.

1. Secuencia de la proteína (FASTA)

Pasos en uniprot.org:

Buscar el código P07766 y verificar que el título diga CD3E_HUMAN (Homo sapiens).
Botón Download → formato FASTA (canonical).
Guardar como datos/referencias/CD3E_proteina.fasta.
>sp|P07766|CD3E_HUMAN T-cell surface glycoprotein CD3 epsilon chain OS=Homo sapiens OX=9606 GN=CD3E PE=1 SV=2
MQSGTHWRVLGLCLLSVGVWGQDGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQ
HNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCE
NCMEMDVMSVATIVIVDICITGGLLLLVYYWSKNRKAKAKPVTRGAGAGGRQRGQNKERP
PPVPNPDYEPIRKGQRDLYSGLNQRRI

Cómo saber que es la canónica: el código P07766 no tiene sufijo (las isoformas son P07766-2, -3...), mide 207 aminoácidos y empieza con M (metionina).

Lectura del encabezado:

Parte	Significado
sp	Swiss-Prot: entrada revisada manualmente
P07766	Código de UniProt
CD3E_HUMAN	Nombre de la entrada
OS=Homo sapiens	Organismo
OX=9606	Código taxonómico del humano
GN=CD3E	Nombre del gen
PE=1	Existencia de la proteína confirmada experimentalmente
SV=2	Versión 2 de la secuencia
2. Regiones funcionales (CSV)

Las regiones están repartidas en tres secciones de la ficha de UniProt. En cada una se busca la tabla Features (columnas Type, Position(s) y Description):

Sección de UniProt	Qué se encontró
Subcellular Location	Extracelular (23–126), transmembrana (127–152), citoplasmático (153–207)
PTM/Processing	Péptido señal (1–22), cadena madura (23–207), puente disulfuro (C49↔C98)
Family & Domains	Dominio Ig-like (32–112), región de unión a NUMB (175–192), ITAM (178–205), región rica en prolinas (179–186)

Se creó el archivo datos/regiones_funcionales/regiones.csv:

gen,uniprot,tipo,descripcion,inicio,fin
CD3E,P07766,peptido_senal,,1,22
CD3E,P07766,dominio_topologico,extracelular,23,126
CD3E,P07766,dominio,Ig-like,32,112
CD3E,P07766,transmembrana,helice,127,152
CD3E,P07766,dominio_topologico,citoplasmatico,153,207
CD3E,P07766,region,union_NUMB,175,192
CD3E,P07766,dominio,ITAM,178,205
CD3E,P07766,region,rica_en_prolinas,179,186

Decisiones tomadas:

No se incluyó la fila Chain (23–207), porque es "todo menos el péptido señal" y ya está cubierta por las otras regiones.
No se incluyó la región Disordered (161–207), porque es una predicción automática.
El puente disulfuro y los residuos modificados no son regiones sino posiciones puntuales. Quedan pendientes para un archivo aparte (sitios_importantes.csv).
Las regiones se superponen (por ejemplo, el ITAM está dentro del dominio citoplasmático). Esto es normal, y significa que una variante puede caer en varias regiones a la vez. Hay que tenerlo en cuenta al definir el "molde" de la variante.
3. Conceptos aprendidos

PTM (post-translational modifications): modificaciones postraduccionales, es decir, los cambios que sufre la proteína después de sintetizarse (cortes, puentes disulfuro, fosforilaciones, etc.).

Péptido señal (1–22): una "etiqueta" que dirige la proteína al retículo endoplásmico para que termine en la membrana, y que luego se corta. También define la orientación: el extremo N-terminal queda afuera de la célula y la cola con el ITAM queda adentro. Es importante para el proyecto por dos razones:

Una variante en el péptido señal puede impedir que el complejo llegue a la membrana, aunque esa parte no esté en la proteína madura. No se deben descartar variantes en 1–22.
Numeración: ClinVar y la nomenclatura HGVS cuentan desde la metionina 1 (precursor), pero algunos artículos cuentan desde la proteína madura (sin las 22 posiciones). La herramienta usará siempre la numeración del precursor (desde Met1).

ITAM (Immunoreceptor Tyrosine-based Activation Motif): el "interruptor" que transmite la señal del TCR al interior del linfocito. Tiene el patrón YxxL/I – (6 a 8 aa) – YxxL/I. En CD3ε:

...D Y188 E P I191 R K G Q R D L Y199 S G L202 N...

Cuando el TCR reconoce un antígeno, la enzima Lck fosforila las dos tirosinas. A ellas se une ZAP-70 y se activa el linfocito. En el complejo completo hay 10 ITAM: uno en cada CD3ε (×2), CD3δ y CD3γ, y tres en cada CD3ζ (×2). Las variantes que cambien una de las tirosinas del ITAM deben marcarse con prioridad alta.

4. Primera visualización de las proteínas

Para revisar que el CSV quedó bien antes de agregar más proteínas, se creó interfaz/ver_proteinas.py. Dibuja cada proteína como barras horizontales con sus regiones en colores (como la barra de UniProt), y debajo muestra la tabla del CSV.

Se instaló la librería de gráficos (en el ambiente tcrcd3):

bash
pip install plotly

Para correrla (desde la carpeta principal del proyecto, porque el script busca el CSV en datos/regiones_funcionales/):

bash
conda activate tcrcd3
streamlit run interfaz/ver_proteinas.py

Qué revisar en el dibujo:

Topología: las cuatro regiones deben quedar seguidas, sin huecos ni superposiciones, de 1 a 207.
Dominios: el Ig-like debe caer dentro de la zona extracelular, y el ITAM dentro de la citoplasmática.
Al pasar el mouse sobre cada barra se ven sus posiciones, para compararlas con UniProt.

Cuando se agreguen más proteínas al CSV, aparecen solas en el menú de selección sin cambiar el código.

5. Estructura 3D con AlphaFold (CD3ε)
¿Por qué AlphaFold y no PDB?
Fuente	Qué es	Ventaja	Desventaja
AlphaFold	Modelos predichos por computador, uno por proteína	La numeración coincide exactamente con UniProt (y con nuestro CSV)	Muestra cada cadena sola, no el complejo
PDB	Estructuras experimentales (por ejemplo 6JXR, el complejo TCR-CD3 completo)	Se ve cómo encajan todas las cadenas	La numeración y los nombres de las cadenas no siempre coinciden con UniProt

Se decidió empezar con AlphaFold. PDB queda como mejora opcional.

Búsqueda en AlphaFold
Entrar a alphafold.ebi.ac.uk y buscar P07766.
Aparecen dos tipos de resultado:
Monomer, AF-P07766-F1: solo CD3ε. Es el que se usó. Longitud 207 (coincide con el FASTA) y pLDDT promedio de 73.06.
Heterodimer CD3ε + CD3γ: modelo del par ε/γ juntos (dataset de NVIDIA). Más cercano a la biología, pero queda para después.
Cómo leer la página del modelo

pLDDT: es la confianza de AlphaFold en cada aminoácido, de 0 a 100.

Color	pLDDT	Confianza
Azul oscuro	> 90	Muy alta
Azul claro	70–90	Alta
Amarillo	50–70	Baja
Naranja	< 50	Muy baja

Distribución en CD3ε: 28.5% muy alta, 28.5% alta, 28% baja y 15% muy baja.

Regiones reconocidas en el modelo 3D:

Lo que se ve	Color	Región del CSV
"Ovillo" compacto con flechas planas (láminas β)	Azul oscuro	Ig-like (32–112)
Hélice larga	Azul	Transmembrana (127–152)
Hélice suelta pegada al dominio	Amarillo	Péptido señal (1–22), que en la proteína real se corta
Cola larga y suelta	Amarillo/naranja	Citoplasmática (153–207), con el ITAM

PAE (Predicted Aligned Error): el cuadro verde. Indica qué tan seguro está AlphaFold de la posición de una parte respecto a otra (verde oscuro = seguro, blanco = inseguro). La cola citoplasmática (~160–207) sale en blanco: AlphaFold no sabe dónde va, y la forma de "fideo" es solo una de muchas posiciones posibles.

Conclusión: una variante en el Ig-like o en la transmembrana se puede interpretar con la estructura 3D, por ejemplo viendo si queda enterrada o en la superficie. Una variante en la cola citoplasmática no: ahí lo que importa es el CSV (si cae en el ITAM o en una de sus tirosinas). Esto no es un error de AlphaFold, sino que esas colas son naturalmente desordenadas.

El archivo PDB

Se descargó en formato PDB y se guardó como datos/estructuras/AF-P07766-F1-model_v6.pdb. Es un archivo de texto.

Encabezado:

Palabra	Qué contiene
HEADER	Fecha del archivo
TITLE	Qué es (predicción de AlphaFold para CD3ε)
COMPND	Molécula y cadena (solo hay una, la cadena A)
SOURCE	Organismo (Homo sapiens, taxid 9606)
REMARK	Notas y autores. Hay que citar el artículo de AlphaFold: Jumper et al. (2021)

Líneas ATOM: cada una es un átomo.

ATOM      2  CA  MET A   1   -48.204  11.585  23.268  1.00  44.97   C
Columna	Valor	Significado
Nombre del átomo	CA	Carbono alfa
Aminoácido	MET	Metionina
Cadena	A	La única cadena
Número de residuo	1	Posición en la proteína (misma numeración del CSV)
x, y, z	-48.204 11.585 23.268	Coordenadas en Ångström
Ocupación	1.00	Siempre 1.00 en AlphaFold
B-factor	44.97	Aquí AlphaFold guarda el pLDDT
Elemento	C	Carbono

Observaciones:

Los residuos 1, 2 y 3 son MET-GLN-SER (MQS...), igual que el FASTA. El FASTA, el CSV y el modelo 3D usan la misma numeración.
Cada aminoácido tiene varios átomos: N, CA, C y O forman el esqueleto, y el resto (CB, CG, SD...) es la cadena lateral.
Todos los átomos de un mismo aminoácido tienen el mismo pLDDT. Los primeros (44.97, 55.12, 59.91) son bajos porque están en el péptido señal.
Visor 3D en Streamlit

Se instaló la librería de dibujo 3D (en el ambiente tcrcd3):

bash
pip install py3Dmol

Se creó interfaz/ver_estructura.py, que:

Lee el archivo PDB como texto.
Saca el pLDDT de cada aminoácido: recorre las líneas ATOM, se queda con los CA y lee el B-factor. En el formato PDB cada dato está en columnas fijas, por eso se usa linea[60:66].
Muestra la estructura en 3D coloreada con los 4 colores de AlphaFold.
Muestra un gráfico del pLDDT por posición.

Para correrlo (desde la carpeta principal del proyecto; necesita internet la primera vez):

bash
conda activate tcrcd3
streamlit run interfaz/ver_estructura.py

En el visor se arrastra con el mouse para rotar y se usa la rueda o el trackpad para acercar.

Resultado: el visor se ve igual que en la página de AlphaFold, y el gráfico muestra la confianza alta en el Ig-like y la transmembrana, y baja en el péptido señal y la cola. Las bajadas a ~70 dentro del Ig-like probablemente son los lazos (loops) entre las láminas β.

# 23 de Octubre del 2026

Se completó el dataset de regiones con las otras tres proteínas y el visor 3D se adaptó para las cuatro.

1. Isoformas: siempre la canónica

Algunas proteínas (CD3D, CD247) tienen isoformas: versiones distintas de la proteína que salen del mismo gen por splicing alternativo (al procesar el ARN mensajero se incluyen u omiten algunos exones). UniProt elige una como canónica, y todas las anotaciones de la ficha se refieren a ella.

Al descargar en UniProt:

FASTA (canonical): solo la canónica. ✅ Esta es la que se usa.
FASTA (canonical & isoform): todas juntas. ❌

Cómo comprobar que es la canónica: el archivo tiene un solo encabezado (>sp|P04234|CD3D_HUMAN, sin -2), y la longitud coincide con la de la sección Sequence de UniProt.

Decisión del proyecto: usar siempre la isoforma canónica de UniProt. Cuando se descarguen las secuencias de nucleótidos de NCBI y las variantes de ClinVar, el transcrito debe corresponder a esta misma isoforma; si no, las posiciones no van a cuadrar.

2. Cómo decidir qué filas de UniProt usar

La tabla Features mezcla cosas muy distintas. Para cada fila se hacen dos preguntas:

Pregunta 1: ¿la posición es un rango o un punto?

Rango (1-22, 23-182): es una región → puede ir a regiones.csv.
Punto (52) o par (46↔87): es un sitio → va a sitios_importantes.csv (pendiente).

Pregunta 2 (solo para sitios): ¿la fuente es confiable? (columna Source)

UniProt: revisado por un curador, con publicaciones de respaldo → ✅ se guarda.
"Large scale data" (PTMeXchange, PRIDE): experimentos masivos de proteómica, menos seguros (las etiquetas Gold, Silver y Bronze indican cuánta evidencia tienen) → ❌ por ahora no.

Además:

Chain no se incluye: es "todo menos el péptido señal" y no agrega información.
Las anotaciones automáticas (Disordered, Compositional bias, con la etiqueta Automatic Annotation) no se incluyen.
3. Dominio, región y motivo
Tipo	Tamaño típico	Qué es	Ejemplo
Domain (dominio)	50 a 150+ aa	Unidad que se pliega sola en una estructura 3D estable	Ig-like de CD3ε (32–112)
Region (región)	Variable	Tramo con una función o propiedad anotada, que no necesariamente se pliega	Unión a NUMB, rica en prolinas
Motif (motivo)	2 a ~15 aa	Patrón corto de secuencia con función conocida	Di-leucina de CD3γ (153–154)

El ITAM es por su nombre un motivo (la "M" es Motif), pero UniProt lo clasifica como Domain.

En el CSV, la categoría de los motivos se escribe motif (en inglés, decisión del equipo). Hay que escribirla siempre igual (minúscula, sin espacios), porque el código trata como categorías distintas motif, Motif o motivo.

4. Reglas de formato del CSV
Columna gen: usar el símbolo oficial del gen (HGNC), no el nombre común de la proteína. Por eso es CD247 y no CD3Z. Es el que usan ClinVar y NCBI.
Columna descripcion del ITAM: escribir ITAM, sin número (no ITAM_1). El código busca exactamente esa palabra para el color rojo y para encontrar las tirosinas. El orden por posición ya indica cuál es el 1, el 2 y el 3.
Ordenar las filas de cada proteína por posición, empezando por el péptido señal.
5. CSV completo

datos/regiones_funcionales/regiones.csv:

gen,uniprot,tipo,descripcion,inicio,fin
CD3E,P07766,peptido_senal,,1,22
CD3E,P07766,dominio_topologico,extracelular,23,126
CD3E,P07766,dominio,Ig-like,32,112
CD3E,P07766,transmembrana,helice,127,152
CD3E,P07766,dominio_topologico,citoplasmatico,153,207
CD3E,P07766,region,union_NUMB,175,192
CD3E,P07766,dominio,ITAM,178,205
CD3E,P07766,region,rica_en_prolinas,179,186
CD3D,P04234,peptido_senal,,1,21
CD3D,P04234,dominio_topologico,extracelular,22,105
CD3D,P04234,transmembrana,helice,106,126
CD3D,P04234,dominio_topologico,citoplasmatico,127,171
CD3D,P04234,dominio,ITAM,138,166
CD3G,P09693,peptido_senal,,1,22
CD3G,P09693,dominio_topologico,extracelular,23,116
CD3G,P09693,dominio,Ig-like,37,94
CD3G,P09693,transmembrana,helice,117,137
CD3G,P09693,dominio_topologico,citoplasmatico,138,182
CD3G,P09693,dominio,ITAM,149,177
CD3G,P09693,motif,di_leucina,153,154
CD247,P20963,peptido_senal,,1,21
CD247,P20963,dominio_topologico,extracelular,22,30
CD247,P20963,transmembrana,helice,31,51
CD247,P20963,dominio_topologico,citoplasmatico,52,164
CD247,P20963,dominio,ITAM,61,89
CD247,P20963,dominio,ITAM,100,128
CD247,P20963,dominio,ITAM,131,159
Resumen por proteína
Proteína	Gen	UniProt	Longitud	ITAM	Particularidad
CD3ε	CD3E	P07766	207	1	Ig-like, región rica en prolinas, unión a NUMB
CD3δ	CD3D	P04234	171	1	Sin Ig-like anotado en UniProt; glicosilada
CD3γ	CD3G	P09693	182	1	Ig-like; motivo de di-leucina (internalización del TCR)
CD3ζ	CD247	P20963	164	3	Casi sin parte extracelular (9 aa); la cola es casi toda ITAM

Observaciones:

El péptido señal no mide siempre lo mismo (21 o 22 aa); cada proteína tiene el suyo.
CD3ε y CD3γ se parecen mucho (son parientes evolutivos): péptido señal, extracelular con Ig-like, transmembrana y un ITAM. CD3δ tiene la misma arquitectura, pero UniProt no anota su Ig-like.
CD3ζ es distinta: está especializada en señalizar, no en reconocer. Aporta 6 de los 10 ITAM del complejo.
6. Sitios importantes encontrados (para sitios_importantes.csv)

Solo los de fuente UniProt:

Proteína	Sitio	Qué es
CD3ε	C49 ↔ C98	Puente disulfuro (dentro del Ig-like)
CD3δ	C37 ↔ C73	Puente disulfuro
CD3δ	N38, N74	Glicosilación (N-linked)
CD3δ	Y149, Y160	Fosfotirosinas: son las dos tirosinas del ITAM (confirmación experimental)
CD3γ	C46 ↔ C87	Puente disulfuro (dentro del Ig-like)
CD3γ	N52, N92	Glicosilación (N-linked)
CD3γ	S145	Fosfoserina
CD3γ	S148	Fosfoserina por PKC: controla la internalización del TCR
CD247	58, 64, 72, 83, 111, 123, 142, 153	Modified residue: falta confirmar en la columna Description qué modificación es

CD3γ S148 y la numeración: en muchos artículos este sitio aparece como S126, porque cuentan desde la proteína madura (148 − 22 del péptido señal = 126). Igual con la di-leucina: L153/L154 en UniProt corresponden a L131/L132 en la literatura. Cuando PKC fosforila la serina, el motivo de di-leucina queda expuesto y la célula retira el TCR de la superficie.

CD247: las posiciones curadas caen en pares dentro de cada ITAM (72/83, 111/123, 142/153), así que muy probablemente son sus seis tirosinas.

7. Ajuste en ver_proteinas.py: carril de motivos

El motivo de di-leucina aparecía en un carril llamado "Otros". Se agregó una línea al diccionario CARRILES:

python
CARRILES = {
    "peptido_senal": "Topología",
    "dominio_topologico": "Topología",
    "transmembrana": "Topología",
    "dominio": "Dominios",
    "motif": "Motivos",
}

Las cuatro proteínas se revisaron en el dibujo de barras: la topología de cada una es continua y llega hasta su longitud. La pequeña separación oscura entre regiones vecinas no es un error de datos; es porque cada barra se dibuja de inicio a fin.

8. Visor 3D para las cuatro proteínas

Se descargaron de AlphaFold los modelos Monomer de las otras tres proteínas (revisando que la Sequence length coincidiera con el CSV) y se guardaron en datos/estructuras/:

Proteína	Archivo
CD3ε	AF-P07766-F1-model_v6.pdb
CD3δ	AF-P04234-F1-model_v6.pdb
CD3γ	AF-P09693-F1-model_v6.pdb
CD3ζ	AF-P20963-F1-model_v6.pdb

Cambios en interfaz/ver_estructura.py:

Menú de proteínas que sale del CSV, igual que en ver_proteinas.py. Debajo aparece el nombre, el gen y el código de UniProt.
Búsqueda automática del archivo: toma el código de UniProt del CSV y busca AF-<código>-F1-model_v*.pdb. El * es un comodín que acepta cualquier versión del modelo.
Avisos:
Si falta el modelo, aparece un recuadro rojo diciendo qué archivo falta y dónde debería estar.
Si el modelo no mide lo mismo que la proteína en el CSV, aparece una advertencia amarilla (posible isoforma equivocada).

Tirosinas del ITAM que debe encontrar el visor:

Proteína	Tirosinas esperadas
CD3ε	Y188, Y199
CD3δ	Y149, Y160
CD3γ	Y160, Y171
CD3ζ	Y72, Y83, Y111, Y123, Y142, Y153

Si coinciden, todo el flujo funciona de punta a punta: lo anotado a mano en UniProt, el CSV y el modelo 3D usan la misma numeración.

Pendientes
 Completar el dataset de regiones con CD3D, CD3G y CD247.
 Primer intento de estructura 3D con AlphaFold (CD3ε): visor coloreado por pLDDT.
 Visor 3D nivel 2: colorear la estructura con las regiones del CSV.
 Descargar los modelos de AlphaFold de las otras tres proteínas y adaptar el visor a las cuatro.
 Confirmar en el visor que las tirosinas encontradas coinciden con la tabla de la sección 8.
 Descargar los FASTA de proteína que falten en datos/referencias/.
 Confirmar en UniProt qué son los Modified residue de CD247.
 Crear sitios_importantes.csv con la tabla de la sección 6.
 Juntar COLORES y nombre_region (repetidos en ver_proteinas.py y ver_estructura.py) en un solo archivo compartido.
 Descargar las secuencias de nucleótidos de los cuatro genes desde NCBI (transcrito que corresponda a la isoforma canónica).
 Armar una lista de variantes reales conocidas desde ClinVar (datos/ejemplos/variantes_conocidas.csv).
 Acordar con Bioinformática el "molde" de la variante (analisis/variante.py). Tener en cuenta que una variante puede caer en varias regiones a la vez.
 Armar una lista de variantes reales conocidas desde ClinVar (datos/ejemplos/variantes_conocidas.csv).
 Acordar con Bioinformática el "molde" de la variante (analisis/variante.py).

