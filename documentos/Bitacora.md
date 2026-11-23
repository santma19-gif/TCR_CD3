# 21 de Octubre del 2026

Por el momento no se va a iniciar nada que tenga que ver con el repositorio en GitHub, esto es más como una manera de ir adelantando de forma paralela el proyecto.

---

Se crea la carpeta general de `Proyecto_TCR_CD3` en VS code, además se incluyen las subcarpetas para tener ordenado los datos:

- `analisis`
- `datos`
- `documentos`
- `interfaz`
- `pruebas`

Además se creó el ambiente en conda llamado `tcrcd3`, con Python 3.11, Streamlit, pandas y Biopython.

Se creó una primera interfaz en `interfaz/app.py` funcionando: Una página con título, un botón para subir el FASTA, un filtro por gen y una tabla con variantes de ejemplo.

```python
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

Para detener la interfaz se presiona `Control + C` en la terminal (cerrar la pestaña del navegador no la detiene).

---

# 22 de Octubre del 2026

Se empieza a construir el **dataset de referencia**: la información confiable (tomada de UniProt) que la herramienta va a consultar. Se decidió descargar los datos **a mano desde la página de UniProt** y no con código, para entender qué significa cada dato. Por ahora solo se hizo **CD3ε (CD3E)**, para revisar que el procedimiento funcione antes de repetirlo con las demás.

---

## 1. Secuencia de la proteína (FASTA)

Pasos en [uniprot.org](https://www.uniprot.org):

1. Buscar el código `P07766` y verificar que el título diga **CD3E_HUMAN** (*Homo sapiens*).
2. Botón **Download** → formato **FASTA (canonical)**.
3. Guardar como `datos/referencias/CD3E_proteina.fasta`.

```
>sp|P07766|CD3E_HUMAN T-cell surface glycoprotein CD3 epsilon chain OS=Homo sapiens OX=9606 GN=CD3E PE=1 SV=2
MQSGTHWRVLGLCLLSVGVWGQDGNEEMGGITQTPYKVSISGTTVILTCPQYPGSEILWQ
HNDKNIGGDEDDKNIGSDEDHLSLKEFSELEQSGYYVCYPRGSKPEDANFYLYLRARVCE
NCMEMDVMSVATIVIVDICITGGLLLLVYYWSKNRKAKAKPVTRGAGAGGRQRGQNKERP
PPVPNPDYEPIRKGQRDLYSGLNQRRI
```

Cómo saber que es la canónica: el código `P07766` no tiene sufijo (las isoformas son `P07766-2`, `-3`...), mide **207 aminoácidos** y empieza con `M` (metionina).

Lectura del encabezado:

| Parte | Significado |
|---|---|
| `sp` | Swiss-Prot: entrada revisada manualmente |
| `P07766` | Código de UniProt |
| `CD3E_HUMAN` | Nombre de la entrada |
| `OS=Homo sapiens` | Organismo |
| `OX=9606` | Código taxonómico del humano |
| `GN=CD3E` | Nombre del gen |
| `PE=1` | Existencia de la proteína confirmada experimentalmente |
| `SV=2` | Versión 2 de la secuencia |

---

## 2. Regiones funcionales (CSV)

Las regiones están repartidas en tres secciones de la ficha de UniProt. En cada una se busca la tabla **Features** (columnas *Type*, *Position(s)* y *Description*):

| Sección de UniProt | Qué se encontró |
|---|---|
| **Subcellular Location** | Extracelular (23–126), transmembrana (127–152), citoplasmático (153–207) |
| **PTM/Processing** | Péptido señal (1–22), cadena madura (23–207), puente disulfuro (C49↔C98) |
| **Family & Domains** | Dominio Ig-like (32–112), región de unión a NUMB (175–192), ITAM (178–205), región rica en prolinas (179–186) |

Se creó el archivo `datos/regiones_funcionales/regiones.csv`:

```
gen,uniprot,tipo,descripcion,inicio,fin
CD3E,P07766,peptido_senal,,1,22
CD3E,P07766,dominio_topologico,extracelular,23,126
CD3E,P07766,dominio,Ig-like,32,112
CD3E,P07766,transmembrana,helice,127,152
CD3E,P07766,dominio_topologico,citoplasmatico,153,207
CD3E,P07766,region,union_NUMB,175,192
CD3E,P07766,dominio,ITAM,178,205
CD3E,P07766,region,rica_en_prolinas,179,186
```

Decisiones tomadas:

- **No se incluyó** la fila *Chain* (23–207), porque es "todo menos el péptido señal" y ya está cubierta por las otras regiones.
- **No se incluyó** la región *Disordered* (161–207), porque es una predicción automática.
- El **puente disulfuro** y los **residuos modificados** no son regiones sino posiciones puntuales. Quedan pendientes para un archivo aparte (`sitios_importantes.csv`).
- Las regiones **se superponen** (por ejemplo, el ITAM está dentro del dominio citoplasmático). Esto es normal, y significa que una variante puede caer en varias regiones a la vez. Hay que tenerlo en cuenta al definir el "molde" de la variante.

---

## 3. Conceptos aprendidos

**PTM (*post-translational modifications*):** modificaciones postraduccionales, es decir, los cambios que sufre la proteína después de sintetizarse (cortes, puentes disulfuro, fosforilaciones, etc.).

**Péptido señal (1–22):** una "etiqueta" que dirige la proteína al retículo endoplásmico para que termine en la membrana, y que luego se corta. También define la orientación: el extremo N-terminal queda afuera de la célula y la cola con el ITAM queda adentro. Es importante para el proyecto por dos razones:

- Una variante en el péptido señal puede impedir que el complejo llegue a la membrana, aunque esa parte no esté en la proteína madura. **No se deben descartar variantes en 1–22.**
- **Numeración:** ClinVar y la nomenclatura HGVS cuentan desde la metionina 1 (precursor), pero algunos artículos cuentan desde la proteína madura (sin las 22 posiciones). **La herramienta usará siempre la numeración del precursor (desde Met1).**

**ITAM (*Immunoreceptor Tyrosine-based Activation Motif*):** el "interruptor" que transmite la señal del TCR al interior del linfocito. Tiene el patrón `YxxL/I – (6 a 8 aa) – YxxL/I`. En CD3ε:

```
...D Y188 E P I191 R K G Q R D L Y199 S G L202 N...
```

Cuando el TCR reconoce un antígeno, la enzima Lck fosforila las dos tirosinas. A ellas se une ZAP-70 y se activa el linfocito. En el complejo completo hay **10 ITAM**: uno en cada CD3ε (×2), CD3δ y CD3γ, y tres en cada CD3ζ (×2). **Las variantes que cambien una de las tirosinas del ITAM deben marcarse con prioridad alta.**

---

## 4. Primera visualización de las proteínas

Para revisar que el CSV quedó bien antes de agregar más proteínas, se creó `interfaz/ver_proteinas.py`. Dibuja cada proteína como barras horizontales con sus regiones en colores (como la barra de UniProt), y debajo muestra la tabla del CSV.

Se instaló la librería de gráficos (en el ambiente `tcrcd3`):

```bash
pip install plotly
```

Para correrla (desde la carpeta principal del proyecto, porque el script busca el CSV en `datos/regiones_funcionales/`):

```bash
conda activate tcrcd3
streamlit run interfaz/ver_proteinas.py
```

Qué revisar en el dibujo:

- **Topología:** las cuatro regiones deben quedar seguidas, sin huecos ni superposiciones, de 1 a 207.
- **Dominios:** el Ig-like debe caer dentro de la zona extracelular, y el ITAM dentro de la citoplasmática.
- Al pasar el mouse sobre cada barra se ven sus posiciones, para compararlas con UniProt.

Cuando se agreguen más proteínas al CSV, aparecen solas en el menú de selección sin cambiar el código.

---

## 5. Estructura 3D con AlphaFold (CD3ε)

### ¿Por qué AlphaFold y no PDB?

| Fuente | Qué es | Ventaja | Desventaja |
|---|---|---|---|
| **AlphaFold** | Modelos predichos por computador, uno por proteína | La numeración coincide exactamente con UniProt (y con nuestro CSV) | Muestra cada cadena sola, no el complejo |
| **PDB** | Estructuras experimentales (por ejemplo `6JXR`, el complejo TCR-CD3 completo) | Se ve cómo encajan todas las cadenas | La numeración y los nombres de las cadenas no siempre coinciden con UniProt |

Se decidió empezar con **AlphaFold**. PDB queda como mejora opcional.

### Búsqueda en AlphaFold

1. Entrar a [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk) y buscar `P07766`.
2. Aparecen dos tipos de resultado:
   - **Monomer, `AF-P07766-F1`:** solo CD3ε. **Es el que se usó.** Longitud 207 (coincide con el FASTA) y pLDDT promedio de 73.06.
   - **Heterodimer CD3ε + CD3γ:** modelo del par ε/γ juntos (dataset de NVIDIA). Más cercano a la biología, pero queda para después.

### Cómo leer la página del modelo

**pLDDT:** es la confianza de AlphaFold en cada aminoácido, de 0 a 100.

| Color | pLDDT | Confianza |
|---|---|---|
| Azul oscuro | > 90 | Muy alta |
| Azul claro | 70–90 | Alta |
| Amarillo | 50–70 | Baja |
| Naranja | < 50 | Muy baja |

Distribución en CD3ε: 28.5% muy alta, 28.5% alta, 28% baja y 15% muy baja.

Regiones reconocidas en el modelo 3D:

| Lo que se ve | Color | Región del CSV |
|---|---|---|
| "Ovillo" compacto con flechas planas (láminas β) | Azul oscuro | Ig-like (32–112) |
| Hélice larga | Azul | Transmembrana (127–152) |
| Hélice suelta pegada al dominio | Amarillo | Péptido señal (1–22), que en la proteína real se corta |
| Cola larga y suelta | Amarillo/naranja | Citoplasmática (153–207), con el ITAM |

**PAE (*Predicted Aligned Error*):** el cuadro verde. Indica qué tan seguro está AlphaFold de la posición de una parte respecto a otra (verde oscuro = seguro, blanco = inseguro). La cola citoplasmática (~160–207) sale en blanco: AlphaFold no sabe dónde va, y la forma de "fideo" es solo una de muchas posiciones posibles.

**Conclusión:** una variante en el Ig-like o en la transmembrana se puede interpretar con la estructura 3D, por ejemplo viendo si queda enterrada o en la superficie. Una variante en la cola citoplasmática no: ahí lo que importa es el CSV (si cae en el ITAM o en una de sus tirosinas). Esto no es un error de AlphaFold, sino que esas colas son naturalmente desordenadas.

### El archivo PDB

Se descargó en formato **PDB** y se guardó como `datos/estructuras/AF-P07766-F1-model_v6.pdb`. Es un archivo de texto.

**Encabezado:**

| Palabra | Qué contiene |
|---|---|
| `HEADER` | Fecha del archivo |
| `TITLE` | Qué es (predicción de AlphaFold para CD3ε) |
| `COMPND` | Molécula y cadena (solo hay una, la **cadena A**) |
| `SOURCE` | Organismo (*Homo sapiens*, taxid 9606) |
| `REMARK` | Notas y autores. Hay que citar el artículo de AlphaFold: **Jumper et al. (2021)** |

**Líneas `ATOM`:** cada una es un átomo.

```
ATOM      2  CA  MET A   1   -48.204  11.585  23.268  1.00  44.97   C
```

| Columna | Valor | Significado |
|---|---|---|
| Nombre del átomo | `CA` | Carbono alfa |
| Aminoácido | `MET` | Metionina |
| Cadena | `A` | La única cadena |
| Número de residuo | `1` | Posición en la proteína (**misma numeración del CSV**) |
| x, y, z | `-48.204 11.585 23.268` | Coordenadas en Ångström |
| Ocupación | `1.00` | Siempre 1.00 en AlphaFold |
| **B-factor** | **`44.97`** | **Aquí AlphaFold guarda el pLDDT** |
| Elemento | `C` | Carbono |

Observaciones:

- Los residuos 1, 2 y 3 son MET-GLN-SER (`MQS...`), igual que el FASTA. **El FASTA, el CSV y el modelo 3D usan la misma numeración.**
- Cada aminoácido tiene varios átomos: `N`, `CA`, `C` y `O` forman el **esqueleto**, y el resto (`CB`, `CG`, `SD`...) es la **cadena lateral**.
- Todos los átomos de un mismo aminoácido tienen el mismo pLDDT. Los primeros (44.97, 55.12, 59.91) son bajos porque están en el péptido señal.

### Visor 3D en Streamlit

Se instaló la librería de dibujo 3D (en el ambiente `tcrcd3`):

```bash
pip install py3Dmol
```

Se creó `interfaz/ver_estructura.py`, que:

1. Lee el archivo PDB como texto.
2. Saca el pLDDT de cada aminoácido: recorre las líneas `ATOM`, se queda con los `CA` y lee el B-factor. En el formato PDB cada dato está en columnas fijas, por eso se usa `linea[60:66]`.
3. Muestra la estructura en 3D coloreada con los 4 colores de AlphaFold.
4. Muestra un gráfico del pLDDT por posición.

Para correrlo (desde la carpeta principal del proyecto; necesita internet la primera vez):

```bash
conda activate tcrcd3
streamlit run interfaz/ver_estructura.py
```

En el visor se arrastra con el mouse para rotar y se usa la rueda o el trackpad para acercar.

**Resultado:** el visor se ve igual que en la página de AlphaFold, y el gráfico muestra la confianza alta en el Ig-like y la transmembrana, y baja en el péptido señal y la cola. Las bajadas a ~70 dentro del Ig-like probablemente son los lazos (*loops*) entre las láminas β.

### Visor 3D, nivel 2: colorear con el CSV

Se agregaron dos controles arriba del visor:

1. **"Colorear por":**
   - *Confianza (pLDDT)*: como antes.
   - *Regiones del CSV*: usa los mismos colores de `ver_proteinas.py`, para que el dibujo de barras y el modelo 3D coincidan. Abajo aparece una leyenda.
2. **"Resaltar las tirosinas del ITAM":** dibuja las tirosinas como "palitos" (sus átomos) con una etiqueta. Las posiciones no están escritas a mano: el código busca en el PDB los aminoácidos `TYR` que caen dentro de las posiciones del ITAM del CSV. En CD3ε encontró **Y188 y Y199**.

Detalles del código:

- **Orden de pintado:** como las regiones se superponen, primero se pinta la topología y después los dominios. Lo último que se pinta queda "encima", y por eso el ITAM (rojo) se ve sobre la zona citoplasmática (verde). Las filas de tipo `region` no se pintan, para no tapar el ITAM.
- **`setStyle` y `addStyle`:** `setStyle` **reemplaza** el estilo de esos aminoácidos y `addStyle` lo **agrega** sin borrar el anterior. Por eso las tirosinas se ven como palitos encima de la cinta.

Al comparar los dos modos se ve que la zona del ITAM coincide con la cola de baja confianza: para interpretar variantes ahí sirve el CSV, no la forma 3D.

---

# 23 de Octubre del 2026

Se completó el **dataset de regiones** con las otras tres proteínas y el visor 3D se adaptó para las cuatro.

---

## 1. Isoformas: siempre la canónica

Algunas proteínas (CD3D, CD247) tienen **isoformas**: versiones distintas de la proteína que salen del mismo gen por **splicing alternativo** (al procesar el ARN mensajero se incluyen u omiten algunos exones). UniProt elige una como **canónica**, y **todas las anotaciones de la ficha se refieren a ella**.

Al descargar en UniProt:

- **FASTA (canonical):** solo la canónica. ✅ **Esta es la que se usa.**
- **FASTA (canonical & isoform):** todas juntas. ❌

Cómo comprobar que es la canónica: el archivo tiene **un solo encabezado** (`>sp|P04234|CD3D_HUMAN`, sin `-2`), y la longitud coincide con la de la sección *Sequence* de UniProt.

**Decisión del proyecto: usar siempre la isoforma canónica de UniProt.** Cuando se descarguen las secuencias de nucleótidos de NCBI y las variantes de ClinVar, el transcrito debe corresponder a esta misma isoforma; si no, las posiciones no van a cuadrar.

---

## 2. Cómo decidir qué filas de UniProt usar

La tabla *Features* mezcla cosas muy distintas. Para cada fila se hacen dos preguntas:

**Pregunta 1: ¿la posición es un rango o un punto?**

- **Rango** (`1-22`, `23-182`): es una **región** → puede ir a `regiones.csv`.
- **Punto** (`52`) o **par** (`46↔87`): es un **sitio** → va a `sitios_importantes.csv` (pendiente).

**Pregunta 2 (solo para sitios): ¿la fuente es confiable?** (columna *Source*)

- **UniProt:** revisado por un curador, con publicaciones de respaldo → ✅ se guarda.
- **"Large scale data"** (PTMeXchange, PRIDE): experimentos masivos de proteómica, menos seguros (las etiquetas Gold, Silver y Bronze indican cuánta evidencia tienen) → ❌ por ahora no.

**Además:**

- **Chain** no se incluye: es "todo menos el péptido señal" y no agrega información.
- Las anotaciones automáticas (*Disordered*, *Compositional bias*, con la etiqueta *Automatic Annotation*) no se incluyen.

---

## 3. Dominio, región y motivo

| Tipo | Tamaño típico | Qué es | Ejemplo |
|---|---|---|---|
| **Domain** (dominio) | 50 a 150+ aa | Unidad que **se pliega sola** en una estructura 3D estable | Ig-like de CD3ε (32–112) |
| **Region** (región) | Variable | Tramo con una función o propiedad anotada, que no necesariamente se pliega | Unión a NUMB, rica en prolinas |
| **Motif** (motivo) | 2 a ~15 aa | **Patrón corto de secuencia** con función conocida | Di-leucina de CD3γ (153–154) |

El ITAM es por su nombre un motivo (la "M" es *Motif*), pero UniProt lo clasifica como *Domain*.

**En el CSV, la categoría de los motivos se escribe `motif`** (en inglés, decisión del equipo). Hay que escribirla siempre igual (minúscula, sin espacios), porque el código trata como categorías distintas `motif`, `Motif` o `motivo`.

---

## 4. Reglas de formato del CSV

- **Columna `gen`: usar el símbolo oficial del gen (HGNC)**, no el nombre común de la proteína. Por eso es `CD247` y no `CD3Z`. Es el que usan ClinVar y NCBI.
- **Columna `descripcion` del ITAM: escribir `ITAM`, sin número** (no `ITAM_1`). El código busca exactamente esa palabra para el color rojo y para encontrar las tirosinas. El orden por posición ya indica cuál es el 1, el 2 y el 3.
- Ordenar las filas de cada proteína por posición, empezando por el péptido señal.

---

## 5. CSV completo

`datos/regiones_funcionales/regiones.csv`:

```
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
```

### Resumen por proteína

| Proteína | Gen | UniProt | Longitud | ITAM | Particularidad |
|---|---|---|---|---|---|
| CD3ε | CD3E | P07766 | 207 | 1 | Ig-like, región rica en prolinas, unión a NUMB |
| CD3δ | CD3D | P04234 | 171 | 1 | Sin Ig-like anotado en UniProt; glicosilada |
| CD3γ | CD3G | P09693 | 182 | 1 | Ig-like; motivo de di-leucina (internalización del TCR) |
| CD3ζ | CD247 | P20963 | 164 | 3 | Casi sin parte extracelular (9 aa); la cola es casi toda ITAM |

Observaciones:

- **El péptido señal no mide siempre lo mismo** (21 o 22 aa); cada proteína tiene el suyo.
- **CD3ε y CD3γ se parecen mucho** (son parientes evolutivos): péptido señal, extracelular con Ig-like, transmembrana y un ITAM. CD3δ tiene la misma arquitectura, pero UniProt no anota su Ig-like.
- **CD3ζ es distinta:** está especializada en **señalizar**, no en reconocer. Aporta 6 de los 10 ITAM del complejo.

---

## 6. Sitios importantes encontrados (para `sitios_importantes.csv`)

Solo los de fuente **UniProt**:

| Proteína | Sitio | Qué es |
|---|---|---|
| CD3ε | C49 ↔ C98 | Puente disulfuro (dentro del Ig-like) |
| CD3δ | C37 ↔ C73 | Puente disulfuro |
| CD3δ | N38, N74 | Glicosilación (*N-linked*) |
| CD3δ | Y149, Y160 | Fosfotirosinas: **son las dos tirosinas del ITAM** (confirmación experimental) |
| CD3γ | C46 ↔ C87 | Puente disulfuro (dentro del Ig-like) |
| CD3γ | N52, N92 | Glicosilación (*N-linked*) |
| CD3γ | S145 | Fosfoserina |
| CD3γ | S148 | Fosfoserina por PKC: controla la **internalización del TCR** |
| CD247 | 58, 64, 72, 83, 111, 123, 142, 153 | *Modified residue*: falta confirmar en la columna *Description* qué modificación es |

**CD3γ S148 y la numeración:** en muchos artículos este sitio aparece como **S126**, porque cuentan desde la proteína madura (148 − 22 del péptido señal = 126). Igual con la di-leucina: L153/L154 en UniProt corresponden a **L131/L132** en la literatura. Cuando PKC fosforila la serina, el motivo de di-leucina queda expuesto y la célula retira el TCR de la superficie.

**CD247:** las posiciones curadas caen en pares dentro de cada ITAM (72/83, 111/123, 142/153), así que muy probablemente son sus seis tirosinas.

---

## 7. Ajuste en `ver_proteinas.py`: carril de motivos

El motivo de di-leucina aparecía en un carril llamado "Otros". Se agregó una línea al diccionario `CARRILES`:

```python
CARRILES = {
    "peptido_senal": "Topología",
    "dominio_topologico": "Topología",
    "transmembrana": "Topología",
    "dominio": "Dominios",
    "motif": "Motivos",
}
```

Las cuatro proteínas se revisaron en el dibujo de barras: la topología de cada una es continua y llega hasta su longitud. La pequeña separación oscura entre regiones vecinas no es un error de datos; es porque cada barra se dibuja de `inicio` a `fin`.

---

## 8. Visor 3D para las cuatro proteínas

Se descargaron de AlphaFold los modelos **Monomer** de las otras tres proteínas (revisando que la *Sequence length* coincidiera con el CSV) y se guardaron en `datos/estructuras/`:

| Proteína | Archivo |
|---|---|
| CD3ε | `AF-P07766-F1-model_v6.pdb` |
| CD3δ | `AF-P04234-F1-model_v6.pdb` |
| CD3γ | `AF-P09693-F1-model_v6.pdb` |
| CD3ζ | `AF-P20963-F1-model_v6.pdb` |

Cambios en `interfaz/ver_estructura.py`:

- **Menú de proteínas** que sale del CSV, igual que en `ver_proteinas.py`. Debajo aparece el nombre, el gen y el código de UniProt.
- **Búsqueda automática del archivo:** toma el código de UniProt del CSV y busca `AF-<código>-F1-model_v*.pdb`. El `*` es un comodín que acepta cualquier versión del modelo.
- **Avisos:**
  - Si falta el modelo, aparece un recuadro rojo diciendo qué archivo falta y dónde debería estar.
  - Si el modelo no mide lo mismo que la proteína en el CSV, aparece una advertencia amarilla (posible isoforma equivocada).

Tirosinas del ITAM que debe encontrar el visor:

| Proteína | Tirosinas esperadas |
|---|---|
| CD3ε | Y188, Y199 |
| CD3δ | Y149, Y160 |
| CD3γ | Y160, Y171 |
| CD3ζ | Y72, Y83, Y111, Y123, Y142, Y153 |

Si coinciden, todo el flujo funciona de punta a punta: lo anotado a mano en UniProt, el CSV y el modelo 3D usan la misma numeración.

---

# 23 de Noviembre del 2026

Se retoma el proyecto después de la revisión con el profesor. Se crea el repositorio en GitHub y se construye la primera tabla de variantes reales desde ClinVar.

---

## 1. Sugerencia del profesor

El profesor propone que la herramienta **no empiece desde un FASTA ya armado, sino desde datos crudos de secuenciación masiva**, como en un hospital donde se secuencia a un paciente con **Oxford Nanopore**:

1. Buscar datos públicos de Nanopore de humano (genoma o transcriptoma).
2. Preprocesarlos, quedarse solo con los 4 genes, mapearlos y detectar las variantes.
3. Pasar esas variantes a lo que ya se tiene para generar el reporte.

Advertencias del profesor:

- La calidad de Nanopore puede no ser tan buena para detectar polimorfismos. Conviene hacerlo **por mapeo** e incluir un **reporte del mapeo**.
- Lo que más le interesa es la **estandarización del análisis**: qué filtros se usan para extraer las variantes y cómo se seleccionan los datos.
- *Machine learning* es opcional. Su sugerencia: integrar la información genética con otros factores de riesgo (por ejemplo, la población de origen) usando un **método bayesiano**.

Flujo propuesto:

| Paso | Qué hace | Herramientas típicas |
|---|---|---|
| 1. Datos | Lecturas de Nanopore de humano | SRA/ENA (FASTQ) |
| 2. Control de calidad | Revisar y filtrar las lecturas | NanoPlot, chopper |
| 3. Mapeo | Alinear contra el genoma de referencia | minimap2, samtools |
| 4. Filtrado por genes | Quedarse con las lecturas de los 4 genes | `samtools view` + archivo BED |
| 5. Reporte del mapeo | Cobertura por gen | samtools, mosdepth |
| 6. Llamado de variantes | Detectar SNV e indels y filtrarlos | Clair3 o bcftools |
| 7. Traducción a proteína | Coordenada genómica → `p.Tyr188Cys` | VEP o código propio |
| 8. Anotación y reporte | Región, ITAM, ClinVar | Lo que ya existe |
| 9. (Opcional) Riesgo | Puntaje bayesiano | — |

Decisiones:

- **Usar genoma (WGS), no transcriptoma:** los genes CD3 solo se expresan en linfocitos T. En el genoma están en cualquier muestra.
- **Mapear contra el genoma completo (GRCh38)** y después extraer la región. Mapear solo contra los 4 genes genera variantes falsas por lecturas de otras partes del genoma. CD3D, CD3E y CD3G están juntos en 11q23.3 y CD247 en 1q24.2.
- Para desarrollar no hace falta descargar genomas completos: muchos proyectos publican BAM/CRAM indexados, y `samtools view <URL> <región>` descarga solo esa región.

---

## 2. División del trabajo: el VCF como "contrato"

El proyecto tiene dos partes: **Bioinformática** (procesar los datos de secuenciación, encontrar los genes, llamar las variantes) y **Programación** (automatización, anotación, reporte y *machine learning*). Para que Programación no tenga que esperar a Bioinformática, se acuerda que **la salida de Bioinformática es un archivo VCF**:

- **Entregable 1:** un VCF de los 4 genes en coordenadas **GRCh38**.
- **Entregable 2:** la cobertura por gen (para poder decir "no evaluable" si es baja, en vez de "sin variantes").
- **Archivo compartido:** BED con las coordenadas de los 4 genes.
- **Transcritos:** el MANE Select de cada gen.

Mientras tanto, Programación trabaja con un **VCF de prueba** hecho a partir de ClinVar. Cuando Bioinformática entregue, solo se cambia el archivo.

---

## 3. Avance de Bioinformática (7 de noviembre)

María José documentó en *Inicio de Proyecto 7 Nov*:

- Ambiente `bioinfo` (Python 3.11, pandas, Biopython, requests, BLAST).
- Descarga por código desde la API de UniProt del FASTA y el JSON de los 4 genes.
- Descarga de `variant_summary.txt.gz` de ClinVar y conteo con `awk` de las variantes de los 4 genes en GRCh38: **927 variantes** (68 patogénicas, 500 benignas, 307 VUS, 31 en conflicto, 21 sin clasificar).
- Pendiente de su parte: fenotipos desde OMIM (requiere registro y clave de API; mientras tanto, la columna `PhenotypeList` de ClinVar ya trae la enfermedad).

---

## 4. Repositorio en GitHub

Se creó el repositorio **público** [santma19-gif/TCR_CD3](https://github.com/santma19-gif/TCR_CD3), con licencia MIT. Por ahora solo trabaja Santiago; más adelante se puede agregar a María José como colaboradora (*Settings → Collaborators*), y en ese momento cada uno trabajará en su propia rama.

Se mantuvo la estructura de carpetas que ya existía (`datos/`, `analisis/`, `interfaz/`, `pruebas/`, `documentos/`), para no tener que cambiar las rutas de los scripts. Se agregó `resultados/` para lo que genere el código.

Archivos base en la raíz:

| Archivo | Para qué |
|---|---|
| `README.md` | Descripción, estructura, instalación y uso |
| `.gitignore` | Evita subir archivos pesados (`*.gz`, `*.fastq`, `*.bam`...), `resultados/`, `.vscode/` y `.env` (claves de API) |
| `environment.yml` | Un solo ambiente `tcrcd3` para todo el proyecto |

Comandos usados para crearlo:

```bash
git init
git branch -M main
git remote add origin https://github.com/santma19-gif/TCR_CD3.git
git add .
git commit -m "Estructura inicial del proyecto"
git push -u origin main
```

Para subir cambios de ahí en adelante:

```bash
git add .
git commit -m "Qué cambió"
git push
```

Lecciones:

- El archivo se debe llamar **`.gitignore`, con punto**. Al descargarlo quedó como `gitignore` y por eso aparecían `.DS_Store` y `.vscode/` en `git status`. Se corrigió con `mv gitignore .gitignore`.
- En GitHub no se suben los datos pesados (límite de 100 MB por archivo): se sube el comando que los descarga.

---

## 5. Variantes conocidas desde ClinVar (`analisis/variantes_conocidas.py`)

Se descargó ClinVar (437 MB, no se sube a GitHub):

```bash
curl -o datos/clinvar/variant_summary.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
```

El script, para cada variante de CD3D, CD3E, CD3G y CD247 en GRCh38:

1. Separa el nombre de ClinVar (`NM_000733.4(CD3E):c.520C>T (p.Arg174Ter)`) en transcrito, cambio en el ADN (`c.`), cambio en la proteína (`p.`) y posición.
2. Agrupa la clasificación en Patogénica, Benigna, VUS, Conflicto o Sin clasificar (mismo criterio que Bioinformática).
3. Clasifica la consecuencia: missense, nonsense, frameshift, sinónima, inframe_indel, pérdida del inicio, pérdida del stop. Las que no cambian la proteína se separan en splicing (±1/±2 del borde del exón), intrónica, UTR_5, UTR_3 y solo_genomica.
4. Cruza la posición con `regiones.csv` (puede caer en varias regiones) y marca si cambia una **tirosina del ITAM**.
5. Compara el aminoácido de ClinVar con el FASTA canónico (`revision_fasta`).

Genera, en `datos/ejemplos/`:

- `variantes_conocidas.csv`: la tabla anotada.
- `variantes_conocidas.vcf`: el **VCF de prueba** (927 variantes).

Para correrlo:

```bash
conda activate tcrcd3
python analisis/variantes_conocidas.py
```

### Resultados

| Gen | Patogénicas | Benignas | VUS | Conflicto | Sin clasificar | Total |
|---|---|---|---|---|---|---|
| CD3D | 17 | 120 | 82 | 5 | 6 | 230 |
| CD3E | 20 | 181 | 86 | 11 | 7 | 305 |
| CD3G | 22 | 73 | 69 | 13 | 5 | 182 |
| CD247 | 9 | 126 | 70 | 2 | 3 | 210 |
| **Total** | **68** | **500** | **307** | **31** | **21** | **927** |

Coincide con el conteo de Bioinformática.

| Consecuencia | Patogénicas | Benignas | VUS |
|---|---|---|---|
| nonsense | 22 | 0 | 3 |
| splicing | 19 | 0 | 3 |
| frameshift | 17 | 0 | 3 |
| pérdida del inicio | 4 | 0 | 3 |
| missense | 3 | 10 | **231** |
| intrónica | 2 | 298 | 19 |
| solo_genomica | 1 | 5 | 0 |
| sinónima | 0 | 178 | 4 |

- **Variantes que cambian una tirosina del ITAM: 7.**
- **Revisión contra el FASTA canónico: 293 `ok`, ningún desajuste.** Las variantes de CD3D y CD3E en ClinVar usan la misma numeración que el FASTA, `regiones.csv` y AlphaFold. Las 634 `sin_revisar` son las que no cambian la proteína, más las de CD3G y CD247 (faltan sus FASTA).
- **Transcritos que usa ClinVar:** `NM_000732.6` (CD3D, también `.4`), `NM_000733.4` (CD3E, también `.3`), `NM_000073.3` (CD3G) y `NM_198053.3` (CD247). Las versiones viejas no afectan la numeración. Falta confirmar en NCBI que sean los MANE Select.

### Qué dicen estos números

1. **La mayoría de las patogénicas destruyen la proteína** (nonsense, frameshift, splicing). Son las fáciles de interpretar.
2. **El aporte de la herramienta está en las 231 missense VUS.** Para esas, saber si caen en una tirosina del ITAM, en la transmembrana o en el Ig-like es justo lo que `regiones.csv` y el visor 3D pueden aportar.
3. **No hay datos para entrenar un modelo de *machine learning*:** solo hay 3 missense patogénicas y 10 benignas. Lo defendible es un **puntaje bayesiano o por reglas** (región, ITAM, frecuencia poblacional, conservación, criterios ACMG).
4. **19 de las 68 patogénicas son de *splicing*.** Si la herramienta solo mirara regiones de la proteína, se perdería casi un tercio de las patogénicas. El reporte debe anotar también si la variante cae en un sitio de *splicing*.

### Casos raros encontrados y cómo se resolvieron

- **`p.Ter208=` en CD3E** salía como `fuera_de_rango`: es un cambio sinónimo en el **codón de parada**, que está en la posición 208 (la proteína mide 207). Se corrigió la revisión para que lo cuente como `ok`.
- **6 variantes descritas solo en coordenadas genómicas** (`NC_000011.10:g.118343086G>A`): ClinVar no dice qué le pasa al transcrito. Se clasifican como `solo_genomica`. Una es **patogénica**: la deleción de 17 bases `NC_000011.10:g.118349744_118349760del` en CD3G. Se podrá ubicar (exón, intrón o splicing) cuando se tengan las coordenadas de los exones.
- Se agregó la columna `nombre_clinvar` al CSV para revisar casos raros sin volver a leer el archivo de 437 MB.

---

## Pendientes

### Programación

- [x] Crear el repositorio en GitHub.
- [x] Tabla de variantes conocidas desde ClinVar y VCF de prueba (`analisis/variantes_conocidas.py`).
- [ ] Mostrar las variantes en la interfaz: marcas sobre el dibujo de barras (color según el grupo clínico), filtros por consecuencia y resaltado en el visor 3D. Revisar primero qué hace `interfaz/ver_variantes.py` y documentarlo.
- [ ] Lector de VCF genérico (para cuando llegue el VCF de Bioinformática).
- [ ] Generar el reporte por variante.
- [ ] Esqueleto de la automatización (script o Snakemake) con los pasos de Bioinformática.
- [ ] Puntaje bayesiano o por reglas (en lugar de *machine learning*). Buscar tablas de prevalencia por población y frecuencias en gnomAD.
- [ ] Juntar `COLORES` y `nombre_region` (repetidos en `ver_proteinas.py` y `ver_estructura.py`) en un solo archivo compartido.

### Datos

- [ ] Descargar los FASTA de proteína de **CD3G** y **CD247** en `datos/referencias/`.
- [ ] Confirmar en NCBI que los transcritos de ClinVar son los **MANE Select**.
- [ ] Crear `sitios_importantes.csv` con la tabla del 23 de octubre, sección 6.
- [ ] Confirmar en UniProt qué son los *Modified residue* de CD247.
- [ ] Confirmar en el visor que las tirosinas encontradas coinciden con la tabla del 23 de octubre, sección 8.

### Bioinformática

- [ ] Acordar con María José el "contrato" del VCF (GRCh38, cobertura por gen, BED, transcritos).
- [ ] Buscar datos públicos de Oxford Nanopore de humano (genoma). Candidato para validar: HG002 de Genome in a Bottle.
- [ ] Definir y documentar los filtros del llamado de variantes (profundidad, calidad, frecuencia alélica).
- [ ] Coordenadas de exones y sitios de *splicing* de los 4 genes (archivo BED).
- [ ] Fenotipos desde OMIM.

### Por definir

- [ ] Unificar los ambientes `tcrcd3` y `bioinfo` en `environment.yml`.
- [ ] Confirmar con el profesor que no hay problema con que el repositorio sea público.