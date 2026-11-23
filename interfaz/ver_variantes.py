import math
import re
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from plotly.subplots import make_subplots

st.title("Variantes de ClinVar sobre las proteínas CD3")
st.caption("Lee `datos/ejemplos/variantes_conocidas.csv`, generado por "
           "`analisis/variantes_conocidas.py`.")

ARCHIVO_VARIANTES = Path("datos/ejemplos/variantes_conocidas.csv")

# ---------------------------------------------------------------
# 1. Leer las regiones y elegir la proteína
# ---------------------------------------------------------------
regiones = pd.read_csv("datos/regiones_funcionales/regiones.csv")

gen = st.selectbox("Proteína", regiones["gen"].unique())
regiones_gen = regiones[regiones["gen"] == gen]

# ---------------------------------------------------------------
# 2. Leer las variantes ya procesadas
# ---------------------------------------------------------------
if not ARCHIVO_VARIANTES.exists():
    st.error(f"No se encontró `{ARCHIVO_VARIANTES}`. Primero corre:\n\n"
             "`python analisis/variantes_conocidas.py`")
    st.stop()

todas = pd.read_csv(ARCHIVO_VARIANTES)
todas_gen = todas[todas["gen"] == gen]

# El review status viene en texto; lo convertimos a estrellas
ESTRELLAS = {
    "no assertion criteria provided": 0,
    "criteria provided, single submitter": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}

# Altura del "chupetín" según el grupo clínico
NIVEL = {
    "Benigna": 1,
    "Sin_clasificar": 1.5,
    "VUS": 2,
    "Conflicto": 2.5,
    "Patogenica": 3,
}

# Tipos de variante que cambian un aminoácido concreto (se pueden dibujar)
CON_POSICION = ["missense", "nonsense", "frameshift", "inframe_indel",
                "perdida_inicio", "perdida_stop", "sinonima"]

variantes = pd.DataFrame({
    "cambio_p": todas_gen["hgvs_p"],
    "posicion_aa": todas_gen["posicion_proteina"],
    "aa_ref": todas_gen["aa_ref"],
    "consecuencia": todas_gen["consecuencia"],
    "grupo": todas_gen["grupo"],
    "clasificacion": todas_gen["ClinicalSignificance"],
    "estrellas": todas_gen["ReviewStatus"].map(ESTRELLAS),
    "tirosina_ITAM": todas_gen["tirosina_ITAM"],
    "condicion": todas_gen["PhenotypeList"],
    "cambio_c": todas_gen["hgvs_c"],
    "clinvar_id": todas_gen["VariationID"],
})

# ---------------------------------------------------------------
# 3. Filtros
# ---------------------------------------------------------------
# Nombre que se muestra en pantalla para cada valor del CSV
ETIQUETAS_GRUPO = {
    "Patogenica": "Patogénica",
    "Conflicto": "En conflicto",
    "VUS": "Incierta (VUS)",
    "Sin_clasificar": "Sin clasificar",
    "Benigna": "Benigna",
}
ETIQUETAS_CONSECUENCIA = {
    "missense": "Missense",
    "nonsense": "Nonsense",
    "frameshift": "Frameshift",
    "inframe_indel": "Inserción/deleción sin corrimiento",
    "perdida_inicio": "Pérdida del inicio",
    "perdida_stop": "Pérdida del stop",
    "sinonima": "Sinónima",
}
GRUPOS_POR_DEFECTO = ["Patogenica", "Conflicto", "VUS"]


def casillas(titulo, etiquetas, marcadas):
    """Una casilla por opción en la barra lateral. Devuelve las opciones marcadas."""
    st.sidebar.subheader(titulo)
    return [valor for valor, texto in etiquetas.items()
            if st.sidebar.checkbox(texto, value=valor in marcadas, key=f"{titulo}_{valor}")]


st.sidebar.header("Filtros")
grupos_elegidos = casillas("Grupo clínico", ETIQUETAS_GRUPO, GRUPOS_POR_DEFECTO)
consecuencias_elegidas = casillas("Consecuencia", ETIQUETAS_CONSECUENCIA,
                                  [c for c in CON_POSICION if c != "sinonima"])
st.sidebar.subheader("Otros")
solo_itam = st.sidebar.checkbox("Solo las que cambian una tirosina del ITAM")

# Las variantes sin cambio en la proteína (splicing, intrónicas, UTR) no se
# pueden ubicar en el dibujo: se cuentan aparte
sin_posicion = variantes[variantes["posicion_aa"].isna()]
variantes = variantes.dropna(subset=["posicion_aa"])
variantes["posicion_aa"] = variantes["posicion_aa"].astype(int)

variantes = variantes[variantes["grupo"].isin(grupos_elegidos)
                      & variantes["consecuencia"].isin(consecuencias_elegidas)]
if solo_itam:
    variantes = variantes[variantes["tirosina_ITAM"]]

if variantes.empty:
    st.warning("Ninguna variante cumple los filtros elegidos.")
    st.stop()

# ---------------------------------------------------------------
# 4. Relacionar cada variante con su región
# ---------------------------------------------------------------
def nombre_region(fila):
    """El nombre de cada fila del CSV (igual que en ver_proteinas.py)."""
    if fila["tipo"] in ("peptido_senal", "transmembrana"):
        return fila["tipo"]
    return fila["descripcion"]


def regiones_de(posicion):
    """Todas las regiones del CSV que contienen esa posición."""
    dentro = regiones_gen[(regiones_gen["inicio"] <= posicion) & (regiones_gen["fin"] >= posicion)]
    return ", ".join(nombre_region(fila) for _, fila in dentro.iterrows())


variantes["regiones"] = variantes["posicion_aa"].apply(regiones_de)

# ---------------------------------------------------------------
# 5. Dibujo: chupetines arriba, regiones abajo
# ---------------------------------------------------------------
COLORES = {
    "peptido_senal": "#9e9e9e",
    "extracelular": "#4c9be8",
    "transmembrana": "#f2b134",
    "citoplasmatico": "#7cc47f",
    "Ig-like": "#1f5fa8",
    "ITAM": "#d64545",
}
COLORES_CONSECUENCIA = {
    "missense": "#e8a33d",
    "nonsense": "#b5179e",
    "frameshift": "#7209b7",
    "inframe_indel": "#3a86ff",
    "perdida_inicio": "#d00000",
    "perdida_stop": "#6a994e",
    "sinonima": "#bbbbbb",
}
CARRILES = {
    "peptido_senal": "Topología",
    "dominio_topologico": "Topología",
    "transmembrana": "Topología",
    "dominio": "Dominios",
    "motif": "Motivos",
}

# Variantes en la misma posición y grupo se separan un poco en altura,
# para que no queden una encima de otra
variantes = variantes.sort_values(["posicion_aa", "grupo", "cambio_p"])
variantes["altura"] = (variantes["grupo"].map(NIVEL)
                       + 0.12 * variantes.groupby(["posicion_aa", "grupo"]).cumcount())

fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                    row_heights=[0.55, 0.45], vertical_spacing=0.04)

# 5a. Los "palitos" de cada chupetín (líneas verticales grises)
palitos_x, palitos_y = [], []
for _, v in variantes.iterrows():
    palitos_x += [v["posicion_aa"], v["posicion_aa"], None]
    palitos_y += [0, v["altura"], None]
fig.add_trace(go.Scatter(x=palitos_x, y=palitos_y, mode="lines",
                         line=dict(color="#cccccc", width=1),
                         hoverinfo="skip", showlegend=False), row=1, col=1)

# 5b. Las "cabezas": color según la consecuencia; estrella si toca una tirosina del ITAM
for consecuencia, grupo in variantes.groupby("consecuencia"):
    fig.add_trace(go.Scatter(
        x=grupo["posicion_aa"],
        y=grupo["altura"],
        mode="markers",
        marker=dict(size=[16 if t else 10 for t in grupo["tirosina_ITAM"]],
                    symbol=["star" if t else "circle" for t in grupo["tirosina_ITAM"]],
                    color=COLORES_CONSECUENCIA.get(consecuencia, "#888888"),
                    line=dict(width=1, color="white")),
        name=ETIQUETAS_CONSECUENCIA.get(consecuencia, consecuencia),
        customdata=grupo[["cambio_p", "clasificacion", "estrellas", "regiones", "condicion"]],
        hovertemplate=("<b>%{customdata[0]}</b><br>%{customdata[1]} "
                       "(%{customdata[2]} ★)<br>Región: %{customdata[3]}"
                       "<br>%{customdata[4]}<extra></extra>"),
    ), row=1, col=1)

# 5c. Las barras de regiones (igual que en ver_proteinas.py)
for _, fila in regiones_gen.iterrows():
    nombre = nombre_region(fila)
    carril = fila["descripcion"] if fila["tipo"] == "region" else CARRILES.get(fila["tipo"], "Otros")
    fig.add_trace(go.Bar(
        y=[carril], base=[fila["inicio"]], x=[fila["fin"] - fila["inicio"]],
        orientation="h", text=nombre, textposition="inside", insidetextanchor="middle",
        marker_color=COLORES.get(nombre, "#bdbdbd"), showlegend=False,
        hovertemplate=f"<b>{nombre}</b><br>Posición {fila['inicio']}–{fila['fin']}<extra></extra>",
    ), row=2, col=1)

numero_carriles = regiones_gen.apply(
    lambda f: f["descripcion"] if f["tipo"] == "region" else CARRILES.get(f["tipo"], "Otros"), axis=1
).nunique()

fig.update_layout(barmode="overlay", height=380 + 40 * numero_carriles,
                  margin=dict(t=20), legend_title_text="Consecuencia")
fig.update_yaxes(row=1, col=1, range=[0, variantes["altura"].max() + 0.5],
                 tickvals=list(NIVEL.values()), ticktext=[ETIQUETAS_GRUPO[g] for g in NIVEL])
fig.update_yaxes(row=2, col=1, autorange="reversed")
fig.update_xaxes(row=2, col=1, title_text="Posición (aminoácido)",
                 range=[0, regiones_gen["fin"].max() + 5])

st.plotly_chart(fig)
st.caption("★ = la variante cambia una tirosina del ITAM.")

# ---------------------------------------------------------------
# 6. Resumen y tablas
# ---------------------------------------------------------------
st.write(f"**{len(variantes)}** variantes ubicadas en la proteína con los filtros elegidos.")

if not sin_posicion.empty:
    with st.expander(f"{len(sin_posicion)} variantes de {gen} sin cambio en la proteína "
                     "(no se pueden dibujar)"):
        st.write("Variantes de *splicing*, intrónicas, en las UTR o descritas solo en "
                 "coordenadas genómicas. Las de *splicing* pueden ser patogénicas aunque "
                 "no cambien un aminoácido.")
        st.dataframe(pd.crosstab(sin_posicion["consecuencia"], sin_posicion["grupo"]))
        st.dataframe(sin_posicion[["cambio_c", "consecuencia", "grupo", "clasificacion",
                                   "condicion", "clinvar_id"]], hide_index=True)

st.subheader("Variantes por región")
conteo = []
for _, fila in regiones_gen.iterrows():
    dentro = variantes[(variantes["posicion_aa"] >= fila["inicio"]) & (variantes["posicion_aa"] <= fila["fin"])]
    conteo.append({"región": nombre_region(fila),
                   "posiciones": f"{fila['inicio']}–{fila['fin']}",
                   "variantes": len(dentro)})
st.dataframe(pd.DataFrame(conteo), hide_index=True)

st.subheader("Tabla de variantes")
st.dataframe(variantes.drop(columns=["altura"]), hide_index=True)

# ===============================================================
# 7. VER UNA VARIANTE EN 3D
# ===============================================================
st.header("Ver una variante en la estructura 3D")

# ---------------------------------------------------------------
# 7a. Elegir la variante
# ---------------------------------------------------------------
variantes_ordenadas = variantes.reset_index(drop=True)
etiquetas = [f"{'★ ' if v['tirosina_ITAM'] else ''}{v['cambio_p']}  ·  {v['consecuencia']}  ·  {v['grupo']}"
             for _, v in variantes_ordenadas.iterrows()]
eleccion = st.selectbox("Variante", range(len(etiquetas)), format_func=lambda i: etiquetas[i])
variante = variantes_ordenadas.iloc[eleccion]
posicion = int(variante["posicion_aa"])   # int() normal de Python, para que py3Dmol lo acepte

# ---------------------------------------------------------------
# 7b. Leer el modelo de AlphaFold (igual que en ver_estructura.py)
# ---------------------------------------------------------------
codigo_uniprot = regiones_gen["uniprot"].iloc[0]
encontrados = sorted(Path("datos/estructuras").glob(f"AF-{codigo_uniprot}-F1-model_v*.pdb"))
if not encontrados:
    st.warning(f"No se encontró el modelo de AlphaFold de {gen} en `datos/estructuras/`.")
    st.stop()
pdb = encontrados[-1].read_text()

# Guardamos, por cada aminoácido: su nombre, su pLDDT, la posición de su
# carbono alfa (CA) y la de todos sus átomos (para medir distancias).
aminoacidos, plddt, carbono_alfa, atomos = {}, {}, {}, {}
for linea in pdb.splitlines():
    if not linea.startswith("ATOM"):
        continue
    residuo = int(linea[22:26])
    xyz = (float(linea[30:38]), float(linea[38:46]), float(linea[46:54]))
    atomos.setdefault(residuo, []).append(xyz)
    if linea[12:16].strip() == "CA":
        aminoacidos[residuo] = linea[17:20]
        plddt[residuo] = float(linea[60:66])
        carbono_alfa[residuo] = xyz

longitud = max(aminoacidos)

# ---------------------------------------------------------------
# 7c. Casos que no se pueden mostrar en 3D
# ---------------------------------------------------------------
if posicion > longitud:
    # p.Ter208=: el codón de parada está justo después del último aminoácido
    st.info(f"{variante['cambio_p']} está en el codón de parada (posición {posicion}); "
            "no corresponde a ningún aminoácido del modelo.")
    st.stop()

# ---------------------------------------------------------------
# 7d. Comprobar que el aminoácido de ClinVar es el mismo del modelo
# ---------------------------------------------------------------
original = str(variante["aa_ref"]).upper()
if aminoacidos.get(posicion) != original:
    st.error(f"ClinVar dice que en la posición {posicion} hay {original}, pero en el modelo "
             f"hay {aminoacidos.get(posicion)}. La numeración no coincide.")
    st.stop()

# Propiedades básicas de cada aminoácido, para explicar el cambio
PROPIEDADES = {
    "ALA": "Alanina, pequeña e hidrofóbica", "ARG": "Arginina, carga positiva",
    "ASN": "Asparagina, polar", "ASP": "Aspartato, carga negativa",
    "CYS": "Cisteína, puede formar puentes disulfuro", "GLN": "Glutamina, polar",
    "GLU": "Glutamato, carga negativa", "GLY": "Glicina, la más pequeña y flexible",
    "HIS": "Histidina, puede tener carga positiva", "ILE": "Isoleucina, hidrofóbica",
    "LEU": "Leucina, hidrofóbica", "LYS": "Lisina, carga positiva",
    "MET": "Metionina, hidrofóbica", "PHE": "Fenilalanina, hidrofóbica y aromática",
    "PRO": "Prolina, rígida; rompe hélices y láminas", "SER": "Serina, polar",
    "THR": "Treonina, polar", "TRP": "Triptófano, grande y aromático",
    "TYR": "Tirosina, aromática; se fosforila en los ITAM", "VAL": "Valina, hidrofóbica",
}


def distancia(a, b):
    return math.dist(a, b)


# ---------------------------------------------------------------
# 7e. Funciones para dibujar
# ---------------------------------------------------------------
def crear_vista():
    """Un visor 3D con la proteína pintada con los colores de sus regiones."""
    vista = py3Dmol.view(width=430, height=430)
    vista.addModel(pdb, "pdb")
    vista.setBackgroundColor("white")
    vista.setStyle({}, {"cartoon": {"color": "#bdbdbd"}})
    # Primero la topología y luego los dominios, para que queden encima
    for tipo in ["peptido_senal", "dominio_topologico", "transmembrana", "dominio"]:
        for _, fila in regiones_gen[regiones_gen["tipo"] == tipo].iterrows():
            posiciones = list(range(fila["inicio"], fila["fin"] + 1))
            vista.setStyle({"resi": posiciones},
                           {"cartoon": {"color": COLORES.get(nombre_region(fila), "#bdbdbd")}})
    return vista


def marcar_variante(vista, estilo):
    """Dibuja el aminoácido de la variante en palitos magenta, con su nombre."""
    vista.addStyle({"resi": posicion}, {"stick": estilo})
    vista.addLabel(variante["cambio_p"], {"fontSize": 14, "fontColor": "black",
                                          "backgroundColor": "white", "backgroundOpacity": 0.8},
                   {"resi": posicion, "atom": "CA"})


def mostrar_leyenda(incluir_perdida=False):
    """Cuadritos de color con el nombre de cada región de esta proteína."""
    partes = []
    for _, fila in regiones_gen.iterrows():
        nombre = nombre_region(fila)
        if fila["tipo"] == "region" or nombre not in COLORES:
            continue   # las 'region' y los 'motif' no se pintan en 3D
        partes.append(f"<span style='color:{COLORES[nombre]}; font-size:20px'>■</span> "
                      f"{nombre} ({fila['inicio']}–{fila['fin']})")
    if incluir_perdida:
        partes.append("<span style='color:#d0d0d0; font-size:20px'>■</span> parte perdida")
    partes.append("<span style='color:#c51b8a; font-size:20px'>■</span> variante")
    st.markdown("<br>".join(partes), unsafe_allow_html=True)


def mapa_lineal(clave, corte=None):
    """Barra de regiones (como en ver_proteinas.py) con la variante marcada.
    Si 'corte' tiene un número, sombrea desde ahí hasta el final.
    'clave' es un nombre único: Streamlit lo pide cuando hay dos gráficos iguales."""
    mapa = go.Figure()
    for _, fila in regiones_gen.iterrows():
        if fila["tipo"] == "region":
            continue
        nombre = nombre_region(fila)
        mapa.add_trace(go.Bar(
            y=[CARRILES.get(fila["tipo"], "Otros")], base=[fila["inicio"]],
            x=[fila["fin"] - fila["inicio"]], orientation="h",
            marker_color=COLORES.get(nombre, "#bdbdbd"), showlegend=False,
            hovertemplate=f"<b>{nombre}</b><br>Posición {fila['inicio']}–{fila['fin']}<extra></extra>",
        ))
    if corte is not None:
        # rectángulo gris semitransparente sobre la parte perdida
        mapa.add_vrect(x0=corte, x1=longitud, fillcolor="#888888", opacity=0.55, line_width=0,
                       annotation_text="se pierde", annotation_position="top right")
    mapa.add_vline(x=posicion, line_color="#c51b8a", line_width=3,
                   annotation_text=variante["cambio_p"], annotation_position="top left")
    mapa.update_layout(barmode="overlay", height=210, margin=dict(t=40, b=40),
                       xaxis_title="Posición (aminoácido)", xaxis_range=[0, longitud + 5],
                       yaxis=dict(autorange="reversed"))
    st.plotly_chart(mapa, key=clave)


# ---------------------------------------------------------------
# 7f. Dos pestañas: la variante y la proteína normal
# ---------------------------------------------------------------
es_truncante = variante["consecuencia"] in ("nonsense", "frameshift")
pestana_variante, pestana_normal = st.tabs(["🧬 Variante", "✅ Proteína normal"])

# ======================== PESTAÑA 1: VARIANTE ========================
with pestana_variante:
    vista = crear_vista()
    columna_visor, columna_info = st.columns([3, 2])

    if es_truncante:
        # ---------------- OPCIÓN 1: variante truncante ----------------
        # Todo desde la posición de la variante hasta el final se pierde (o se
        # altera, en un frameshift): lo pintamos como un "fantasma" gris.
        perdidos = list(range(posicion, longitud + 1))
        vista.setStyle({"resi": perdidos}, {"cartoon": {"color": "#e0e0e0", "opacity": 0.35}})
        marcar_variante(vista, {"colorscheme": "magentaCarbon"})
        vista.zoomTo()

        conservado = posicion - 1
        with columna_info:
            st.subheader("Variante truncante")
            st.write(f"**{variante['cambio_p']}** ({variante['consecuencia']}) · {variante['clasificacion']}")
            st.write(f"Se conservan los aminoácidos **1–{conservado}** "
                     f"({conservado / longitud:.0%} de la proteína). "
                     f"Desde el **{posicion}** hasta el **{longitud}** la proteína se pierde"
                     + (" o queda alterada." if variante["consecuencia"] == "frameshift" else "."))

            st.write("**Regiones afectadas:**")
            for _, fila in regiones_gen.iterrows():
                if fila["inicio"] >= posicion:
                    st.write(f"- ❌ {nombre_region(fila)} ({fila['inicio']}–{fila['fin']}): se pierde completa")
                elif fila["fin"] >= posicion:
                    st.write(f"- ⚠️ {nombre_region(fila)} ({fila['inicio']}–{fila['fin']}): se pierde en parte")

            st.caption("La parte gris transparente es lo que faltaría. Muchas veces la proteína truncada "
                       "ni siquiera se produce, porque la célula degrada el ARN defectuoso "
                       "(*nonsense-mediated decay*).")

    else:
        # ------- OPCIÓN 2: cambio puntual (missense, inframe, inicio, sinónima) -------
        # Vecinos: aminoácidos con algún átomo a menos de 5 Å de algún átomo del residuo
        vecinos = []
        for otro, sus_atomos in atomos.items():
            if otro == posicion:
                continue
            if any(distancia(a, b) < 5 for a in atomos[posicion] for b in sus_atomos):
                vecinos.append(otro)

        # Indicador aproximado de qué tan enterrado está: cuántos carbonos alfa
        # hay a menos de 10 Å de su carbono alfa
        cercanos = sum(1 for otro, ca in carbono_alfa.items()
                       if otro != posicion and distancia(ca, carbono_alfa[posicion]) < 10)
        if cercanos >= 18:
            entorno = "enterrado (dentro de la proteína)"
        elif cercanos >= 12:
            entorno = "parcialmente enterrado"
        else:
            entorno = "expuesto (en la superficie)"

        vista.addStyle({"resi": vecinos}, {"stick": {"colorscheme": "grayCarbon", "radius": 0.15}})
        marcar_variante(vista, {"colorscheme": "magentaCarbon", "radius": 0.3})
        vista.zoomTo({"resi": [posicion] + vecinos})

        nuevo = re.search(r"p\.[A-Za-z]{3}\d+([A-Z][a-z]{2})", str(variante["cambio_p"]))
        nuevo = nuevo.group(1).upper() if nuevo else None

        titulos = {"missense": "Variante missense", "inframe_indel": "Inserción/deleción sin corrimiento",
                   "perdida_inicio": "Pérdida del codón de inicio", "sinonima": "Variante sinónima"}
        with columna_info:
            st.subheader(titulos.get(variante["consecuencia"], "Variante"))
            st.write(f"**{variante['cambio_p']}** · {variante['clasificacion']}")
            st.write(f"- Original: {PROPIEDADES.get(original, original)}")
            if nuevo and nuevo != original:
                st.write(f"- Nuevo: {PROPIEDADES.get(nuevo, nuevo)}")
            st.write(f"- Región: {variante['regiones'] or 'sin región anotada'}")
            st.write(f"- Entorno: **{entorno}** ({cercanos} aminoácidos cercanos)")
            st.write(f"- Vecinos a menos de 5 Å: {len(vecinos)}")
            st.write(f"- Confianza del modelo (pLDDT): **{plddt[posicion]:.0f}**")

            if variante["consecuencia"] == "perdida_inicio":
                st.warning("Se pierde la metionina inicial: la proteína puede no producirse, "
                           "o empezar en otra metionina más adelante.")
            if variante["consecuencia"] == "sinonima":
                st.info("No cambia el aminoácido. En general no afecta la proteína, salvo que "
                        "altere el *splicing*.")
            if plddt[posicion] < 70:
                st.warning("La confianza del modelo en esta zona es baja: el entorno 3D no es confiable.")
            if original == "CYS":
                st.info("Es una cisteína: si forma un puente disulfuro, el cambio podría romperlo.")
            if variante["tirosina_ITAM"]:
                st.info("Es una tirosina del ITAM: sin ella no se puede fosforilar ni unir ZAP-70.")

            st.caption("El entorno se calcula sobre la cadena sola. En el complejo real, un residuo "
                       "de la superficie puede quedar enterrado en el contacto con otras cadenas.")

    with columna_visor:
        components.html(vista._make_html(), height=450)
        mostrar_leyenda(incluir_perdida=es_truncante)

    mapa_lineal("mapa_variante", corte=posicion if es_truncante else None)

# ===================== PESTAÑA 2: PROTEÍNA NORMAL =====================
with pestana_normal:
    vista_normal = crear_vista()
    # Marcamos dónde está la variante, para comparar con la otra pestaña
    marcar_variante(vista_normal, {"colorscheme": "magentaCarbon"})
    vista_normal.zoomTo()

    columna_visor, columna_info = st.columns([3, 2])
    with columna_visor:
        components.html(vista_normal._make_html(), height=450)
        mostrar_leyenda()
    with columna_info:
        st.subheader(f"{gen} sin la variante")
        st.write(f"Proteína completa: **{longitud}** aminoácidos.")
        st.write(f"En magenta se marca la posición **{posicion}**, donde está {variante['cambio_p']}, "
                 f"para comparar con la pestaña de la variante.")
        st.write("Regiones en esa posición: " + (variante["regiones"] or "ninguna anotada"))

    mapa_lineal("mapa_normal")
