from pathlib import Path

import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

CARPETA_ESTRUCTURAS = Path("datos/estructuras")
ARCHIVO_REGIONES = "datos/regiones_funcionales/regiones.csv"

# Nombre de cada proteína para mostrar en el título
NOMBRES = {
    "CD3E": "CD3ε",
    "CD3D": "CD3δ",
    "CD3G": "CD3γ",
    "CD247": "CD3ζ",
}

# ---------------------------------------------------------------
# 1. Leer el CSV y elegir la proteína
# ---------------------------------------------------------------
st.title("Estructura 3D del complejo CD3 (AlphaFold)")

regiones = pd.read_csv(ARCHIVO_REGIONES)

gen = st.selectbox("Proteína", regiones["gen"].unique())
regiones_gen = regiones[regiones["gen"] == gen]
codigo_uniprot = regiones_gen["uniprot"].iloc[0]

st.subheader(f"{NOMBRES.get(gen, gen)} · gen {gen} · UniProt {codigo_uniprot}")

# ---------------------------------------------------------------
# 2. Buscar el archivo PDB de esa proteína
# ---------------------------------------------------------------
# Los archivos de AlphaFold se llaman AF-<código UniProt>-F1-model_v<versión>.pdb
# Buscamos con un comodín (*) en la versión, por si AlphaFold la actualiza.
encontrados = sorted(CARPETA_ESTRUCTURAS.glob(f"AF-{codigo_uniprot}-F1-model_v*.pdb"))

if not encontrados:
    st.error(f"No se encontró el modelo de {gen}. Debe haber un archivo "
             f"`AF-{codigo_uniprot}-F1-model_vX.pdb` en `{CARPETA_ESTRUCTURAS}/`.")
    st.stop()   # detiene la página aquí

archivo_pdb = encontrados[-1]   # si hay varias versiones, usa la más reciente
st.caption(f"Archivo: `{archivo_pdb.name}`")

pdb = archivo_pdb.read_text()

# ---------------------------------------------------------------
# 3. Sacar del PDB el pLDDT y el nombre de cada aminoácido
# ---------------------------------------------------------------
# En el formato PDB cada dato ocupa columnas fijas:
#   columnas 13-16 -> nombre del átomo     -> linea[12:16]
#   columnas 18-20 -> nombre del residuo   -> linea[17:20]
#   columnas 23-26 -> número del residuo   -> linea[22:26]
#   columnas 61-66 -> B-factor (= pLDDT)   -> linea[60:66]
plddt = {}
aminoacidos = {}
for linea in pdb.splitlines():
    if linea.startswith("ATOM") and linea[12:16].strip() == "CA":
        residuo = int(linea[22:26])
        plddt[residuo] = float(linea[60:66])
        aminoacidos[residuo] = linea[17:20]

# Comprobación: el modelo debe medir lo mismo que la proteína en el CSV
longitud_csv = regiones_gen["fin"].max()
if len(plddt) != longitud_csv:
    st.warning(f"El modelo tiene {len(plddt)} aminoácidos, pero según el CSV "
               f"{gen} mide {longitud_csv}. Revisa que sea la isoforma canónica.")

# ---------------------------------------------------------------
# 4. Colores
# ---------------------------------------------------------------
def color_confianza(valor):
    """Los mismos 4 colores que usa la página de AlphaFold."""
    if valor > 90:
        return "#0053D6"   # muy alta
    elif valor > 70:
        return "#65CBF3"   # alta
    elif valor > 50:
        return "#FFDB13"   # baja
    else:
        return "#FF7D45"   # muy baja


# Los mismos colores de ver_proteinas.py, para que ambos dibujos coincidan
COLORES = {
    "peptido_senal": "#9e9e9e",
    "extracelular": "#4c9be8",
    "transmembrana": "#f2b134",
    "citoplasmatico": "#7cc47f",
    "Ig-like": "#1f5fa8",
    "ITAM": "#d64545",
}


def nombre_region(fila):
    """El nombre de cada fila del CSV (igual que en ver_proteinas.py)."""
    if fila["tipo"] in ("peptido_senal", "transmembrana"):
        return fila["tipo"]
    return fila["descripcion"]


# Orden en que se pintan las regiones. Como se superponen, las que se
# pintan después quedan "encima": primero la topología y luego los
# dominios (así el ITAM se ve sobre la zona citoplasmática).
ORDEN_DE_PINTADO = ["peptido_senal", "dominio_topologico", "transmembrana", "dominio"]

# ---------------------------------------------------------------
# 5. Opciones para el usuario
# ---------------------------------------------------------------
modo = st.radio("Colorear por", ["Confianza (pLDDT)", "Regiones del CSV"], horizontal=True)
mostrar_tirosinas = st.checkbox("Resaltar las tirosinas del ITAM")

# ---------------------------------------------------------------
# 6. Crear el visor 3D
# ---------------------------------------------------------------
vista = py3Dmol.view(width=700, height=500)
vista.addModel(pdb, "pdb")

if modo == "Confianza (pLDDT)":
    for residuo, valor in plddt.items():
        vista.setStyle({"resi": residuo}, {"cartoon": {"color": color_confianza(valor)}})
else:
    # Todo gris primero, por si queda alguna zona sin anotar
    vista.setStyle({}, {"cartoon": {"color": "#bdbdbd"}})
    for tipo in ORDEN_DE_PINTADO:
        for _, fila in regiones_gen[regiones_gen["tipo"] == tipo].iterrows():
            posiciones = list(range(fila["inicio"], fila["fin"] + 1))
            color = COLORES.get(nombre_region(fila), "#bdbdbd")
            vista.setStyle({"resi": posiciones}, {"cartoon": {"color": color}})

# Buscar las tirosinas canónicas de cada ITAM: las del patrón YxxL/I
# (tirosina, dos aminoácidos cualquiera y luego leucina o isoleucina).
# Otras tirosinas que caen dentro del rango del ITAM no se fosforilan.
tirosinas_itam = []
for _, fila in regiones_gen[regiones_gen["descripcion"] == "ITAM"].iterrows():
    for residuo in range(fila["inicio"], fila["fin"] + 1):
        if aminoacidos.get(residuo) == "TYR" and aminoacidos.get(residuo + 3) in ("LEU", "ILE"):
            tirosinas_itam.append(residuo)

if mostrar_tirosinas:
    # addStyle agrega sin borrar la cinta; "stick" dibuja los átomos como palitos
    vista.addStyle({"resi": tirosinas_itam}, {"stick": {"colorscheme": "magentaCarbon"}})
    for residuo in tirosinas_itam:
        vista.addLabel(f"Y{residuo}",
                       {"fontSize": 14, "fontColor": "black", "backgroundColor": "white",
                        "backgroundOpacity": 0.8},
                       {"resi": residuo, "atom": "OH"})

vista.setBackgroundColor("white")
vista.zoomTo()

# Mostrar el visor dentro de la página de Streamlit
components.html(vista._make_html(), height=520)

# ---------------------------------------------------------------
# 7. Leyenda
# ---------------------------------------------------------------
if modo == "Confianza (pLDDT)":
    st.caption("Azul oscuro: muy alta (>90) · Azul claro: alta (70–90) · "
               "Amarillo: baja (50–70) · Naranja: muy baja (<50)")
else:
    partes = [f"<span style='color:{color}; font-size:20px'>■</span> {nombre}"
              for nombre, color in COLORES.items()]
    st.markdown(" &nbsp; ".join(partes), unsafe_allow_html=True)

if mostrar_tirosinas:
    st.write("Tirosinas del ITAM encontradas:", ", ".join(f"Y{r}" for r in tirosinas_itam))

# ---------------------------------------------------------------
# 8. Gráfico del pLDDT a lo largo de la proteína
# ---------------------------------------------------------------
st.subheader("Confianza (pLDDT) por posición")
st.line_chart(pd.Series(plddt, name="pLDDT"), x_label="Posición (aminoácido)", y_label="pLDDT")