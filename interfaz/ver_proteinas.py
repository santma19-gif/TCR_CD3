import pandas as pd
import plotly.graph_objects as go
import streamlit as st 

st.title("Regiones funcionales del complejo CD3")

# Leer el CSV 
regiones = pd.read_csv("datos/regiones_funcionales/regiones.csv")

# Elegir la proteína que queremos ver 
gen = st.selectbox("Proteína", regiones["gen"].unique())
datos_gen = regiones[regiones["gen"]== gen ]

# Configuración del dibujo 
CARRILES = {
    "peptido_senal": "Topología",
    "dominio_topologico": "Topología",
    "transmembrana": "Topología",
    "dominio": "Dominios",
    "motif": "Motivos",
}
 
# Color de cada región; las que no estén aquí salen en gris
COLORES = {
    "peptido_senal": "#9e9e9e",
    "extracelular": "#4c9be8",
    "transmembrana": "#f2b134",
    "citoplasmatico": "#7cc47f",
    "Ig-like": "#1f5fa8",
    "ITAM": "#d64545",
}

def carril(fila):
    """En qué fila del dibujo va cada región."""
    if fila["tipo"] == "region":
        return fila["descripcion"]
    return CARRILES.get(fila["tipo"], "Otros")
 
 
def nombre_region(fila):
    """El nombre que se muestra en el dibujo para cada fila del CSV."""
    if fila["tipo"] in ("peptido_senal", "transmembrana"):
        return fila["tipo"]
    return fila["descripcion"]

# Dibujar una barra horizontal por cada fila del CSV

fig = go.Figure()
 
for _, fila in datos_gen.iterrows():
    nombre = nombre_region(fila)
    fig.add_trace(go.Bar(
        y=[carril(fila)],                          # en qué carril va
        base=[fila["inicio"]],                     # dónde empieza
        x=[fila["fin"] - fila["inicio"]],          # cuánto mide
        orientation="h",
        name=nombre,
        text=nombre,
        textposition="inside",
        insidetextanchor="middle",
        marker_color=COLORES.get(nombre, "#bdbdbd"),
        hovertemplate=f"<b>{nombre}</b><br>Posición {fila['inicio']}–{fila['fin']}<extra></extra>",
    ))
 
# Orden de los carriles, de abajo hacia arriba
otras_regiones = list(datos_gen.loc[datos_gen["tipo"] == "region", "descripcion"])
orden_carriles = ["Topología", "Dominios"] + otras_regiones
 
fig.update_layout(
    barmode="overlay",
    showlegend=False,
    height=120 + 45 * len(orden_carriles),
    margin=dict(t=20),
    xaxis_title="Posición (aminoácido)",
    xaxis_range=[0, datos_gen["fin"].max() + 5],
    yaxis=dict(categoryorder="array", categoryarray=orden_carriles, autorange="reversed"),
)
 
st.plotly_chart(fig)

# 5. Mostrar la tabla debajo, para comparar con el dibujo

st.subheader("Datos del CSV")
st.dataframe(datos_gen)