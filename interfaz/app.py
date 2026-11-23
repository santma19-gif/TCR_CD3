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