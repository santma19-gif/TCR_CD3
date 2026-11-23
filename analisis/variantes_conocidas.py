"""
Variantes conocidas del complejo TCR-CD3 desde ClinVar.

Lee variant_summary.txt.gz de ClinVar, se queda con CD3D, CD3E, CD3G y CD247
en GRCh38, clasifica cada variante y la cruza con regiones.csv.

Genera:
  datos/ejemplos/variantes_conocidas.csv   tabla anotada
  datos/ejemplos/variantes_conocidas.vcf   VCF de prueba (reemplaza al de Bioinformática
                                           mientras el pipeline de Nanopore no esté listo)

Uso (desde la carpeta principal del proyecto):
  python analisis/variantes_conocidas.py

Antes hay que descargar ClinVar (no se sube a GitHub por su tamaño):
  wget -P datos/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
"""

import re
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Rutas y parámetros
# ---------------------------------------------------------------------------
CLINVAR = Path("datos/clinvar/variant_summary.txt.gz")
REGIONES = Path("datos/regiones_funcionales/regiones.csv")
REFERENCIAS = Path("datos/referencias")
ESTRUCTURAS = Path("datos/estructuras")
SALIDA = Path("datos/ejemplos")

GENES = ["CD3D", "CD3E", "CD3G", "CD247"]
ENSAMBLAJE = "GRCh38"

# Solo estas columnas de ClinVar (el archivo completo tiene ~40)
COLUMNAS = [
    "VariationID", "GeneSymbol", "Type", "Name", "ClinicalSignificance",
    "ReviewStatus", "PhenotypeList", "Assembly", "Chromosome",
    "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
]

# Código de 3 letras → 1 letra (para comparar con el FASTA)
AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}


# ---------------------------------------------------------------------------
# 1. Leer y filtrar ClinVar
# ---------------------------------------------------------------------------
def leer_clinvar(ruta):
    """Lee ClinVar por pedazos y se queda solo con los 4 genes en GRCh38."""
    pedazos = []
    for pedazo in pd.read_csv(ruta, sep="\t", usecols=COLUMNAS, dtype=str,
                              chunksize=200_000):
        filtro = (pedazo["Assembly"] == ENSAMBLAJE) & pedazo["GeneSymbol"].isin(GENES)
        pedazos.append(pedazo[filtro])
    return pd.concat(pedazos, ignore_index=True)


# ---------------------------------------------------------------------------
# 2. Separar el nombre de la variante
#    Ejemplo: "NM_000733.4(CD3E):c.520C>T (p.Arg174Ter)"
# ---------------------------------------------------------------------------
PATRON_NOMBRE = re.compile(
    r"^(?P<transcrito>[NX]M_\d+\.\d+)\([^)]*\):(?P<hgvs_c>c\.[^ ]+)"
    r"(?: \((?P<hgvs_p>p\.[^)]+)\))?"
)
PATRON_PROTEINA = re.compile(r"^p\.(?P<aa_ref>[A-Z][a-z]{2})(?P<posicion>\d+)")


def separar_nombre(nombre):
    """Devuelve transcrito, cambio en el ADN (c.) y cambio en la proteína (p.)."""
    m = PATRON_NOMBRE.match(str(nombre))
    if not m:
        return pd.Series({"transcrito": None, "hgvs_c": None, "hgvs_p": None})
    return pd.Series(m.groupdict())


def posicion_proteina(hgvs_p):
    """De 'p.Tyr188Cys' saca ('Tyr', 188). Si no aplica, (None, None)."""
    if not isinstance(hgvs_p, str):
        return pd.Series({"aa_ref": None, "posicion_proteina": None})
    m = PATRON_PROTEINA.match(hgvs_p)
    if not m:
        return pd.Series({"aa_ref": None, "posicion_proteina": None})
    return pd.Series({"aa_ref": m["aa_ref"], "posicion_proteina": int(m["posicion"])})


# ---------------------------------------------------------------------------
# 3. Clasificar
# ---------------------------------------------------------------------------
def grupo_clinico(significado):
    """Agrupa la clasificación de ClinVar en 5 grupos (mismo criterio que Bioinformática)."""
    s = str(significado)
    if "onflicting" in s:
        return "Conflicto"
    if "enign" in s:
        return "Benigna"
    if "athogenic" in s:
        return "Patogenica"
    if "ncertain" in s:
        return "VUS"
    return "Sin_clasificar"


def consecuencia(hgvs_c, hgvs_p):
    """Tipo de cambio, leyendo la nomenclatura HGVS."""
    p = hgvs_p if isinstance(hgvs_p, str) else ""
    c = hgvs_c if isinstance(hgvs_c, str) else ""

    if p:
        if p.endswith("=") or "=" in p:
            return "sinonima"
        if "fs" in p:
            return "frameshift"
        if "ext" in p:
            return "perdida_stop"  # ej. p.Ter208GlnextTer5: se pierde el codón de parada
        if "Ter" in p or p.endswith("*"):
            return "nonsense"
        if p.startswith("p.Met1"):
            return "perdida_inicio"
        if re.search(r"del|ins|dup", p):
            return "inframe_indel"
        return "missense"

    # Sin cambio en la proteína: mirar el cambio en el ADN
    if re.match(r"c\.-", c):
        return "UTR_5"
    if re.match(r"c\.\*", c):
        return "UTR_3"
    if re.search(r"[+-][12](?!\d)", c):
        return "splicing"          # ±1 o ±2 del borde del exón
    if re.search(r"\d[+-]\d", c):
        return "intronica"
    return "otra"


# ---------------------------------------------------------------------------
# 4. Cruzar con regiones.csv y con el FASTA
# ---------------------------------------------------------------------------
def secuencia_desde_pdb(ruta):
    """Lee la secuencia del modelo de AlphaFold (un aminoácido por carbono alfa)."""
    tres_a_una = {k.upper(): v for k, v in AA3.items()}
    residuos = {}
    for linea in ruta.read_text().splitlines():
        if linea.startswith("ATOM") and linea[12:16].strip() == "CA":
            residuos[int(linea[22:26])] = tres_a_una.get(linea[17:20], "X")
    return "".join(residuos[i] for i in sorted(residuos))


def leer_secuencias(regiones):
    """
    Secuencia de cada proteína → {gen: secuencia}.
    Usa el FASTA de datos/referencias/; si falta, la saca del modelo de AlphaFold
    (misma numeración, ya verificada). Devuelve también de dónde salió cada una.
    """
    secuencias, origen = {}, {}
    for gen in GENES:
        fasta = REFERENCIAS / f"{gen}_proteina.fasta"
        if fasta.exists():
            lineas = fasta.read_text().splitlines()
            secuencias[gen] = "".join(l.strip() for l in lineas if not l.startswith(">"))
            origen[gen] = "FASTA"
            continue
        uniprot = regiones.loc[regiones["gen"] == gen, "uniprot"]
        if uniprot.empty:
            continue
        modelos = sorted(ESTRUCTURAS.glob(f"AF-{uniprot.iloc[0]}-F1-model_v*.pdb"))
        if modelos:
            secuencias[gen] = secuencia_desde_pdb(modelos[-1])
            origen[gen] = "AlphaFold"
    return secuencias, origen


def tirosinas_canonicas_itam(regiones, secuencias):
    """
    Tirosinas de los ITAM que siguen el patrón YxxL/I → {gen: {posiciones}}.
    Un ITAM es YxxL/I – (6 a 8 aa) – YxxL/I: solo esas tirosinas se fosforilan.
    Otras tirosinas que caen dentro del rango del ITAM (como la Y64 de CD3ζ) no cuentan.
    """
    resultado = {}
    for _, fila in regiones[regiones["descripcion"] == "ITAM"].iterrows():
        sec = secuencias.get(fila["gen"])
        if sec is None:
            continue
        for pos in range(fila["inicio"], fila["fin"] + 1):
            # pos está numerada desde 1; en Python la cadena empieza en 0
            if sec[pos - 1] == "Y" and pos + 2 < len(sec) and sec[pos + 2] in "LI":
                resultado.setdefault(fila["gen"], set()).add(pos)
    return resultado


def anotar_regiones(fila, regiones):
    """Lista de regiones del CSV en las que cae la posición (pueden ser varias)."""
    pos = fila["posicion_proteina"]
    if pd.isna(pos):
        return ""
    r = regiones[(regiones["gen"] == fila["GeneSymbol"]) &
                 (regiones["inicio"] <= pos) & (regiones["fin"] >= pos)]
    nombres = [d if isinstance(d, str) and d else t
               for t, d in zip(r["tipo"], r["descripcion"])]
    return ";".join(nombres)


def tirosina_itam(fila, tirosinas):
    """True si la variante cambia una de las tirosinas canónicas (YxxL/I) de un ITAM."""
    pos = fila["posicion_proteina"]
    if pd.isna(pos) or fila["aa_ref"] != "Tyr":
        return False
    if str(fila["hgvs_p"]).endswith("="):
        return False  # sinónima: la tirosina no cambia
    return int(pos) in tirosinas.get(fila["GeneSymbol"], set())


def revisar_fasta(fila, secuencias):
    """
    Compara el aminoácido de referencia de ClinVar con el FASTA canónico.
    Si no coinciden, el transcrito de ClinVar no corresponde a la isoforma
    canónica (o la numeración está corrida).
    """
    sec = secuencias.get(fila["GeneSymbol"])
    pos = fila["posicion_proteina"]
    if sec is None or pd.isna(pos) or fila["aa_ref"] not in AA3:
        return "sin_revisar"
    pos = int(pos)
    if fila["aa_ref"] == "Ter" and pos == len(sec) + 1:
        return "ok"  # el codón de parada va justo después del último aminoácido
    if pos > len(sec):
        return "fuera_de_rango"
    return "ok" if sec[pos - 1] == AA3[fila["aa_ref"]] else "NO_COINCIDE"


# ---------------------------------------------------------------------------
# 5. Escribir el VCF de prueba
# ---------------------------------------------------------------------------
def escribir_vcf(tabla, ruta):
    """VCF mínimo (formato 4.2) con las variantes que tienen coordenadas."""
    v = tabla.dropna(subset=["PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"])
    v = v[(v["PositionVCF"] != "-1") & (v["ReferenceAlleleVCF"] != "na")
          & (v["AlternateAlleleVCF"] != "na")].copy()
    v["pos"] = v["PositionVCF"].astype(int)
    v["chrom"] = "chr" + v["Chromosome"]
    v = v.sort_values(["chrom", "pos"])

    with open(ruta, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##reference={ENSAMBLAJE}\n")
        f.write("##source=ClinVar variant_summary (variantes_conocidas.py)\n")
        for chrom in sorted(v["chrom"].unique()):
            f.write(f"##contig=<ID={chrom}>\n")
        f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gen">\n')
        f.write('##INFO=<ID=GRUPO,Number=1,Type=String,Description="Grupo clinico (ClinVar)">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for _, x in v.iterrows():
            f.write(f"{x['chrom']}\t{x['pos']}\tclinvar{x['VariationID']}\t"
                    f"{x['ReferenceAlleleVCF']}\t{x['AlternateAlleleVCF']}\t.\tPASS\t"
                    f"GENE={x['GeneSymbol']};GRUPO={x['grupo']}\n")
    return len(v)


# ---------------------------------------------------------------------------
# Programa principal
# ---------------------------------------------------------------------------
def main():
    if not CLINVAR.exists():
        raise SystemExit(f"No se encontró {CLINVAR}. Descárgalo primero (ver el inicio del script).")

    tabla = leer_clinvar(CLINVAR)
    print(f"Variantes en {', '.join(GENES)} ({ENSAMBLAJE}): {len(tabla)}")

    tabla = pd.concat([tabla, tabla["Name"].apply(separar_nombre)], axis=1)
    tabla = pd.concat([tabla, tabla["hgvs_p"].apply(posicion_proteina)], axis=1)
    tabla["posicion_proteina"] = tabla["posicion_proteina"].astype("Int64")
    tabla["grupo"] = tabla["ClinicalSignificance"].apply(grupo_clinico)
    tabla["consecuencia"] = [consecuencia(c, p) for c, p in zip(tabla["hgvs_c"], tabla["hgvs_p"])]
    # Variantes que ClinVar describe solo en coordenadas genómicas
    # (ej. "NC_000011.10:g.118343086G>A"): no hay cambio en el transcrito para leer
    solo_g = tabla["Name"].str.startswith("NC_", na=False) & tabla["hgvs_c"].isna()
    tabla.loc[solo_g, "consecuencia"] = "solo_genomica"

    regiones = pd.read_csv(REGIONES)
    secuencias, origen = leer_secuencias(regiones)
    tirosinas = tirosinas_canonicas_itam(regiones, secuencias)
    tabla["regiones"] = tabla.apply(anotar_regiones, axis=1, regiones=regiones)
    tabla["tirosina_ITAM"] = tabla.apply(tirosina_itam, axis=1, tirosinas=tirosinas)
    tabla["revision_fasta"] = tabla.apply(revisar_fasta, axis=1, secuencias=secuencias)

    columnas_salida = [
        "GeneSymbol", "VariationID", "Name", "transcrito", "hgvs_c", "hgvs_p",
        "aa_ref", "posicion_proteina", "consecuencia", "grupo",
        "ClinicalSignificance", "ReviewStatus", "regiones", "tirosina_ITAM",
        "revision_fasta", "PhenotypeList", "Chromosome", "PositionVCF",
        "ReferenceAlleleVCF", "AlternateAlleleVCF",
    ]
    SALIDA.mkdir(parents=True, exist_ok=True)
    tabla[columnas_salida].rename(columns={"GeneSymbol": "gen", "Name": "nombre_clinvar"}).to_csv(
        SALIDA / "variantes_conocidas.csv", index=False)
    n_vcf = escribir_vcf(tabla, SALIDA / "variantes_conocidas.vcf")

    # ---- Resumen en la terminal ----
    print("\nGrupo clínico por gen:")
    print(pd.crosstab(tabla["GeneSymbol"], tabla["grupo"], margins=True, margins_name="Total"))
    print("\nConsecuencia por grupo clínico:")
    print(pd.crosstab(tabla["consecuencia"], tabla["grupo"], margins=True, margins_name="Total"))
    print("\nTirosinas canónicas de los ITAM (YxxL/I):")
    for gen in GENES:
        posiciones = ", ".join(f"Y{p}" for p in sorted(tirosinas.get(gen, [])))
        print(f"  {gen}: {posiciones or 'sin secuencia'}  (secuencia de {origen.get(gen, '—')})")
    print(f"Variantes que cambian una de esas tirosinas: {tabla['tirosina_ITAM'].sum()}")
    itam = tabla[tabla["tirosina_ITAM"]]
    for _, x in itam.iterrows():
        print(f"  {x['GeneSymbol']} {x['hgvs_p']} · {x['grupo']}")
    print("\nRevisión contra el FASTA canónico:")
    print(tabla["revision_fasta"].value_counts().to_string())
    transcritos = tabla.dropna(subset=["transcrito"]).groupby("GeneSymbol")["transcrito"].unique()
    print("\nTranscritos usados por ClinVar:")
    for gen, ts in transcritos.items():
        print(f"  {gen}: {', '.join(ts)}")
    print(f"\nGuardado: {SALIDA / 'variantes_conocidas.csv'}")
    print(f"Guardado: {SALIDA / 'variantes_conocidas.vcf'} ({n_vcf} variantes)")


if __name__ == "__main__":
    main()
