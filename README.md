# Variantes del complejo TCR-CD3

Herramienta para identificar y caracterizar variantes en los genes **CD3D, CD3E, CD3G y CD247** y relacionarlas con patologías inmunológicas reportadas.

Proyecto de la materia *Fundamentos de Programación para Ciencias Biológicas*.
Autores: Santiago Matta Muñoz, María José Oquendo Chavarría.

> Herramienta con fines informativos y de apoyo al análisis bioinformático. No sustituye una evaluación ni un diagnóstico clínico.

## Flujo general

```
Lecturas Nanopore (FASTQ)
  → control de calidad → mapeo a GRCh38 → extracción de los 4 genes
  → llamado de variantes (VCF)
  → anotación: cambio de proteína, región funcional, ITAM, ClinVar
  → reporte
```

## Estructura

| Carpeta | Contenido |
|---|---|
| `datos/referencias/` | Secuencias de proteína de UniProt (FASTA) |
| `datos/regiones_funcionales/` | `regiones.csv` (curado a mano desde UniProt) |
| `datos/estructuras/` | Modelos de AlphaFold (`.pdb`) |
| `datos/ejemplos/` | Variantes conocidas y VCF de prueba |
| `datos/clinvar/` | `variant_summary.txt.gz` (no se sube; ver "Datos") |
| `analisis/` | Código de análisis: lectura del VCF, anotación, reporte |
| `interfaz/` | Aplicación en Streamlit |
| `pruebas/` | Pruebas del código |
| `documentos/` | Bitácora y anteproyecto |
| `resultados/` | Archivos generados por el código (no se suben) |

Convenciones:
- Los datos descargados no se editan a mano. Lo que se filtra o limpia se genera con código en `resultados/`.
- Se usa la isoforma canónica de UniProt y la numeración desde Met1.
- Coordenadas genómicas en GRCh38.

## Instalación

```bash
conda env create -f environment.yml
conda activate tcrcd3
```

## Datos

ClinVar (no se incluye en el repositorio por su tamaño):

```bash
wget -P datos/clinvar https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
```

## Uso

Desde la carpeta principal del proyecto:

```bash
streamlit run interfaz/app.py
streamlit run interfaz/ver_proteinas.py
streamlit run interfaz/ver_estructura.py
```

## Referencias

- Jumper J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*.
- Richards S. et al. (2015). Standards and guidelines for the interpretation of sequence variants. *Genet Med* 17(5):405-24.
- de Saint Basile G. et al. (2004). Severe combined immunodeficiency caused by deficiency in either the delta or the epsilon subunit of CD3. *J Clin Invest* 114(10):1512-7.
