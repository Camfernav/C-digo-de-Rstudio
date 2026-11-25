# C-digo-de-Rstudio
Código de R-studio para proyecto de genética de Arabidopsis thaliana
Bases de datos obtenidas de: 
1- Glucosinolatos de semilla: Katz et al. (2021) 
 - Artículo: https://pmc.ncbi.nlm.nih.gov/articles/PMC8205490
 - Supplementary file 1 (elife-67784-supp1.xlsx, 547.7 KB).

2- Lista de accesiones: 1001 Genomes
 - Enlace: https://1001genomes.org/1135/accessions
 - Seleccionamos Download CSV y el archivo trae ID, nombre, número CS, latitud, longitud y ancestría de las 1.135 accesiones.

3- Variantes de MYB28: VCF Subset del 1001 Genomes
 - Enlace: https://tools.1001genomes.org/vcfsubset
 - En la sección de "Select Strains" seleccionar todas proseguir con Next, en la parte de "Select Loci" en Regions poner Chr5:24689069..24693187. Este es el gen AT5G61420 completo más 2 kb de promotor (permite ver todos los lugares donde una variante podría cambiar lo que hace MYB28).
 - En "Select Formats VCF Source" seleccionar Full genome SnpEFF VCF (este trae el tipo de cada variante)
 - En "File format" seleccionamos VCF y después Download now
