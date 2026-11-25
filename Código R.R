


# =============================================================================
# SECCIÓN 0. PREPARACIÓN: carpeta de trabajo, paquetes y archivos
# =============================================================================

setwd("C:/Users/Camfe/Downloads")          # carpeta que CONTIENE a "datos" (en R se usa / y no \)
cat("Carpeta de trabajo:", getwd(), "\n")  # muestra en qué carpeta está trabajando R

# Instala los paquetes solo si no están instalados todavía
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")  # para leer Excel
if (!requireNamespace("dplyr",  quietly = TRUE)) install.packages("dplyr")   # para manejar tablas

library(readxl)   # carga el paquete que lee archivos .xlsx
library(dplyr)    # carga el paquete para filtrar, unir y resumir tablas

dir.create("resultados", showWarnings = FALSE)   # crea la carpeta "resultados" (si ya existe, no hace nada)

print(list.files("datos"))   # muestra los archivos que hay dentro de "datos"

# Busca cada archivo por su extensión dentro de la carpeta "datos"
archivo_katz <- list.files("datos", pattern = "\\.xlsx$", full.names = TRUE)[1]  # primer Excel
archivo_acc  <- list.files("datos", pattern = "\\.csv$",  full.names = TRUE)[1]  # primer CSV
vcfs         <- list.files("datos", pattern = "\\.vcf$",  full.names = TRUE)     # todos los VCF
archivo_vcf  <- vcfs[grepl("snpeff", vcfs, ignore.case = TRUE)][1]              # el que dice "snpeff"

print(c(katz = archivo_katz, accesiones = archivo_acc, vcf = archivo_vcf))  # muestra qué archivo tomó

# Si falta alguno, detiene el script con un mensaje claro
if (is.na(archivo_katz)) stop("No hay ningún .xlsx en la carpeta 'datos'")
if (is.na(archivo_acc))  stop("No hay ningún .csv en la carpeta 'datos'")
if (is.na(archivo_vcf))  stop("No hay un .vcf con 'snpeff' en el nombre en 'datos'")

cat("✔ SECCIÓN 0 OK: archivos encontrados\n\n")


# =============================================================================
# SECCIÓN 1. FENOTIPO: glucosinolatos (GSL) de semilla, Katz et al. (2021)
# =============================================================================

print(excel_sheets(archivo_katz))   # muestra los nombres de las hojas del Excel

katz <- read_excel(archivo_katz, sheet = "D. GSLs_data emmeans")  # hoja D: una fila por accesión

# Lista de los 20 GSL alifáticos (según la hoja A); se excluyen los 3 indólicos
gsl_alifaticos <- c("3OHP", "3MSO", "4OHB", "2-OH-3-Butenyl", "4MSO", "Allyl",
                    "5MSO", "3-Butenyl", "Branched", "OHPentenyl", "6MSO",
                    "3MT", "7MSO", "4Pentenyl", "4MT", "8MSO", "7MT", "BZO",
                    "8MT", "8MTder")

gsl_c4 <- c("4OHB", "2-OH-3-Butenyl", "4MSO", "3-Butenyl", "4MT")  # GSL de 4 carbonos

recortar <- function(x) pmax(x, 0)   # función: cambia valores negativos (ruido) por 0

katz <- katz %>%                                                     # toma la tabla katz y...
  mutate(
    gsl_alif = rowSums(across(all_of(gsl_alifaticos), recortar)) * 1000,  # suma de alifáticos, en nmol/mg
    gsl_4mso = recortar(`4MSO`) * 1000,                                    # glucorafanina, en nmol/mg
    frac_c4  = rowSums(across(all_of(gsl_c4), recortar)) * 1000 / gsl_alif, # proporción de GSL de 4 carbonos
    MAM      = factor(`MAM status`),                                       # estado del locus MAM (C3 o C4)
    AOP      = factor(`AOP status`)                                        # estado del locus AOP (Alk, OH, MSO)
  )

cat("Accesiones con fenotipo:", nrow(katz), "(esperado: 797)\n")    # cuántas filas tiene la tabla
cat("✔ SECCIÓN 1 OK\n\n")


# =============================================================================
# SECCIÓN 2. LISTA DE ACCESIONES DEL 1001 GENOMES
# =============================================================================

acc <- read.csv(archivo_acc, colClasses = "character")   # lee el CSV; todo como texto

# La columna CS_Number trae un enlace ".../CS76790"; se extrae solo "CS76790"
acc$CS <- regmatches(acc$CS_Number, regexpr("CS[0-9]+", acc$CS_Number))

acc <- acc %>% select(ID, CS, Admixture_Group)   # se quedan solo las 3 columnas necesarias

cat("Accesiones del 1001 Genomes:", nrow(acc), "(esperado: 1135)\n")
cat("✔ SECCIÓN 2 OK\n\n")


# =============================================================================
# SECCIÓN 3. EMPAREJAMIENTO GENOTIPO–FENOTIPO Y GRUPOS EUROPEOS
# =============================================================================

muestras <- inner_join(katz, acc, by = "CS")   # une las dos tablas; solo quedan accesiones en AMBAS
cat("Accesiones emparejadas:", nrow(muestras), "(esperado: 791)\n")

# Países de cada grupo geográfico (códigos que usa el archivo de Katz)
norte  <- c("SWE", "FIN", "NOR", "DEN", "EST", "LTU")
centro <- c("GER", "UK", "FRA", "CZE", "NED", "AUT", "SVK", "SUI", "POL", "BEL")
sur    <- c("ESP", "ITA", "POR", "GRC", "BUL", "SRB", "ROU", "CRO")

muestras <- muestras %>%
  mutate(grupo = case_when(Country %in% norte  ~ "Norte",     # si el país está en norte -> "Norte"
                           Country %in% centro ~ "Centro",    # si está en centro -> "Centro"
                           Country %in% sur    ~ "Sur",       # si está en sur -> "Sur"
                           TRUE                ~ "Excluida")) # cualquier otro país -> "Excluida"

print(table(muestras$Country[muestras$grupo == "Excluida"]))  # muestra qué países se excluyen (62)

muestras <- muestras %>%
  filter(grupo != "Excluida") %>%                                     # quita las excluidas
  mutate(grupo = factor(grupo, levels = c("Norte", "Centro", "Sur")),  # orden fijo de los grupos
         Admixture_Group = factor(Admixture_Group))                    # grupo genético como factor

print(table(muestras$grupo))   # esperado: Norte 171, Centro 278, Sur 280
cat("Accesiones europeas:", nrow(muestras), "(esperado: 729)\n")
cat("✔ SECCIÓN 3 OK\n\n")


# =============================================================================
# SECCIÓN 4. LEER EL VCF Y CONVERTIR GENOTIPOS A 0 / 1 / 2
# =============================================================================

lineas <- readLines(archivo_vcf)   # lee el VCF completo como líneas de texto

encabezado <- strsplit(lineas[grep("^#CHROM", lineas)], "\t")[[1]]  # nombres de columnas (línea #CHROM)

vcf <- read.table(text = lineas[!startsWith(lineas, "#")],  # solo las líneas de datos (sin "#")
                  sep = "\t",                                # separadas por tabulador
                  colClasses = "character",                  # todo como texto
                  col.names = encabezado,                    # con los nombres del encabezado
                  check.names = FALSE, quote = "", comment.char = "")

colnames(vcf)[1] <- "CHROM"            # renombra "#CHROM" a "CHROM"
vcf$POS <- as.integer(vcf$POS)         # la posición se convierte a número

vcf <- vcf[nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1, ]   # deja solo SNV bialélicas (quita indels)
cat("SNV bialélicas:", nrow(vcf), "(esperado: 199)\n")

if (!all(muestras$ID %in% colnames(vcf))) stop("Hay accesiones sin columna en el VCF")  # verificación

gt <- as.matrix(vcf[, muestras$ID])       # matriz de genotipos: filas = SNV, columnas = accesiones
gt <- sub(":.*", "", gt)                  # de "0|0:40:12" deja solo "0|0"
gt <- gsub("|", "/", gt, fixed = TRUE)    # cambia "0|0" por "0/0"

codigo <- c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)   # cuántas copias del alelo ALT

dosis <- matrix(codigo[gt],                               # traduce cada genotipo a 0, 1 o 2 (o NA)
                nrow = nrow(gt),
                dimnames = list(vcf$POS, muestras$ID))    # filas = posición, columnas = accesión

cat("Matriz de dosis:", nrow(dosis), "SNV x", ncol(dosis), "accesiones (esperado: 199 x 729)\n")
cat("✔ SECCIÓN 4 OK\n\n")


# =============================================================================
# SECCIÓN 5. ANOTACIÓN (SnpEff) Y FILTROS DE CALIDAD
# =============================================================================

# Función: del campo INFO extrae el efecto sobre el transcrito AT5G61420.2 (proteína de 366 aa)
efecto_tx2 <- function(info) {
  eff <- sub(".*EFF=([^;]*).*", "\\1", info)                                   # toma el texto de EFF=
  entradas <- regmatches(eff, gregexpr("[A-Za-z0-9_]+\\([^)]*\\)", eff))[[1]]  # separa cada efecto
  e <- entradas[grepl("AT5G61420.2", entradas, fixed = TRUE)]                  # se queda con el del .2
  if (length(e) == 0) return(c("intergenica", ""))                             # si no hay, es intergénica
  tipo <- sub("\\(.*", "", e[1])                                               # nombre del efecto
  aa <- regmatches(e[1], regexpr("p\\.[A-Za-z]{3}[0-9]+[A-Za-z]{3}", e[1]))    # cambio de aminoácido
  c(tipo, ifelse(length(aa) == 0, "", aa))                                     # devuelve tipo y cambio
}

anot <- t(sapply(vcf$INFO, efecto_tx2, USE.NAMES = FALSE))   # aplica la función a cada SNV

# Traducción de los nombres de SnpEff al español
traduce <- c(upstream_gene_variant   = "aguas arriba (promotor)",
             missense_variant        = "cambio de sentido",
             synonymous_variant      = "sinónima",
             stop_gained             = "sin sentido",
             intron_variant          = "intrónica",
             `3_prime_UTR_variant`   = "3'UTR",
             `5_prime_UTR_variant`   = "5'UTR",
             downstream_gene_variant = "aguas abajo",
             intergenica             = "intergénica")

variantes <- data.frame(
  pos     = vcf$POS,                            # posición en el cromosoma 5
  ref     = vcf$REF,                            # base de referencia
  alt     = vcf$ALT,                            # base alternativa
  tipo    = unname(traduce[anot[, 1]]),         # tipo de variante en español
  cambio  = anot[, 2],                          # cambio de aminoácido (si aplica)
  llamada = rowMeans(!is.na(dosis)),            # proporción de accesiones con genotipo
  q       = rowMeans(dosis, na.rm = TRUE) / 2   # frecuencia del alelo ALT
)
variantes$maf <- pmin(variantes$q, 1 - variantes$q)   # frecuencia del alelo menor (MAF)

polim <- variantes %>% filter(maf > 0)   # solo las que varían en las 729 accesiones
cat("SNV polimórficas:", nrow(polim), "(esperado: 166)\n")
print(table(polim$tipo))                 # cuántas hay de cada tipo

comunes <- polim %>% filter(llamada >= 0.95, maf >= 0.05)   # filtro: llamada >= 95 % y MAF >= 5 %
cat("SNV comunes para asociación:", nrow(comunes), "(esperado: 10)\n")

dos_c <- dosis[as.character(comunes$pos), , drop = FALSE]   # matriz de dosis solo de las comunes
cat("✔ SECCIÓN 5 OK\n\n")


# =============================================================================
# SECCIÓN 6. P1 · ASOCIACIÓN DE CADA SNV CON GSL ALIFÁTICOS TOTALES
# =============================================================================

muestras$y <- log10(muestras$gsl_alif)   # fenotipo en log10 para reducir asimetría

asoc <- do.call(rbind, lapply(seq_len(nrow(comunes)), function(i) {   # repite para cada SNV común
  d <- muestras                                  # copia de la tabla
  d$x <- dos_c[i, ] / 2                          # genotipo: 0 = REF/REF, 1 = ALT/ALT
  m0 <- lm(y ~ x, data = d)                      # modelo SIN covariables
  m1 <- lm(y ~ x + Admixture_Group + MAM + AOP, data = d)   # modelo AJUSTADO por estructura y loci
  c0 <- summary(m0)$coefficients["x", ]          # resultados del genotipo, modelo sin covariables
  c1 <- summary(m1)$coefficients["x", ]          # resultados del genotipo, modelo ajustado
  data.frame(pos = comunes$pos[i], tipo = comunes$tipo[i], maf = round(comunes$maf[i], 3),
             p_sin_cov  = c0[4],                 # p sin covariables
             beta       = c1[1],                 # efecto del homocigoto ALT (log10)
             ee         = c1[2],                 # error estándar
             p          = c1[4],                 # p ajustada
             cambio_pct = (10^c1[1] - 1) * 100)  # efecto expresado como % de cambio
}))

# Número efectivo de pruebas (las SNV están ligadas, no son independientes)
G <- apply(t(dos_c), 2, function(g) { g[is.na(g)] <- mean(g, na.rm = TRUE); g })  # rellena NA con la media
ev <- eigen(cor(G), only.values = TRUE)$values   # valores propios de la correlación entre SNV
m_ef <- sum(ev)^2 / sum(ev^2)                    # número efectivo de pruebas
alfa <- 0.05 / m_ef                              # alfa corregido (Bonferroni efectivo)
cat("Pruebas efectivas:", round(m_ef, 1), "- alfa corregido:", round(alfa, 4), "(esperado: 3.8 y 0.013)\n")

asoc$FDR <- p.adjust(asoc$p, method = "BH")   # corrección de Benjamini-Hochberg
asoc <- asoc %>% arrange(p)                   # ordena de menor a mayor p
print(asoc, digits = 3)                       # Tabla 1 del reporte

m_cov <- lm(y ~ Admixture_Group + MAM + AOP, data = muestras)   # modelo solo con covariables
cat("R2 de las covariables:", round(summary(m_cov)$r.squared, 3), "(esperado: 0.073)\n")
print(round(tapply(muestras$gsl_alif, muestras$grupo, mean), 2))  # GSL medio por región
print(round(quantile(muestras$gsl_alif, c(0.25, 0.5, 0.75)), 1))  # mediana y cuartiles
cat("✔ SECCIÓN 6 OK\n\n")


# =============================================================================
# SECCIÓN 7. P1 · VARIANTES RARAS (PRUEBA DE CARGA) Y PODER ESTADÍSTICO
# =============================================================================

pos_mis <- polim$pos[polim$tipo == "cambio de sentido"]        # posiciones de cambio de sentido
d_mis <- dosis[as.character(pos_mis), , drop = FALSE]          # sus genotipos
d_mis[is.na(d_mis)] <- 0                                       # sin llamada = se cuenta como REF
muestras$portadora <- as.numeric(colSums(d_mis >= 1) > 0)      # 1 si tiene al menos una variante
cat("Accesiones portadoras:", sum(muestras$portadora), "(esperado: 87)\n")

carga <- lm(y ~ portadora + Admixture_Group + MAM + AOP, data = muestras)  # modelo de carga
print(summary(carga)$coefficients["portadora", ])                          # efecto y p (esperado p = 0.59)
print(round(tapply(muestras$gsl_alif, muestras$portadora, mean), 2))       # medias: no portadora / portadora

sigma <- summary(m_cov)$sigma   # desviación estándar residual
n <- nrow(muestras)             # tamaño de muestra
for (f in c(0.05, 0.10)) {      # para MAF de 5 % y 10 %
  ee  <- sigma / sqrt(n * f * (1 - f))                  # error estándar esperado del efecto
  mde <- (qnorm(1 - alfa / 2) + qnorm(0.80)) * ee        # efecto mínimo detectable (80 % de poder)
  cat("MAF", f, "-> efecto mínimo detectable:", round((10^mde - 1) * 100, 1), "%\n")  # 12.5 % y 8.9 %
}
cat("✔ SECCIÓN 7 OK\n\n")


# =============================================================================
# SECCIÓN 8. P1 · ANÁLISIS SECUNDARIO: 4MSOB Y FRACCIÓN DE GSL C4
# =============================================================================

secundario <- function(variable) {   # repite el modelo ajustado con otra variable respuesta
  out <- do.call(rbind, lapply(seq_len(nrow(comunes)), function(i) {
    d <- muestras; d$x <- dos_c[i, ] / 2; d$v <- variable                     # genotipo y respuesta
    co <- summary(lm(v ~ x + Admixture_Group + MAM + AOP, data = d))$coefficients["x", ]
    data.frame(pos = comunes$pos[i], tipo = comunes$tipo[i], beta = co[1], p = co[4])
  }))
  out$FDR <- p.adjust(out$p, method = "BH")   # corrección FDR
  out[order(out$p), ]                          # ordena por p
}

sec_4mso <- secundario(log10(muestras$gsl_4mso + 0.01))   # glucorafanina (+0.01 por los ceros)
sec_c4   <- secundario(muestras$frac_c4)                  # fracción de GSL de 4 carbonos
print(sec_4mso, digits = 3)   # esperado: ninguna asociación
print(sec_c4,   digits = 3)   # esperado: 24692910 con beta = +0.031 y p = 1.9e-4

muestras$g_24692910 <- dosis["24692910", ]                # genotipo de la SNV del bloque centroeuropeo
print(muestras %>%
        filter(!is.na(g_24692910)) %>%                    # quita las accesiones sin llamada
        group_by(MAM, g_24692910) %>%                     # agrupa por MAM y genotipo
        summarise(frac_c4 = round(mean(frac_c4), 3),      # fracción C4 promedio
                  n = n(), .groups = "drop"))             # Tabla 2 del reporte
cat("✔ SECCIÓN 8 OK\n\n")


# =============================================================================
# SECCIÓN 9. P3 · FRECUENCIAS, HARDY-WEINBERG, F Y FST
# =============================================================================

hardy_weinberg <- function(g) {      # recibe genotipos 0/1/2 de un grupo
  g   <- g[!is.na(g)]                # quita los faltantes
  nAA <- sum(g == 0)                 # homocigotos REF
  nAa <- sum(g == 1)                 # heterocigotos
  naa <- sum(g == 2)                 # homocigotos ALT
  N   <- nAA + nAa + naa             # total de individuos
  q   <- (nAa + 2 * naa) / (2 * N)   # frecuencia del alelo ALT
  p   <- 1 - q                       # frecuencia del alelo REF
  Ho  <- nAa / N                     # heterocigosis observada
  He  <- 2 * p * q                   # heterocigosis esperada bajo HW
  obs <- c(nAA, nAa, naa)            # genotipos observados
  esp <- c(p^2, 2 * p * q, q^2) * N  # genotipos esperados bajo HW
  chi <- sum(((obs - esp)^2 / esp)[esp > 0])   # chi-cuadrada (gl = 1)
  data.frame(AA = nAA, Aa = nAa, aa = naa, N = N, p = p, q = q, Ho = Ho, He = He,
             F = ifelse(He > 0, 1 - Ho / He, NA),          # coeficiente de endogamia
             E_AA = esp[1], E_Aa = esp[2], E_aa = esp[3],  # esperados
             chi2 = chi, p_valor = pchisq(chi, 1, lower.tail = FALSE))
}

grupos <- c("Norte", "Centro", "Sur")

tabla_p3 <- do.call(rbind, lapply(seq_len(nrow(comunes)), function(i) {   # para cada SNV común
  g  <- dos_c[i, ]                                                         # sus genotipos
  pg <- lapply(grupos, function(k) hardy_weinberg(g[muestras$grupo == k])) # HW en cada grupo
  qs <- sapply(pg, `[[`, "q")                                              # frecuencia ALT por grupo
  ns <- sapply(pg, `[[`, "N")                                              # tamaño por grupo
  qb <- sum(qs * ns) / sum(ns)                                             # frecuencia promedio ponderada
  HT <- 2 * qb * (1 - qb)                                                  # heterocigosis total esperada
  HS <- sum(2 * qs * (1 - qs) * ns) / sum(ns)                              # heterocigosis dentro de grupos
  alt <- sapply(pg, function(x) x$Aa + 2 * x$aa)                           # alelos ALT por grupo
  ct  <- suppressWarnings(chisq.test(rbind(alt, 2 * ns - alt)))            # chi2 entre grupos (gl = 2)
  tot <- hardy_weinberg(g)                                                 # HW con las 729 juntas
  data.frame(pos = comunes$pos[i], tipo = comunes$tipo[i],
             q_Norte = qs[1], q_Centro = qs[2], q_Sur = qs[3],
             F = tot$F, chi2_HW = tot$chi2,
             FST = (HT - HS) / HT,                                         # diferenciación entre grupos
             chi2_grupos = unname(ct$statistic), p_grupos = ct$p.value)
}))
tabla_p3 <- tabla_p3 %>% arrange(desc(FST))   # ordena de mayor a menor FST
print(tabla_p3, digits = 3)                   # Tabla 3 del reporte

cat("Heterocigotos:", sum(dosis == 1, na.rm = TRUE), "de", sum(!is.na(dosis)), "genotipos\n")  # esperado 0

q119 <- dosis["24690219", ]                           # SNV Q119Q para el cálculo completo
print(hardy_weinberg(q119), digits = 4)               # total: 631 / 0 / 65, chi2 = 696
for (k in grupos) {                                   # por grupo
  cat("--", k, "--\n")
  print(hardy_weinberg(q119[muestras$grupo == k]), digits = 4)
}
cat("✔ SECCIÓN 9 OK\n\n")


# =============================================================================
# SECCIÓN 10. P3 · DESEQUILIBRIO DE LIGAMIENTO ENTRE SNV COMUNES
# =============================================================================

ld <- cor(t(dos_c), use = "pairwise.complete.obs")^2   # r2 entre cada par de SNV
print(round(ld, 2))                                    # matriz de r2
# Bloque norteño: 24690219, 24690646, 24692594 (r2 = 0.83–0.98)
# Bloque centroeuropeo: 24690018, 24692597, 24692910, 24693146 (r2 = 1.00)
cat("✔ SECCIÓN 10 OK\n\n")


# =============================================================================
# SECCIÓN 11. P2 · CRUCE MENDELIANO (SIMULADO: no son datos reales)
# =============================================================================

f2 <- c(REF_REF = 54, REF_ALT = 97, ALT_ALT = 49)          # genotipos observados en 200 plantas F2
prueba <- chisq.test(f2, p = c(0.25, 0.50, 0.25))          # chi2 contra la proporción 1:2:1
print(prueba$expected)                                     # esperados: 50, 100, 50
print(prueba)                                              # chi2 = 0.43, gl = 2, p = 0.81

media_REF <- 0.230                                         # fracción C4 del parental REF/REF
media_ALT <- 0.347                                         # fracción C4 del parental ALT/ALT
media_F1  <- 0.303                                         # fracción C4 de la F1
m <- (media_REF + media_ALT) / 2                           # punto medio parental
a <- (media_ALT - media_REF) / 2                           # efecto aditivo
d <- media_F1 - m                                          # desviación por dominancia
t_d <- d / 0.015                                           # t de d = 0 (EE de la simulación)
cat("m =", round(m, 3), "| a =", round(a, 3), "| d =", round(d, 3),
    "| d/a =", round(d / a, 2), "| t =", round(t_d, 2),
    "| p =", round(2 * pt(-abs(t_d), df = 57), 2), "\n")   # d/a = 0.25, p = 0.35 -> aditividad
cat("✔ SECCIÓN 11 OK\n\n")


# =============================================================================
# SECCIÓN 12. GUARDAR TABLAS Y FIGURAS EN "resultados/"
# =============================================================================

write.csv(polim,    "resultados/variantes_MYB28.csv",   row.names = FALSE)   # catálogo de variantes
write.csv(asoc,     "resultados/asociacion_P1.csv",     row.names = FALSE)   # Tabla 1
write.csv(sec_c4,   "resultados/secundario_fracC4.csv", row.names = FALSE)   # análisis secundario
write.csv(tabla_p3, "resultados/frecuencias_P3.csv",    row.names = FALSE)   # Tabla 3

png("resultados/fig_frecuencias_P3.png", width = 1600, height = 900, res = 150)  # abre la figura
mat <- t(as.matrix(tabla_p3[, c("q_Norte", "q_Centro", "q_Sur")]))              # frecuencias por grupo
colnames(mat) <- tabla_p3$pos                                                   # nombre = posición
barplot(mat, beside = TRUE, las = 2, ylab = "Frecuencia del alelo ALT (q)",
        legend.text = grupos, args.legend = list(x = "topright", bty = "n"),
        main = "Frecuencias de SNV de MYB28 por región")                        # gráfica de barras
dev.off()                                                                       # cierra y guarda

png("resultados/fig_gsl_por_grupo.png", width = 1200, height = 900, res = 150)
boxplot(gsl_alif ~ grupo, data = muestras, xlab = "",
        ylab = "GSL alifáticos (nmol/mg semilla)",
        main = "GSL alifáticos totales por región")                             # diagrama de caja
dev.off()

print(list.files("resultados"))            # lista lo que se guardó
browseURL(normalizePath("resultados"))     # abre la carpeta en el explorador de Windows
cat("✔ SECCIÓN 12 OK: todo listo\n")