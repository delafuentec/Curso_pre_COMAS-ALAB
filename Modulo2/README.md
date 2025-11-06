# Datos de genotipos: Control de Calidad y Parentesco

## Objetivos
En esta sesión, vamos a:
- repasar la organización de datos de genotipos.
- realizar unos pasos de control de calidad.
- estimar el grado de parentesco entre individuos.
- generar datos listos para llevar a cabo análisis de genética de poblaciones.

## Intrudcción a los datos
Vamos a analizar los datos de genotipos de individuos antiguos y modernos del Cono Sur de Ámerica, y particularmente del Sur de Patagonia.
Existen muchos formatos para este tipo de dato. Vamos a enfocarnos en dos de los más comunes: [</em>plink binary</em>](https://www.cog-genomics.org/plink2/input) y [<em>eigenstrat</em>](https://reich.hms.harvard.edu/software/InputFileFormats).<br>
Ambos formatos, se constituyen de 3 ficheros: 
- uno con la matriz de `genotipos` que puede ser binario, es decir no lisible como texto (`.bed` para <em>plink</em> o `.geno` para <em>packed eigenstrat</em>). 
- uno con la información de las `variantes` analizadas (`.bim` para <em>plink</em> y `.snp` para <em>eigenstrat</em> ).
- uno con la información de los `individuos` analizados (`.fam` para <em>plink</em> y `.ind` para <em>eigenstrat</em>).

Para <em>eigenstrat</em>, se suele añadir el sufijo `.txt` a los 3 ficheros asociados al formato no binario (es decir `<preffile>.geno.txt`, `<preffile>.snp.txt`, `<preffile>.ind.txt`).

### Fichero de anotación de las variantes
Los ficheros `.bim` y `.snp` contienen una linea por variante y 6 columnas con la misma información. La única diferencia es el orden de las 2 primeras columnas:
- ID de la variante: columna 1 para .snp y 2 para .bim
- cromosoma: columna 2 para .snp y 1 para.bim
- posición genética: columna 3. Se trata de la posición de la variante en función de los eventos de recombinación pasados. Para este curso, no es relevante.
- posición física: columna 4. Se trata del par de base en cual se encuentra la variante en la secuencia del cromosoma.
- Alelo 1: columna 5. Suele ser el alelo de menor frecuencia en la muestra analizada
- Alelo 2: columna 6. Ser ser el alelo de mayor frecuencia en la muestra analizada.

Noten entonces que ambos formatos son compatibles solamente para variantes bialélicas.

### Fichero de anotación de las individuos
Los ficheros `.ind` (<em>eigenstrat</em>) contienen 3 columnas, en este orden:
- ID del Individuo
- Sexo (M, F o U)
- Grupo poblacional.

Los ficheros <em>plink</em> son un poco más completos porque están pensados más por estudios de genética epidemiologíca. Contiene 6 columnas:
- ID de la familia del individuo (se suele poner el grupo poblacional en los estudios de antropología molecular )
- ID del individuo
- ID del padre biológico ('0' si no presente en el conjunto de datos)
- ID de la madre biológica ('0' si no presente en el conjunto de datos)
- Sexo biológico  ('1' = hombre, '2' = mujer, '0' = desconocido)
- Valor fenotípico (se suele poner también en esta columna el grupo poblacional en los estudios de antropología molecular ). 

### Fichero de matriz de genotipos
El fichero de genotipos de <em>eigenstrat</em>, el fichero de genotipos se constituye de 1 línea por variantes, con 1 columna por individuo (columna de 1 caracater <strong>no</strong> separadas).<br>
El genotipo de un individuo a una variante se codifica por el número de copias de este individuo del Alelo1 de esta variante como aparece en el archivo `.snp.txt`. 9 significa valor faltante <br>

### Ejemplo

Tenemos 3 individuos (IndA, IndB y IndC) genotipifacos en 4 variantes (snp1, snp2, snp3, snp4).
Archivo `.snp.txt`:
snp1  1 0.1001  100000  A C<br>
snp2  1 0.6162  600000  T G<br>
snp3  2 0.5125  513341  C G<br>
snp4  4 0.1512  251334  G A<br>


Archivo `.ind.txt`:
IndA  M Pop1<br>
IndB  M Pop1<br>
IndC  F Pop2<br>
IndD  M Pop3<br>

Archivo `.geno.txt`:
0192
1122
2222
0211

Permite saber que el individuo IndA de la pop1 tiene el genotipo C/C para snp1, A/C para snp2, Faltante para snp3 y A/A para snp4





Genotypes coded as 0, 1, 2 (number of derived alleles), or 9 (missing)
### Leer los datos
Existen varios paquetes de R que permiten leer y procesar estos formatos, pero vamos a usar <em>dartR.base</em> que permite leer el formato plink.

```
require(dartRverse)
require(dartR.base)
require(ade4)

### you maybe need to change the directory below
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB")
datosGeno<-gl.read.PLINK("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.plink",ind.metafile="PatagoniaDataSetWithOutgroups_ALAB2025/finalSet.metadataPerind.txt")
```

Al leer los datos, en <em>verbose</em>, vemos algunos mensajes (si hay monomorfismos, si hay loci sin datos, etc.).
En algunas versiones de dartR.base, puede aparecer el mensaje:<br>
<em>Warning: Locus metafile not provided, locus metrics will be
        calculated where this is possible</em>.<br>
En este caso, se puede solucionar corriendo lo siguiente.

```### fix allele names
bim <- read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.plink.bim", header = FALSE)
# Columns: CHR SNP CM POS A1 A2
datosGeno@loc.all <- paste(bim$V5, bim$V6, sep="/")```








