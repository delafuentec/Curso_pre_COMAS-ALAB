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
Archivo `.snp.txt`:<br>
snp1  1 0.1001  100000  A C<br>
snp2  1 0.6162  600000  T G<br>
snp3  2 0.5125  513341  C G<br>
snp4  4 0.1512  251334  G A<br>


Archivo `.ind.txt`:<br>
IndA  M Pop1<br>
IndB  M Pop1<br>
IndC  F Pop2<br>
IndD  M Pop3<br>

Archivo `.geno.txt`:<br>
0192<br>
1122<br>
2222<br>
0211<br>

Permite saber que el individuo IndA de la pop1 tiene el genotipo C/C para snp1, A/C para snp2, Faltante para snp3 y A/A para snp4

## Leer los datos
Existen varios paquetes de R que permiten leer y procesar estos formatos, pero vamos a usar <em>dartR.base</em> que permite leer el formato plink con la función `gl.read.PLINK`, a la cual le damos también como input un fichero de metadato para individuo llamado `finalSet.metadataPerind.txt`

```
require(dartRverse)
require(dartR.base)
require(ade4)

### you maybe need to change the directory below
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB")
datosGeno<-gl.read.PLINK("PatagoniaDataSetWithOutgroups_ALAB2025/test.plink")
```

Al leer los datos, en <em>verbose</em>, vemos algunos mensajes (si hay monomorfismos, si hay loci sin datos, etc.).
En algunas versiones de dartR.base, puede aparecer el mensaje:<br>
<em>The slot loc.all, which stores allele name for each locus, is empty. Creating a dummy variable (A/C) to insert in this slot.</em>.<br>
En este caso, se puede solucionar corriendo lo siguiente.

```
### fix allele names
bim <- read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.plink.bim", header = FALSE)
# Columns: CHR SNP CM POS A1 A2
datosGeno@loc.all <- paste(bim$V5, bim$V6, sep="/")
```

Vamos a asignar los nombres de las poblaciones que por ahora están guardadas como "Family" en `ind.metrics`
``` 
datosGeno$pop<-as.factor(datosGeno$other$ind.metrics$Family)
```


## Explorar los datos
Podemos ver la tasa de individuos con datos por locus, y la frecuencia alélica haciendo lo siguiente:
 ```
hist(datosGeno@other$loc.metrics$CallRate,main="Call Rate per Locus")
hist(datosGeno@other$loc.metrics$maf,main="Minor Allele Frequency per Locus",n=100)
```
Vemos que hay locis con un call rate bajo (<0.2) y que son monomórficos (MAF=0).

Podemos ver la tasa de loci con datos por individuo.
 ```
hist(datosGeno@other$ind.metric$Call.rate,main="Call Rate per individual",n=10)
```
En cuanto a la tasa de Heterosigosidad, vemos que son todos los individuos homocigotos:
```
table(datosGeno@other$ind.metrics$Heterozygosity)
```

Esto se debe a que, en el caso del ADN antiguo, la baja cobertura suele obligar a generar datos seudo-haploides. En la práctica, esto significa que para cada posición genómica y cada individuo se selecciona al azar una de las lecturas disponibles, y el nucleótido observado en ella se asigna como un genotipo homocigota para ese individuo.

## Filtrar
Algunos análisis pueden requerir filtrar variantes e individuos por valores faltantes. 
De momento vamos solamente guardar en dos vectores las posiciones monomórficas y las que tienen datos para solo 1 individuo (Call Rate inferior a 2/64 = 0.03125), y luego los individuos con datos para más de 5000 variantes restantes.<br>
Noten que los siguientes comandos, por razón de memoria, no guardan los datos filtrados, y solo generar las listas de SNPs e Individuos que tendríamos que sacar (volveremos alas mismas luego).
```
###filtro de SNPs
SNPs_to_exclude=datosGeno$other$loc.metrics$AlleleID[ datosGeno$other$loc.metrics$CallRate<0.03125 |
                                                      datosGeno$other$loc.metrics$maf==0 ]  
### get number of SNPs still included 
numSNPs=datosGeno$n.loc-length(SNPs_to_exclude)
###get inds metrics 
metricsInds_forKeptSNPs <- gl.drop.loc(datosGeno, loc.list=SNPs_to_exclude)$other$ind.metrics
###add a column of the number of called snps (Call.rate*numSNPs)
metricsInds_forKeptSNPs$CalledSNPs=metricsInds_forKeptSNPs$Call.rate*numSNPs
###save in a vector individuals with less than 50000 SNPs
INDs_to_exclude=metricsInds_forKeptSNPs$id[metricsInds_forKeptSNPs$CalledSNPs< 50000]
```


## Parentesco

Existen varios métodos para detectar individuos aparentados. Vamos a usar [<em>BREADR</em>](https://joss.theoj.org/papers/10.21105/joss.07916). 
Este método es una extensión al gold-standard en el campo ([<em>READv2</em>](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03350-3)), ya que se base en el escore PMR.
Este escore es la proporción de sitios superpuestos con genotipos no coincidentes, es decir la fracción de posiciones del genoma en las que dos individuos tienen datos genotípicos válidos (no faltantes) y sus llamadas de genotipo difieren. En otras palabras, mide qué porcentaje de los sitios comparables entre ambos individuos presentan genotipos distintos, reflejando el grado de discrepancia genética entre ellos.
Las diferencias más importantes con READv2 son:
1. el uso de sitios independientes (separados por recombinación durante la meiosis). Por defecto, el programa usa un locus cada 10,000 par de bases (paramétro ` filter_length`).
2. su implementación en R.

Noten que estos métodos implementan aproximaciones que tomen en cuenta las limitaciones de los datos en ADN antiguo (pseudo-haploides y con baja cobertura). Para analizar datos modernos que no padecen de estas limitaciones, se recomienda usar otros métodos, por ejemplo [<em>KING</em>](https://www.kingrelatedness.com/).
<em>BREADR</em>  usa el formato <em>eigenstrat</em>. Vamos a analizar solamente los individuos antiguos.
Primero se calculan las métricas PMR para cada par de individuo, es la parte más demandante computacionalmente, por lo cual, vamos a leer directamente la tabla de resultado.

```
require(BREADR)
require(tibble)

prefEigenstrat="PatagoniaDataSetWithOutgroups_ALAB2025/Ancient"

## we can read directly the table that would be generated by the "preprocess" which ios 
if(! file.exists("Modulo2/PMR_allAncientTogether.txt")){
  countsPMR <- processEigenstrat(
    indfile = paste(prefEigenstrat,".ind.txt",sep=""),
    snpfile = paste(prefEigenstrat,".snp.txt",sep=""),
    genofile = paste(prefEigenstrat,".geno.txt",sep=""),
    outfile = "Modulo2/PMR_allAncientTogether.txt")
  
}else{
  countsPMR<-as_tibble(read.table("Modulo2/PMR_allAncientTogether.txt",stringsAsFactors = F,header=T,sep="\t"))
}

print(head(countsPMR))

```
Se genera entonces una tabla, con 4 columnas: par analizado | número de variantes | numero de variantes disconcordantes | PMR
Ahora vamos a ver cuales son los individuos supuestamente aparentados, es decir los que tienen un PMR bajo.
Para definir 
```
##estimate relatedness
relatedness_allAncientTogether <- callRelatedness(countsPMR)

### add population label
##define some functions for that first
#This function will be used to split a string containing pairs of individuals
change<-function(string,split,pos){
  return(strsplit(string,split=split)[[1]][pos])
}
#This function will be used to add 2 columns to the relatedness results table: "Pops" (containing pop of each individual) and "Regions" (conaining region of each individual)
addPops<-function(tableRel,metaData){
  tableRel$Ind1=sapply(tableRel$pair,change,split=" - ",pos=1)
  tableRel$Pop1 <- metaData$pop[match(tableRel$Ind1, metaData$id)]
  tableRel$Region1 <- metaData$Region[match(tableRel$Ind1, metaData$id)]
  tableRel$Ind2=sapply(tableRel$pair,change,split=" - ",pos=2)
  tableRel$Pop2 <- metaData$pop[match(tableRel$Ind2, metaData$id)]
  tableRel$Region2 <- metaData$Region[match(tableRel$Ind2, metaData$id)]
  
  return(tableRel)
}
# We read the metadata
meta<-read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.metadataPerind.txt",sep="\t",header=T)
#and we use the above functions to add Pops and Regions
relatedness_allAncientTogether<-addPops(tableRel=relatedness_allAncientTogether,
                                        metaData=meta
                                        )

```
Vamos ahora a ver que resultados obtuvimos. Cuantos pares de individuos aparentados, en qué región y población.

```
###Number of related pairs indivuduals:
related_allAncientTogether<-relatedness_allAncientTogether[ relatedness_allAncientTogether$relationship!="Unrelated",]
nrow(related_allAncientTogether)
table(related_allAncientTogether$relationship)

### we see a lot of related pairs. Let's see in which region
table(related_allAncientTogether$Region1,related_allAncientTogether$Region2)

###Let's see in which populations we have pairs of 1st-degree or Twins/Duplicated
first_allAncientTogether<-related_allAncientTogether[ related_allAncientTogether$relationship %in% c("Same_Twins","First_Degree"),]
table(first_allAncientTogether$Pop1,first_allAncientTogether$Pop2)

### We can already see some temporal/spatial inconsistencies 
```


Vemos entonces inconsistencias. El problema con el análisis de parentesco es que definir si el escore usado para medir la afinidad genética entre individuos es  (PMR en el caso de <em>BREADR</em>) significa parentesco depende de la diversidad genética esperada en la población.
Lo que hace <em>BREADR</em> por defecto es usar la mediana de los PMR observados para todos los pares de individuos analizados.
Sin embargo, en poblaciones con tamaño poblacional reducido, como es el caso de las patagónicas, se espera un escore PMR entre individuos mucho menor que en poblaciones de tamaño más grande, como las de Centro Andes.
Entonces, si comparamos todas las poblaciones a la vez, es probable que se detecte muchos pares de individuos aparentados en poblaciones de Patagonia.
Para evitar estos falsos positivos podemos:
1. realizar las inferencias de parentesco por región pasando al método lo que llaman escore "Promedio de parentesco", pero que en realidad es la mediana observada en todos los pares de individuos posibles en la muestra analizada.
2. pasar al método solo la tabla de escore PMR para el subconjunto de individuos de intéres.
Esta ultima opción es la que vamos a hacer ahora: vamos a ir, región por región, generar la sub-tabla, verificar que hay suficientes pares (digamos por lo menos 5) y en este caso proceder.
Se va a generar una lista, con un elemento por región siendo la tabla de inferencia de parentesco.

```
print(unique(meta$Region))
### regions were actually divided temporaly, we will not do so... 
listRelatednessPerRegion<-list()
for(region in c("SouthPatagonia","CentralChile","Beringia","SouthernNorthAmerica","Brazil","CentralAndes","Pampa")){
  print(paste("analyzing",region))
  tmpRegion<-relatedness_allAncientTogether[ grepl(region,relatedness_allAncientTogether$Region1) & grepl(region,relatedness_allAncientTogether$Region2),c("pair","nsnps","mismatch","pmr")]
  npairs<-nrow(tmpRegion)
  print(paste(npairs,"individual pairs to be analyzed"))
  if(npairs <6){
    print("not enough... skip")
    next
  }
  listRelatednessPerRegion[[region]]<-callRelatedness(tmpRegion)
  print(paste(sum(listRelatednessPerRegion[[region]]$relationship!="Unrelated"),"pairs of related individuals found"))
  
}
```
Vemos que reducimos drastícamente el número de pares de aparentados encontrados. 
vamos a generar una tabla que resume las relaciones al primer grado (en aDNA antiguo se suelen guardar los individuos aparentados al segundo grado). 

```
tableSUM<-c()
for(region in names(listRelatednessPerRegion)){
  ###for each region, keep only lines for 1st degree relationship 
  tmp1st<-listRelatednessPerRegion[[region]]
}



