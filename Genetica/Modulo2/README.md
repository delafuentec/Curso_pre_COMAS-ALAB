# Datos de genotipos: Control de Calidad y Parentesco

## Objetivos
En esta sesión vamos a:  
- repasar la organización de los datos de genotipos,  
- realizar algunos pasos de control de calidad,  
- estimar el grado de parentesco entre individuos, y  
- generar datos listos para realizar análisis de genética de poblaciones.  

---

## Introducción a los datos
Analizaremos datos de genotipos de individuos antiguos y modernos del Cono Sur de América, en particular del sur de la Patagonia.  
Existen muchos formatos para este tipo de datos, pero nos enfocaremos en dos de los más comunes: [**PLINK binary**](https://www.cog-genomics.org/plink2/input) y [**EIGENSTRAT**](https://reich.hms.harvard.edu/software/InputFileFormats).  

Ambos formatos se componen de tres ficheros principales:
- un archivo con la matriz de **genotipos**, que puede ser binario (no legible como texto): `.bed` para *PLINK* o `.geno` para *packed EIGENSTRAT*,  
- un archivo con la información de las **variantes** analizadas (`.bim` para *PLINK* y `.snp` para *EIGENSTRAT*),  
- y un archivo con la información de los **individuos** analizados (`.fam` para *PLINK* y `.ind` para *EIGENSTRAT*).  

En el caso de *EIGENSTRAT*, es común añadir el sufijo `.txt` a los tres ficheros cuando se utilizan en formato no binario (por ejemplo: `<prefijo>.geno.txt`, `<prefijo>.snp.txt`, `<prefijo>.ind.txt`).

---

### Fichero de anotación de las variantes
Los archivos `.bim` y `.snp` contienen una línea por variante y seis columnas con información equivalente.  
La única diferencia entre ellos es el orden de las dos primeras columnas:

| Información         | `.snp`     | `.bim`     |
|---------------------|------------|-------------|
| ID de la variante   | Columna 1  | Columna 2   |
| Cromosoma           | Columna 2  | Columna 1   |
| Posición genética   | Columna 3  | Columna 3   |
| Posición física     | Columna 4  | Columna 4   |
| Alelo 1             | Columna 5  | Columna 5   |
| Alelo 2             | Columna 6  | Columna 6   |

- **Posición genética:** indica la ubicación de la variante según los eventos históricos de recombinación (no será relevante para este curso).  
- **Posición física:** corresponde al par de bases donde se encuentra la variante en la secuencia del cromosoma.  
- **Alelo 1:** suele ser el alelo de menor frecuencia en la muestra.  
- **Alelo 2:** suele ser el alelo de mayor frecuencia.  

> *Atención: Ambos formatos solo son compatibles con variantes **bialélicas**.*

### Fichero de anotación de los individuos

Los ficheros `.ind` (*EIGENSTRAT*) contienen tres columnas, en el siguiente orden:  
1. **ID del individuo**  
2. **Sexo** (`M`, `F` o `U`)  
3. **Grupo poblacional**

Los ficheros de *PLINK* son algo más completos, ya que están diseñados originalmente para estudios de genética epidemiológica. Contienen seis columnas:  
1. **ID de la familia** del individuo (en estudios de antropología molecular suele usarse para indicar el grupo poblacional),  
2. **ID del individuo**,  
3. **ID del padre biológico** (`0` si no está presente en el conjunto de datos),  
4. **ID de la madre biológica** (`0` si no está presente),  
5. **Sexo biológico** (`1` = hombre, `2` = mujer, `0` = desconocido),  
6. **Valor fenotípico** (en estudios de antropología molecular también suele utilizarse para indicar el grupo poblacional).  

>  *En estudios poblacionales o arqueogenéticos, las columnas 1 y 6 de los ficheros `.fam` suelen utilizarse para almacenar información de grupo o contexto arqueológico.*

---

### Fichero de matriz de genotipos

En *EIGENSTRAT*, el archivo de genotipos contiene **una línea por variante** y **una columna por individuo**, con caracteres **no separados por espacios**.  
El genotipo de un individuo para una variante se codifica como el **número de copias del alelo 1** (según el archivo `.snp.txt`).  
El valor `9` indica un **dato faltante**.  

> *Por ejemplo, si el alelo 1 es `A`, el genotipo `0` corresponde a “no tiene A”, `1` a “heterocigota”, `2` a “homocigota A/A” y `9` a dato faltante.*

---

### Ejemplo

Supongamos que tenemos tres individuos (`IndA`, `IndB` y `IndC`) genotipificados en cuatro variantes (`snp1`, `snp2`, `snp3`, `snp4`).

**Archivo `.snp.txt`:**<br>:

snp1 1 0.1001 100000 A C<br>
snp2 1 0.6162 600000 T G<br>
snp3 2 0.5125 513341 C G<br>
snp4 4 0.1512 251334 G A<br>


**Archivo `.ind.txt`:**<br>

IndA M Pop1<br>
IndB M Pop1<br>
IndC F Pop2<br>
IndD M Pop3<br>


**Archivo `.geno.txt`:**<br>

0192<br>
1122<br>
2222<br>
0211<br>


Podemos interpretar que el individuo **IndA**, perteneciente a la **población 1**, presenta los siguientes genotipos:
- `snp1`: C/C  
- `snp2`: A/C  
- `snp3`: faltante  
- `snp4`: A/A


### Posiciones de Transversiones
### ¿Por qué analizamos solo transversiones?

En los análisis de **ADN antiguo**, es común restringirse a las **transversiones**, que son un tipo de cambio entre bases que **no involucra nucleótidos del mismo tipo químico**.  

Existen dos tipos de mutaciones puntuales:  

- **Transiciones**, que ocurren entre bases del mismo tipo:  
  - purinas ↔ purinas (*A ↔ G*)  
  - pirimidinas ↔ pirimidinas (*C ↔ T*)  
- **Transversiones**, que ocurren entre bases de distinto tipo:  
  - purinas ↔ pirimidinas (*A o G ↔ C o T*)

Las **transiciones** son más frecuentes, pero en el ADN antiguo pueden estar **enriquecidas artificialmente** debido al daño post mortem característico de este tipo de material genético.  
En particular, la **desaminación de citosinas** genera un exceso de mutaciones aparentes **C→T y G→A**, lo que puede distorsionar las estimaciones genéticas si no se corrige.  

Por esta razón, al limitar los análisis a **transversiones** —que **no se ven afectadas por estos procesos de daño**— se obtiene una señal genética más confiable y libre de artefactos de degradación.

Es lo que haremos en los modulos 2 y 3.

## Lectura de los datos

Existen varios paquetes de **R** que permiten leer y procesar estos formatos de genotipos.  
En este curso utilizaremos el paquete *dartR.base*, que incluye la función `gl.read.PLINK()`, la cual permite importar directamente archivos en formato **PLINK** en un objeto [`genlight`](https://rdrr.io/cran/adegenet/man/genlight.html).  

Además, esta función admite incorporar un archivo de metadatos por individuo (por ejemplo, `finalSet.metadataPerind.txt`), que contiene información complementaria como la población, la región o el estudio de origen.

```r
require(dartRverse)
require(dartR.base)
require(ade4)

### you maybe need to change the directory below
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB/Genetica/")
pref="PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.TVs.plink"
datosGeno<-gl.read.PLINK(pref)
```

Al leer los datos con la opción `verbose = TRUE`, se muestran varios mensajes informativos (por ejemplo, sobre loci monomórficos o con datos faltantes).  

En algunas versiones de *dartR.base*, puede aparecer el siguiente mensaje:  
*“The slot loc.all, which stores allele name for each locus, is empty. Creating a dummy variable (A/C) to insert in this slot.”*  

Este mensaje indica que la información sobre los alelos de cada locus está vacía, y el programa crea una variable ficticia (“A/C”) para completarla.  
Para corregirlo de forma manual, puede ejecutarse el siguiente comando:

```r
### fix allele names
bim <- read.table(paste(pref,".bim",sep=""), header = FALSE)
# Columns: CHR SNP CM POS A1 A2
datosGeno@loc.all <- paste(bim$V5, bim$V6, sep="/")
```

A continuación, vamos a asignar los nombres de las poblaciones a los individuos.   Esta información se encuentra almacenada en la columna **"Family"** del elemento `ind.metrics` de `datosGeno`. Vamos a transferir esos nombres para que queden correctamente asociados a la variable `pop`.
```r 
datosGeno$pop<-as.factor(datosGeno$other$ind.metrics$Family)
```


## Exploración de los datos

Podemos examinar la calidad general del conjunto de datos observando, por ejemplo:  
- la **proporción de individuos con datos disponibles por locus**, y  
- la **frecuencia alélica** de cada variante.  

Para ello, podemos ejecutar las siguientes funciones:
```r
hist(datosGeno@other$loc.metrics$CallRate,main="Call Rate per Locus")
hist(datosGeno@other$loc.metrics$maf,main="Minor Allele Frequency per Locus",n=100)
```
Observamos que existen loci con una **tasa de llamado** (*call rate*) baja (< 0.2) y muchas posiciones **monomórficas** (MAF = 0). 
> Usamos el panel 1240K que se basa en posiciones encontradas por ser polimorficas en una muestra mundial. Sin embargo, por el hecho que hubo mucha deriva génica durante la historia evolutiva de las poblaciones de Sudamerica (p.ej. con varios efectos fundadores), la diversidad genética en el sub-continente es reducida, particularmente en el Sur de Patagonia. Lo cual llevó a muchas variantesobservadas en el mundo siendo fijadas en las poblaciones que nos interesan. [de la Fuente et al. (2024)](https://pivotscipub.com/hpgg/4/1/0003) lo describieron en detalle.

También podemos explorar la **proporción de loci con datos disponibles por individuo**, para identificar muestras con una alta proporción de datos faltantes.
 ```r
hist(datosGeno@other$ind.metric$Call.rate,main="Call Rate per individual")
```
Ahora podemos miramos la **tasa de heterocigosidad** por cada individuo.
```r
table(datosGeno@other$ind.metrics$Heterozygosity)
```
> Observamos que todos los individuos presentan genotipos **homocigotos**.
> Esto se debe a que, en el caso del ADN antiguo, la baja cobertura suele obligar a generar datos seudo-haploides. En la práctica, esto significa que para cada posición genómica y cada individuo se selecciona al azar una de las lecturas disponibles, y el nucleótido observado en ella se asigna como un genotipo homocigota para ese individuo.




## Parentesco

### Introducción al análisis de parentesco

Para el análisis de parentesco que realizaremos a continuación, **no es necesario filtrar los datos por valores faltantes ni por frecuencia del alelo menor**, por lo que trabajaremos directamente con el conjunto completo de variantes.

Existen diversos métodos para detectar **individuos emparentados**. En este caso, utilizaremos [*BREADR*](https://joss.theoj.org/papers/10.21105/joss.07916), un método desarrollado específicamente para el análisis de parentesco en contextos de ADN antiguo.

*BREADR* extiende el enfoque del método ampliamente utilizado [*READv2*](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03350-3), basado en el **escore PMR** (*Pairwise Mismatch Rate*).  
El PMR representa la **proporción de sitios compartidos entre dos individuos en los que los genotipos no coinciden**. Es decir, mide la fracción de posiciones del genoma en las que ambos individuos tienen datos válidos (no faltantes) pero presentan genotipos distintos. En otras palabras, cuantifica el **grado de discrepancia genética** entre dos individuos.

Las principales diferencias entre *BREADR* y *READv2* son:

1. El uso de **sitios independientes** (separados por recombinación durante la meiosis).  
   Por defecto, el programa utiliza un locus cada 10 000 pares de bases (parámetro `filter_length`).

2. Su **implementación en R**, lo que facilita la integración con otros flujos de análisis.

Estos métodos incorporan **aproximaciones diseñadas para datos de ADN antiguo**, que suelen ser **pseudo-haploides y de baja cobertura**.  
En cambio, para analizar **datos genómicos modernos**, donde estas limitaciones no aplican, se recomienda emplear métodos alternativos como [*KING*](https://www.kingrelatedness.com/).

Finalmente, *BREADR* utiliza el formato de datos **EIGENSTRAT**. En este análisis consideraremos únicamente los **individuos antiguos**.

### Inferencia de los grados de parentesco

#### Generación de la métrica necesaria (PMR)

El primer paso consiste en calcular las **métricas PMR** (*Pairwise Mismatch Rate*) para cada par de individuos.  
Este proceso es **computacionalmente intensivo**, ya que implica comparar los genotipos entre todos los pares posibles del conjunto de datos.  
Por esta razón, en este ejercicio **leeremos directamente la tabla de resultados previamente generada**, en lugar de recalcular las métricas desde cero.

```r
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
Se genera entonces una tabla con cuatro columnas: el par de individuos analizado, el número total de variantes comparadas, el número de variantes discordantes y el valor de PMR.

#### Definir el grado de parentesco

Ahora vamos a identificar qué pares de individuos presentan parentesco y cuál es su grado. Para ello, se utilizan los valores de PMR esperados en la población, estimados a partir de la distribución de PMR de todos los pares analizados.

Se asume que la mayoría de los pares en el conjunto de datos no están emparentados (hasta segundo grado). Bajo esta suposición, la **mediana del PMR**, denotada como $\bar{p}$, constituye una estimación robusta del nivel de discrepancia genética esperado entre individuos no relacionados.

Siguiendo los planteamientos de *READv2*, se define el **valor medio esperado del PMR** para una relación de grado *k* = 0, 1 o 2 como:

\[
p_k = \bar{p} \left(1 - \frac{1}{k + 1}\right).
\]

Una vez definidos estos valores esperados, evaluamos si el PMR observado para un par de individuos es **significativamente menor** que cada uno de los valores esperados para los distintos grados de relación.  
El parentesco asignado será entonces el **menor valor de *k*** para el cual la diferencia **no resulta significativa**, tal como se ilustra en la figura siguiente.

<img width="729" height="344" alt="image" src="https://github.com/user-attachments/assets/d11688e8-3d94-48aa-bb11-d3e9cbef8cb9" />

A continuación, calcularemos el grado de parentesco para todos los pares de individuos incluidos en la muestra.

```r
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
Veamos cuántos pares de individuos aparecen como aparentados y a qué región y población pertenecen.

```r
###Number of related pairs indivuduals:
related_allAncientTogether<-relatedness_allAncientTogether[ relatedness_allAncientTogether$relationship!="Unrelated",]
nrow(related_allAncientTogether)
table(related_allAncientTogether$relationship)

### we see a lot of related pairs. Let's see in which region
table(related_allAncientTogether$Region1,related_allAncientTogether$Region2)

###Let's see in which groups we have pairs of up to 2-degree
second_allAncientTogether<-related_allAncientTogether[ related_allAncientTogether$relationship %in% c("Same_Twins","First_Degree","Second_Degree"),]
table(second_allAncientTogether$Pop1,second_allAncientTogether$Pop2)

### We can already see some temporal/spatial inconsistencies 
```

Vemos entonces que se detectan muchos pares de individuos aparentados en el sur de Patagonia. Sin embargo, al examinar la temporalidad y la ubicación geográfica de los individuos en cada par, aparecen inconsistencias: dadas las distancias temporales y espaciales entre ellos, es imposible que estén realmente emparentados.

> La mayor dificultad en un análisis de parentesco es determinar si el puntaje utilizado para medir la distancia genética —el PMR en el caso de *BREADR*— puede interpretarse efectivamente como evidencia de parentesco o no. La distribución del PMR depende directamente de la diversidad genética esperada en la población.
> En el análisis anterior, indicamos a *BREADR* que utilizara la **mediana de todos los PMR observados** como referencia. Sin embargo, en poblaciones de tamaño reducido, como las patagónicas, se espera que los valores de PMR entre individuos sean naturalmente más bajos que en poblaciones de mayor tamaño, como las de los Andes Centrales.
> Al combinar todas las poblaciones en un único análisis, es probable que *BREADR* identifique un número exagerado de pares como “aparentados” dentro de Patagonia, simplemente porque estas poblaciones experimentaron más deriva génica y presentan menor diversidad genética.
> Para evitar estos falsos positivos, debemos calcular el parentesco utilizando únicamente la tabla de PMR del subconjunto de individuos que comparten un nivel similar de diversidad genética, es decir, aquellos pertenecientes a la misma (meta)población.

En lo que sigue, analizaremos región por región. Para cada región, generaremos la subtabla de PMR considerando únicamente los pares de individuos pertenecientes a ella, verificaremos que haya un número suficiente de pares (al menos 5) y volveremos a inferir el parentesco basándonos en el valor medio esperado del PMR estimado exclusivamente a partir de ese subconjunto.

Para ello, construiremos una lista que contendrá, para cada región, la tabla de inferencias de parentesco correspondiente.

```r
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
  tmpRelatedness<-callRelatedness(tmpRegion)
  ##we can add the Population and Region label for each individual
  tmpRelatedness<-addPops(tableRel=tmpRelatedness,metaData=meta)
  ##we keep the result for that region in an element of our list
  listRelatednessPerRegion[[region]]<-tmpRelatedness
  print(paste(sum(listRelatednessPerRegion[[region]]$relationship!="Unrelated"),"pairs of related individuals found"))
  
}
```
Vemos que el número de pares aparentados disminuye drásticamente. Por ejemplo, en el sur de Patagonia ahora identificamos 20 pares, frente a los más de 300 detectados cuando analizamos a todos los individuos de manera conjunta.

Podemos volver a comprobar si los parentescos inferidos (hasta primer grado) son coherentes con la distribución temporal y geográfica de los individuos.
```r
for(region in names(listRelatednessPerRegion)){
  print(region)
  ###Let's see in which groups we have pairs of 1st-degree or Twins/Duplicated for this region
  tmpRelatedness<-listRelatednessPerRegion[[region]]
  tmpFirst<-tmpRelatedness[ tmpRelatedness$relationship %in% c("Same_Twins","First_Degree"),]
  print(table(tmpFirst$Pop1,tmpFirst$Pop2))
}
```  
Ahora vemos que los parentescos de primer grado son coherentes.  
También podemos comprobar si los parentescos inferidos hasta segundo grado son compatibles con la distribución temporal y geográfica de los individuos.
```r
for(region in names(listRelatednessPerRegion)){
  print(region)
  ###Let's see in which groups we have pairs of 1st-degree or Twins/Duplicated for this region
  tmpRelatedness<-listRelatednessPerRegion[[region]]
  tmpSecond<-tmpRelatedness[ tmpRelatedness$relationship %in% c("Second_Degree"),]
  print(table(tmpSecond$Pop1,tmpSecond$Pop2))
}
```  
Vemos que los resultados son en general bastante coherentes; sin embargo, persisten algunas inconsistencias temporales y varios de los pares de parentesco identificados aquí no fueron reportados en los trabajos originales.  
> Esto puede deberse a que utilizamos un método distinto al de los estudios previos y a que generamos nuevamente los datos pseudo-haploides (es decir, realizando un nuevo muestreo aleatorio de las lecturas).
Es fundamental no interpretar estos resultados de manera estrictamente literal.  
Se recomienda visualizar las salidas de *BREADR* para evaluar los márgenes de error asociados a los valores de PMR.

#### Inspección de los resultados de parentesco
El paquete *BREADR* ofrece distintas funciones para visualizar y evaluar los resultados:

1. **plotLOAF**: Grafica todos los valores observados de PMR (ordenados de menor a mayor) y muestra las clasificaciones de máxima probabilidad posterior mediante colores y formas.
2. **plotSLICE**: Permite representar la información diagnóstica utilizada para clasificar un par específico de individuos (definido por el número de fila o por el nombre del par).

A continuación, examinaremos si existen pares que podrían ser falsos positivos. Para ello, empezamos revisando los valores de PMR más bajos. Veremos que los 20 pares identificados como aparentados no corresponden necesariamente a los pares con PMR más bajo.

```r
### generate table from the result list per  region (and remove  columns 13 to 18 that we added with pop labels)
relatednessPatagonia=listRelatednessPerRegion[["SouthPatagonia"]][,c(1:12)]
###sort the table according to pmr
relatednessPatagonia<-relatednessPatagonia[ order(relatednessPatagonia$pmr),]
### get the ranks of the related pairs 
which(relatednessPatagonia$relationship!="Unrelated")
###to see all the related pairs we need to plot upong rank 22 then
plotLOAF(relatednessPatagonia,N=22)
```

Observamos que el par **AM66 – AM71** presenta un error estándar del PMR muy elevado y, de hecho, su PMR no es significativamente diferente de los niveles esperados para parentesco de **1er grado** ni de **2º grado**.  
Aun así, *BREADR* clasifica este par como de **1er grado**.  
Esta asignación se debe al número limitado de SNPs disponibles en esta comparación.

Analicemos este par con mayor detalle.
 ```r
 ### call rate para ambos:
 datosGeno$other$ind.metrics[datosGeno$other$ind.metrics$id %in% c("AM66","AM71"),]
 ### plotSlice
 plotSLICE(relatednessPatagonia,"AM66 - AM71")
 ```
Para estos dos individuos, observamos que las distribuciones de los valores de PMR esperados para distinguir entre distintos grados de parentesco se solapan de manera considerable, lo que dificulta discriminar correctamente el grado de parentesco para este par de individuos.  
Del mismo modo, al examinar los resultados de PMR para pares clasificados como de **2º grado** (gráfico generado con `plotLOAF`), vemos que varios pares presentan valores por encima del nivel esperado para este grado, y que la significancia estadística entre “2º grado” y “no aparentados” es muy similar. Esto sugiere la posibilidad de **falsos positivos**.

> Todo esto refuerza la importancia de no aceptar ciegamente el resultado de un método sin examinar las métricas subyacentes.  
> También es recomendable evaluar la consistencia de las conclusiones utilizando diferentes métodos y distintos paneles de SNPs.  
> En el caso de poblaciones americanas, una combinación útil es emplear **BREADR** (o **READv2**, muy similar pero incorporando el 3er grado, que debe interpretarse con mucha cautela) junto con **KIN** (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02847-7), que requiere mayor cobertura pero permite distinguir entre hermanos/as y padres-hijos en relaciones de 1er grado.

### Identificar qué individuos remover para evitar pares de aparentados

En lo que sigue, vamos a identificar qué individuos conviene sacar para eliminar los pares de individuos aparentados hasta el **primer grado**. Por cuestiones de tiempo, consideraremos que los pares identificados por **BREADR** son correctos, sin realizar el proceso completo de validación que mostramos anteriormente.

Primero, generamos una tabla que contenga todos estos pares, concatenando los resultados por región. Luego resumimos esta información en una tabla con **una fila por individuo** y cuatro columnas:

1. **id**: identificador del individuo  
2. **NumKin**: número de individuos con parentesco detectado  
3. **ListKin**: lista de dichos individuos (separados por “|”)  
4. **Call.Rate**: tasa de llamado (que recuperamos de `datosGeno$other$ind.metrics`)

> En general es preferible remover también los aparentados de **segundo grado**. Sin embargo, en este ejercicio observamos que muchos parentescos pueden ser falsos positivos, por lo que aquí nos limitamos únicamente a los parentescos de primer grado. En un estudio real, sería necesario validar o invalidar exhaustivamente todos los pares hasta segundo grado antes de generar una lista final libre de aparentados, evitando los falsos positivos.

```r
###keep all pairs
Pairs<-c()
for(region in names(listRelatednessPerRegion)){
  tmp<-listRelatednessPerRegion[[region]]
  tmp<-tmp[ tmp$relationship == "First_Degree",]
  #### next command would be f you want to remove up to 2nd-degree
  ##tmp<-tmp[ tmp$relationship %in% c("Same_Twins","First_Degree","Second_Degree"),]
  
  Pairs<-rbind(Pairs,tmp)
}

###sum up
sumUpPairs<-data.frame(matrix(NA,0,4))
names(sumUpPairs)<-c("id","NumKin","ListKin","Call.Rate")
for(line in c(1:nrow(Pairs))){
  Ind1=Pairs$Ind1[line];
  Ind2=Pairs$Ind2[line];
  if(Ind1 %in% sumUpPairs$id){
    ###case that Ind1 already found in a pair --> Update related line 
    sumUpPairs$NumKin[sumUpPairs$id==Ind1]=as.numeric(sumUpPairs$NumKin[sumUpPairs$id==Ind1]) + 1
    sumUpPairs$ListKin[sumUpPairs$id==Ind1]=paste(sumUpPairs$ListKin[sumUpPairs$id==Ind1],Ind2,sep="|")
  }else{
    ###case that Ind1 first find found in a pair --> add line for this individual
    sumUpPairs=rbind(sumUpPairs,cbind("id"=Ind1,
                                      "NumKin"=1,
                                      "ListKin"=Ind2,
                                      "Call.rate"=datosGeno$other$ind.metrics$Call.rate[metricsInds_forKeptSNPs$id==Ind1]))
  }
  
  if(Ind2 %in% sumUpPairs$id){
    ###case that Ind2 already found in a pair --> Update related line 
    sumUpPairs$NumKin[sumUpPairs$id==Ind2]=as.numeric(sumUpPairs$NumKin[sumUpPairs$id==Ind2]) + 1
    sumUpPairs$ListKin[sumUpPairs$id==Ind2]=paste(sumUpPairs$ListKin[sumUpPairs$id==Ind2],Ind1,sep="|")
  }else{
    ###case that Ind2 first find found in a pair --> add line for this individual
    sumUpPairs=rbind(sumUpPairs,cbind("id"=Ind2,
                                      "NumKin"=1,
                                      "ListKin"=Ind1,
                                      "Call.rate"=metricsInds_forKeptSNPs$Call.rate[metricsInds_forKeptSNPs$id==Ind2]))
  }

}

sumUpPairs$NumKin<-as.numeric(sumUpPairs$NumKin)
sumUpPairs$Call.rate<-as.numeric(sumUpPairs$Call.rate)
```  

Esta tabla permitirá luego eliminar el número mínimo de individuos para evitar cualquier relación hasta el 1er grado, dando prioridad a la eliminación de los individuos involucrados en el mayor número de relaciones por pares y, en caso de empates, a aquellos con menor cantidad de datos (`Call.rate`).
```r
listToRemove<-data.frame(matrix(NA,0,3))
names(listToRemove)<-c("id","NumKin","Call.rate")
sumUpPairs_BU<-sumUpPairs

###until we have related pairs 
while(nrow(sumUpPairs)>0){
    ##get pairs with maximum number of remaining kinship relation
	  tmp<-sumUpPairs[ sumUpPairs$NumKin  == max(sumUpPairs$NumKin),]
	  ##order those pairs according to call.rate
	  tmp<-tmp[order(tmp$Call.rate),]
	  
	  ##we will remove the individual with lowest call.rate
	  id<-tmp$id[1]
	  
	  ##we keep variables that will be saved into the listToRemove table
	  Call.rate<-tmp$Call.rate[1]
	  NumKin<-sumUpPairs_BU$NumKin[sumUpPairs_BU$id==id]
	  listToRemove<-rbind(listToRemove,cbind(id,NumKin,Call.rate))
	  
	  ##we get the loist of individuals related to that removed individual (to update the sumUp table)
	  listAssoc<-strsplit(tmp$ListKin[1],split="\\|")[[1]]
	  
	  ###update the sumUp table:1st remove the line for the removed individual
	  sumUpPairs<-sumUpPairs[ sumUpPairs$id != id,]
	  
	  ###then update the lines for related individuals (deleting the id of the removed individual in the list and uodating the number of kinship relation)
	  for(i in which(sumUpPairs$id %in% listAssoc)){
	    tmp<-strsplit(sumUpPairs$List[i],split="\\|")[[1]]
	    tmp<-tmp[ tmp!=id]
	    sumUpPairs$ListKin[i]<-paste(tmp,collapse="|")
	    sumUpPairs$NumKin[i]<-sumUpPairs$NumKin[i]-1
	  }
	  sumUpPairs<-sumUpPairs[ sumUpPairs$NumKin!=0,]
}

```

## Generar conjunto de datos filtrado

Ahora que ya identificamos las variantes a filtrar (por altos niveles de datos faltantes y/o por ser monomórficas), los individuos con muchos datos faltantes, y la lista de individuos que debemos remover para evitar pares aparentados, podemos generar los ficheros filtrados de **PLINK**.

## Filtrado

Vamos a filtrar los datos eliminando variantes e individuos según su proporción de **datos faltantes**. También eliminaremos las posiciones cuya **frecuencia del alelo menor (MAF)** sea inferior al 1%.

Como primer paso, eliminaremos los SNPs con datos disponibles para **solo un individuo** (es decir, con una tasa de llamado inferior a 2/66).


```r
###filter  SNPs for call.rate (we force to recalculate the metrics and remove monomorhisms)
nLoc_init=datosGeno$n.loc
nInd_init=length(datosGeno$ind.names)
datosGeno_filter<-gl.filter.callrate(datosGeno, threshold=2/66,recalc=TRUE,mono.rm=TRUE) 
nLoc_callRate=datosGeno_filter$n.loc
print(paste("N SNPs before call.rate:",nLoc_init,"; N SNPs after call.rate and monomorphism filtering:",nLoc_callRate))
```
Pueden aparecer uno <em>warnings meassge</em> pero verán que solo se tratan de mensajes asociados a los histogramas que se generan.

Luego, vamos a eliminar los individuos que conservan **más de 10,000 variantes** tras este filtrado inicial y los que identificamos para evitar parentesco en la muestra.
```r
metricsInds_forKeptSNPs <- datosGeno_filter$other$ind.metrics
###add a column of the number of called snps (Call.rate*nLoc_callrate)
metricsInds_forKeptSNPs$CalledSNPs=metricsInds_forKeptSNPs$Call.rate*nLoc_callRate
###save in a vector individuals with less than 10000 SNPs
INDs_to_exclude=metricsInds_forKeptSNPs$id[metricsInds_forKeptSNPs$CalledSNPs< 10000]
print(paste("excluding ",length(INDs_to_exclude),"individuals with less than 10000 SNPs"))
### some may have been already flagged in our kinship analyses
print(paste("of which  ",sum(INDs_to_exclude %in% listToRemove$id)," already flagged to be remived from kinship analyses"))
### update the list to exlude with list from kinship
INDs_to_exclude<-unique(c(INDs_to_exclude,listToRemove$id))
### we can now remove those individuals (and we ask to recalculate)
datosGeno_filter<-gl.drop.ind(x=datosGeno_filter,ind.list=INDs_to_exclude,recalc = TRUE,mono.rm=TRUE)
##
nInd_final=length(datosGeno_filter$ind.names)
nLoc_final=datosGeno_filter$n.loc
print(paste("filtered dataset contains:",nInd_final,"inds and",nLoc_final,"snps"))

```

Ahora podemos escribir los datos en la carpeta del Modulo3, ya que los usaremos para esto.
En los paquetes que ocupamos no hay una función que permita escribir al formato <em>EIGENSTRAT</em> directamente, pero desde el objeto `genLight` que creamos es sencillo.

```r
prefOUT="Modulo3/filteredDataSet"
### Convert genlight genotype data to a numeric matrix
geno_mat <- as.matrix(datosGeno_filter)
# Replace missing data with 9 (per convention)
geno_mat[is.na(geno_mat)] <- 9
###check size
dim(geno_mat)
##write the genotype data. Watchout you have to transpose the matrix
write.table(t(geno_mat), paste(prefOUT,".geno",sep=""),col.names=F,row.names=F,quote=F,sep="")
##write the individual annotation file
indData<-datosGeno_filter$other$ind.metrics[,c("id","Sex","Family")]
indData$Sex[is.na(indData$Sex)]<-"U"
write.table(indData, paste(prefOUT,".ind",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
##write the snp annotation file
snpData<-datosGeno_filter$other$loc.metrics[,c("AlleleID","chromosome","cM","position","allele.1","allele.2")]
write.table(snpData, paste(prefOUT,".snp",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
```
