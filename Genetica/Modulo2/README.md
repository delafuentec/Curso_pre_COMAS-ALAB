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

**Archivo `.snp.txt`:**<br?:

snp1 1 0.1001 100000 A C<br?
snp2 1 0.6162 600000 T G<br?
snp3 2 0.5125 513341 C G<br?
snp4 4 0.1512 251334 G A<br?


**Archivo `.ind.txt`:**<br?

IndA M Pop1<br?
IndB M Pop1<br?
IndC F Pop2<br?
IndD M Pop3<br?


**Archivo `.geno.txt`:**<br?

0192<br?
1122<br?
2222<br?
0211<br?


Podemos interpretar que el individuo **IndA**, perteneciente a la **población 1**, presenta los siguientes genotipos:
- `snp1`: C/C  
- `snp2`: A/C  
- `snp3`: faltante  
- `snp4`: A/A



## Lectura de los datos

Existen varios paquetes de **R** que permiten leer y procesar estos formatos de genotipos.  
En este curso utilizaremos el paquete *dartR.base*, que incluye la función `gl.read.PLINK()`, la cual permite importar directamente archivos en formato **PLINK** en un objeto [`genlight`](https://rdrr.io/cran/adegenet/man/genlight.html).  

Además, esta función admite incorporar un archivo de metadatos por individuo (por ejemplo, `finalSet.metadataPerind.txt`), que contiene información complementaria como la población, la región o el estudio de origen.

```
require(dartRverse)
require(dartR.base)
require(ade4)

### you maybe need to change the directory below
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB/Genetica/")
datosGeno<-gl.read.PLINK("PatagoniaDataSetWithOutgroups_ALAB2025/test.plink")
```

Al leer los datos con la opción `verbose = TRUE`, se muestran varios mensajes informativos (por ejemplo, sobre loci monomórficos o con datos faltantes).  

En algunas versiones de *dartR.base*, puede aparecer el siguiente mensaje:  
*“The slot loc.all, which stores allele name for each locus, is empty. Creating a dummy variable (A/C) to insert in this slot.”*  

Este mensaje indica que la información sobre los alelos de cada locus está vacía, y el programa crea una variable ficticia (“A/C”) para completarla.  
Para corregirlo de forma manual, puede ejecutarse el siguiente comando:

```
### fix allele names
bim <- read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.plink.bim", header = FALSE)
# Columns: CHR SNP CM POS A1 A2
datosGeno@loc.all <- paste(bim$V5, bim$V6, sep="/")
```

A continuación, vamos a asignar los nombres de las poblaciones a los individuos.   Esta información se encuentra almacenada en la columna **"Family"** del elemento `ind.metrics` de `datosGeno`. Vamos a transferir esos nombres para que queden correctamente asociados a la variable `pop`.
``` 
datosGeno$pop<-as.factor(datosGeno$other$ind.metrics$Family)
```


## Exploración de los datos

Podemos examinar la calidad general del conjunto de datos observando, por ejemplo:  
- la **proporción de individuos con datos disponibles por locus**, y  
- la **frecuencia alélica** de cada variante.  

Para ello, podemos ejecutar las siguientes funciones:
```
hist(datosGeno@other$loc.metrics$CallRate,main="Call Rate per Locus")
hist(datosGeno@other$loc.metrics$maf,main="Minor Allele Frequency per Locus",n=100)
```
Observamos que existen loci con una **tasa de llamado** (*call rate*) baja (< 0.2) y posiciones **monomórficas** (MAF = 0).  

También podemos explorar la **proporción de loci con datos disponibles por individuo**, para identificar muestras con una alta proporción de datos faltantes.
 ```
hist(datosGeno@other$ind.metric$Call.rate,main="Call Rate per individual",n=10)
```
En cuanto a la **tasa de heterocigosidad**, observamos que todos los individuos presentan genotipos **homocigotos**.
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

### Introducción al análisis de parentesco
Existen varios métodos para detectar individuos aparentados. Vamos a usar [<em>BREADR</em>](https://joss.theoj.org/papers/10.21105/joss.07916). 
Este método es una extensión del método quizás más usado en el campo ([<em>READv2</em>](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03350-3)), ya que se base en el escore PMR.
Este escore es la proporción de sitios superpuestos con genotipos no coincidentes, es decir la fracción de posiciones del genoma en las que dos individuos tienen datos genotípicos válidos (no faltantes) y sus llamadas de genotipo difieren. En otras palabras, mide qué porcentaje de los sitios comparables entre ambos individuos presentan genotipos distintos, reflejando el grado de discrepancia genética entre ellos.
Las diferencias más importantes con READv2 son:
1. el uso de sitios independientes (separados por recombinación durante la meiosis). Por defecto, el programa usa un locus cada 10,000 par de bases (paramétro ` filter_length`).
2. su implementación en R.

Noten que estos métodos implementan aproximaciones que tomen en cuenta las limitaciones de los datos en ADN antiguo (pseudo-haploides y con baja cobertura). Para analizar datos modernos que no padecen de estas limitaciones, se recomienda usar otros métodos, por ejemplo [<em>KING</em>](https://www.kingrelatedness.com/).
<em>BREADR</em>  usa el formato <em>eigenstrat</em>. Vamos a analizar solamente los individuos antiguos.

### Inferencia de los grados de parentesco
#### Generación de la métrica necesaria (PMR)
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

#### Definir el grado de parentesco
Ahora vamos a ver cuales son los individuos supuestamente aparentados, y el grado de parentesco. Para definir el parentesco entre 2 individuos, se calcula los niveles de PMR esperados en la población, según la distribución de los PMR en todos los pares de individuos analizados.
Para esto se asume que el conjunto de datos muestreados está compuesto principalmente por pares no emparentados (hasta segundo grado), entonces la mediana del PMR,  denotado $\bar{p}$, será una estimación confiable.
Tomando como referencia los planteamientos de <em>READv2</em>, definimos ahora el valor medio esperado del PMR para una relación de grado <em>k</em> = 0, 1, 2 como: <br>
$p_k = \bar{p} \left(1 - \frac{1}{k + 1}\right)$.<br>
Una vez definidos este valor medio esperado del PMR para cada relación de grado, podemos poner a prueba si el PMR para un par de individuos dados es significativamente más pequeño que cada nivel esperado <em>k</em>. El parentesco será entonces el <em>k</em> más bajo que no da significativo, tal como se observa en la figura abajo. <br>
<img width="729" height="344" alt="image" src="https://github.com/user-attachments/assets/d11688e8-3d94-48aa-bb11-d3e9cbef8cb9" />

Vamos ahora calcular el parentesco entre todos los pares de individuos de la muestra, analizados en su conjunto.

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
Veamos cuantos pares de individuos aparentados, en qué región y población.

```
###Number of related pairs indivuduals:
related_allAncientTogether<-relatedness_allAncientTogether[ relatedness_allAncientTogether$relationship!="Unrelated",]
nrow(related_allAncientTogether)
table(related_allAncientTogether$relationship)

### we see a lot of related pairs. Let's see in which region
table(related_allAncientTogether$Region1,related_allAncientTogether$Region2)

###Let's see in which groups we have pairs of 1st-degree or Twins/Duplicated
first_allAncientTogether<-related_allAncientTogether[ related_allAncientTogether$relationship %in% c("Same_Twins","First_Degree"),]
table(first_allAncientTogether$Pop1,first_allAncientTogether$Pop2)

### We can already see some temporal/spatial inconsistencies 
```

Vemos entonces que detectamos muchos pares de aparentados en el Sur de Patagonia. Además mirando la temporalidad y la geografía de cada individuo en un par de individuos aprentados, vemos incosnistencias ya que con esta separación temporal y espacial imposible que sean aparentados.<br>
La mayor dificulta en un análisis de parentesco es definir si el escore usado para medir la distancia genética (PMR en el caso de <em>BREADR</em>) se puede interpretar como parentesco. De hecho, la distribución del PMRdepende de la diversidad genética esperada en la población. <br>
Lo que acabmos de hacer es decir a <em>BREADR</em> de usar la mediana de los PMR observados para todos los pares de individuos analizados. Sin embargo, en poblaciones con tamaño poblacional reducido, como es el caso de las patagónicas, se espera un escore PMR entre individuos mucho menor que en poblaciones de tamaño más grande, como las de Centro Andes. <br>
Entonces, si comparamos todas las poblaciones a la vez, es probable que se detecte muchos pares de individuos aparentados en poblaciones de Patagonia, que son poblaciones que evolucionaron con más deriva génica, por ende que tienen menor diversidad genética.
Para evitar estos falsos positivos podemos pasar al método solo la tabla de escore PMR para el subconjunto de individuos para los cuales esperamos niveles de diversidad genética similares por pertenecer a la misma (meta)población.<br>
En lo que sigue, vamos a ir, región por región, generar la sub-tabla de PMR por los pares de individuos de esta región, verificar que hay suficientes pares (digamos por lo menos 5) y proceder de nuevo a las inferencias de parentesco desde valor medio esperado del PMR para cada relación de grado según estos pares únicamente.<br>
Para esto, vamos a generar una lista que contienen la tabla de inferencia de parentesco por región.

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
  tmpRelatedness<-callRelatedness(tmpRegion)
  print(paste(sum(listRelatednessPerRegion[[region]]$relationship!="Unrelated"),"pairs of related individuals found"))
  ##we can add the Population and Region label for each individual
  tmpRelatedness<-addPops(tableRel=tmpRelatedness,metaData=meta)
  ##we keep the result for that region in an element of our list
  listRelatednessPerRegion[[region]]<-tmpRelatedness
  
  
}
```
Vemos que reducimos drastícamente el número de pares de aparentados encontrados. Por ejemplo, ahora encontramos 20 en el Sur de Patagonia, contra 321 en el analísis con todos los individuos a la vez.
Podemos verificar de nuevo si los resultados de parentesco (hasta el 1er grado) son concordantes con la distribución temporal y geografíca.
```
for(region in names(listRelatednessPerRegion)){
  print(region)
  ###Let's see in which groups we have pairs of 1st-degree or Twins/Duplicated for this region
  tmpRelatedness<-listRelatednessPerRegion[[region]]
  tmpFirst<-tmpRelatedness[ tmpRelatedness$relationship %in% c("Same_Twins","First_Degree"),]
  print(table(tmpFirst$Pop1,tmpFirst$Pop2))
}
```  
Ahora vemos que los parentescos al 1er grado son coherentes.
Podemos verificar también si los resultados de parentesco (hasta el 2ndo grado) son concordantes con la distribución temporal y geografíca.
```
for(region in names(listRelatednessPerRegion)){
  print(region)
  ###Let's see in which groups we have pairs of 1st-degree or Twins/Duplicated for this region
  tmpRelatedness<-listRelatednessPerRegion[[region]]
  tmpSecond<-tmpRelatedness[ tmpRelatedness$relationship %in% c("Second_Degree"),]
  print(table(tmpSecond$Pop1,tmpSecond$Pop2))
}
```  
Vemos que es bastante coherente también, sin embargo muchos pares de parentesco que encontramos en este ejercicio no fueron reportados en los trabajos originales. Esto se debe que usamos otro método, un conjunto de variantes diferentes y se generaron los datos pseudo-haploides analizados nuevamente (es decir rehaciendo un muestreo aleatorio de las lecturas).
Vemos entonces la necesidad de no confiar en nuestros resultados ciegamente. Se recomienda, visualizar los resultados de <em>BREADR</em>, para apreciar los margenes de error asociados a los PMR.<br>

#### Inspección de los resultados de parentesco
El paquete <em>BREADR</em> provee diferentes funciones para esto:
1. plotLOAF: Grafica todos los valores observados de PMR (ordenados de menor a mayor), con las clasificaciones de máxima probabilidad posterior representadas mediante color y forma.
2. plotSLICE: Una función para graficar la información diagnóstica al clasificar un par específico de individuos (definido por el número de fila o el nombre del par).

Veamos entonces si hay pares que pueden ser potencialmente falso positivos. Primero veamos los PMR para los pares con PMR más bajo (los 20 pares de aprentados no son necesariamente los con PMR más bajo como veremos).

```
### generate table from the result list per  region (and remove  columns 13 to 18 that we added with pop labels)
relatednessPatagonia=listRelatednessPerRegion[["SouthPatagonia"]][,c(1:12)]
###sort the table according to pmr
relatednessPatagonia<-relatednessPatagonia[ order(relatednessPatagonia$pmr),]
### get the ranks of the related pairs 
which(relatednessPatagonia$relationship!="Unrelated")
###to see all the related pairs we need to plot upong rank 22 then
plotLOAF(relatednessPatagonia,N=22)
```

Vemos entonces que el par AM66 - AM71 tiene un error estandar asociado al PMR muy grande, y de hecho no es significativamente diferente del nivel esperado para el 1er grado ni el 2o grado.
Sin embargo, <BREADR> da como resultado "1er grado". Esto se debe al número de SNPs limitado en esta comparación. Veamos más en detalle este par.
 ```
 ### call rate para ambos:
 datosGeno$other$ind.metrics[datosGeno$other$ind.metrics$id %in% c("AM66","AM71"),]
 ### plotSlice
 plotSLICE(relatednessPatagonia,"AM66 - AM71")
 ```
Vemos entonces que, a lo mejor, sería más prudente definir este par aparentado  al 2o grado y no al 1er grado.<br> 
Lo mismo, cuando miramos los resultados de PMR para pares de aparentados al 2o grado, vemos que varios pares están por encima del valor esperado para 2o grado, y que la significancia de 2o grado o No aparentados son similares. Entonces podrían ser falsos positivos.
Todo esto, para insistir con la idea de no usar un resultado que retorna un método sin inspeccionar las métricas subyacentes.<br>
Se recomienda también explorar la consistencia de las conclusiones con diferentes métodos y paneles de SNPs. Personalmente, me gusta usar para Ámerica los métodos <em>BREADR</em> (o <em>READv2</em>, muy similar a <em>BREADR</em> pero incorpora grado hasta el 3er grado, usar con mucha precaución hasta este grado) y [<em>KIN</em>](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02847-7) (que requiere mayor cobertura pero permite distinguir entre relación al primer gardo entre hermanos/as y padres/hijos).


## Generación del conjunto de datos filtrado




vamos a generar una tabla que resume las relaciones al primer grado (en aDNA antiguo se suelen guardar los individuos aparentados al segundo grado). 

```
tableSUM<-c()
for(region in names(listRelatednessPerRegion)){
  ###for each region, keep only lines for 1st degree relationship 
  tmp1st<-listRelatednessPerRegion[[region]]
}


## Parentesco
