# Analisis de diversidad y estructura genética

En este módulo realizaremos **dos análisis clásicos en genética de poblaciones humanas**:

1. **Análisis en Componentes Principales (ACP)**  
2. **Estimación de proporciones de componentes genéticos**

> Tradicionalmente, estos análisis se realizan con los programas **smartpca** (del paquete  
> [EIGENSOFT](https://github.com/DReichLab/EIG)) y **ADMIXTURE**  
> (https://dalexander.github.io/admixture/).  
> Sin embargo, el paquete **LEA** de R implementa métodos muy similares —con la ventaja de que  
> el algoritmo utilizado para estimar proporciones de componentes genéticos es **más robusto,  
> menos sensible al ruido y mucho más rápido** que el de *ADMIXTURE*.  

A lo largo del módulo exploraremos cómo ejecutar estos análisis, cómo interpretar sus resultados  
y qué aspectos considerar al trabajar con datos genómicos, especialmente en contextos con  
coberturas desiguales o datos de ADN antiguo.

## Lectura de los datos filtrados

En este paso vamos a cargar los **datos filtrados** que generamos en el módulo anterior.  
Estos archivos ya contienen:

- solo los individuos seleccionados,  
- las variantes con suficiente información,  
- y el formato adecuado para trabajar con las funciones del paquete **LEA**.

A partir de estos datos realizaremos los análisis de ACP y de estimación de componentes genéticos.

```r
require(LEA)
### you maybe need to change the directory below
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB/Genetica/")
pref="Modulo3/filteredDataSet"

```

## Análisis en Componentes Principales (ACP)

El ACP permite **reducir la dimensionalidad** de un conjunto de datos genéticos muy grande, resumiendo la variación en unos pocos componentes principales que capturan la mayor parte de la información.

### Primer ACP a escala continental

Como primer paso, realizaremos un ACP incluyendo **todos los individuos antiguos** presentes en nuestro conjunto de datos filtrado (Módulo anterior).  
Esto nos permitirá visualizar las principales estructuras de variación genética a gran escala.

Antes de correr el análisis, podemos explorar qué individuos están presentes en el dataset leyendo el archivo `.ind` y vinculándolo con nuestros metadatos.
```r
ind<-read.table(paste(pref,".ind",sep=""),stringsAsFactors=F,header=F)
names(ind)<-c("id","sex","pop")
meta=read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.metadataPerind.txt",stringsAsFactors = F,header=T)
### generate table of the number of individual per region
table(meta$Region[ meta$id %in% ind$id])
```

#### Correr el ACP
```r
pc = pca(paste(pref,".geno",sep=""), scale = TRUE)
```
Los resultados del ACP se guardarán automáticamente en una carpeta llamada `<pref>.pca`,  
donde `<pref>` corresponde al prefijo que usamos para los archivos de entrada.
Primero vamos a examinar el porcentaje de varianza explicada por cada componente principal.

```r
### we can read the eigenvalues file
eigenvalues<-read.table(paste(pref,".pca/filteredDataSet.eigenvalues",sep=""),stringsAsFactors=F,header=F)
### let see the porcentage of variance explained per principal components
plot(eigenvalues$V1/sum(eigenvalues$V1)*100, lwd=5, col="mediumseagreen",xlab=("PCs"),ylab="% variance explained")
```
Al visualizar los *eigenvalues*, observamos que la proporción de varianza explicada disminuye conforme aumenta el número de componentes principales (un patrón inherente al método).  
Sin embargo, el patrón no sigue la forma clásica esperada en un ACP: normalmente se observa un “codo” (*elbow*), donde los primeros componentes capturan una gran parte de la variación y, a partir de cierto punto, la varianza explicada por los componentes restantes cae a un nivel casi constante. Más allá de ese “piso”, cada componente adicional aporta muy poca información biológica.

En nuestro caso, este patrón no aparece claramente. Esto se debe a que estamos trabajando con datos de **ADN antiguo**, que suelen presentar:
- **Altos niveles de datos faltantes**, distribuidos de forma heterogénea entre individuos.  
- **Coberturas muy variables**, que afectan la estimación de las covarianzas entre individuos.  
- **Menor señal genética por individuo**, lo que tiende a “aplanar” la distribución de varianza explicada.  

Como consecuencia, los valores propios (*eigenvalues*) no muestran un “codo” definido y la disminución de varianza explicada es más gradual y ruidosa que en un ACP de genomas modernos o completos.

Ahora vamos a graficar las combinaciones de componentes principales de a pares: PC2 vs. PC1, PC4 vs. PC3, …, PC10 vs. PC9.

```r
### we can read the projections of the individuals: one line per individual.
### the individuals are ordered as in the input file (thus reading <pref>.ind.txt we know which column corresponds to which individual)
projections<-read.table(paste(pref,".pca/",strsplit(pref,split="/")[[1]][2],".projections",sep=""),stringsAsFactors=F,header=F)
## add info per ind
projections$id=ind$id
projections$pop=ind$pop
###get metafile for plotting
projections<-merge(projections,meta,by=c("id","pop"))

###let's plot the first 10 PCs
pdf(paste(pref,".PCAWithAll.pdf",sep=""),height=10)
par(mfrow=c(3,2))
forLeg<-unique(projections[,c("Region","pop","Color","Point")])
forLeg<-forLeg[ order(forLeg$Region,forLeg$pop),]
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=forLeg$Point,pt.bg=forLeg$Color,col=ifelse(forLeg$Point<21,forLeg$Color,"black"),
        legend=paste(forLeg$Region,forLeg$pop),ncol=1,cex=0.5,pt.lwd=0.5)
for(i in seq(1,9,2)){
  plot(projections[,paste("V",i,sep="")],projections[,paste("V",i+1,sep="")],
          pch=projections$Point,
          bg=projections$Color,
          col=ifelse(projections$Point<21,projections$Color,"black"),
          xlab=paste("PC",i),
          ylab=paste("PC",i+1))
}
dev.off()
```
Los gráficos están en un pdf: `<pref>.PCAWithAll.pdf`
> Atencion: En smartpca de EIGENSOFT: El archivo llamado eigenvectors (*.evec) no contiene los eigenvectores de los SNPs, sino las proyecciones de los individuos sobre los componentes principales.<br>
> Cada fila = un individuo<br>
> Cada columna = PC1, PC2, …<br>
> Las varianzas de los PCs están en el archivo eigenvalues (*.eval).<br>
> Si quieres conocer la contribución de cada SNP (loadings), smartpca no lo entrega por defecto; se requieren métodos adicionales (lsqproject u otros).<br>
> **Consejo práctico**: para graficar la estructura poblacional, usar el archivo eigenvectors de smartpca como proyecciones de individuos, no como cargas de SNPs.
#### Interpretación de los resultados del primer ACP

El análisis de componentes principales (ACP) nos permite visualizar cómo se relacionan genéticamente los individuos antiguos entre sí. En este caso, los primeros componentes capturan patrones muy coherentes con lo que se conoce sobre el poblamiento de América.

En primer lugar, el ACP separa de manera marcada al individuo Upward Sun River 1 (ASR1), procedente de Beringia. Esto es esperable, ya que este individuo representa un linaje muy basal, es decir, cercano a las poblaciones ancestrales que dieron origen a todos los grupos americanos. Por lo tanto, es lógico que aparezca como el más diferente en el espacio genético.

En segundo lugar, los componentes siguientes distinguen a varios individuos muy tempranos —como Los Rieles, Spirit Cave, Anzick y uno de Lagoa Santa— del resto de Sudamérica. Estos genomas pertenecen a momentos iniciales del poblamiento y conservan señales genéticas que se perdieron o se transformaron en las poblaciones sudamericanas posteriores. La separación que vemos en el ACP refleja justamente esa diversidad temprana: distintos grupos que fueron llegando a América y que aún no habían convergido en la estructura genética que predomina más tarde.

Además, algunos componentes capturan diferencias más sutiles dentro de estos mismos linajes tempranos. Por ejemplo, Anzick y Spirit Cave aparecen separados entre sí, lo que indica que ya existía variación estructurada en Norteamérica desde el inicio del Holoceno. Algo similar ocurre con Los Rieles y el individuo temprano de Lagoa Santa, que muestran distancias genéticas entre sí pese a ser comparables en edad.

Finalmente, uno de los componentes diferencia claramente al individuo Ayayema, del extremo sur de Patagonia. Su posición única en el ACP sugiere que este individuo representa una variación genética específica de esa región austral, probablemente resultado de aislamiento geográfico y dinámicas demográficas propias de Patagonia.

En conjunto, estos patrones nos muestran que la diversidad genética temprana en América fue más compleja de lo que suele imaginarse: no existió un único “grupo fundador”, sino varios linajes tempranos que se diferenciaron entre sí y que dejaron huellas distintas en las poblaciones tardías. El ACP resume esta historia en un espacio bidimensional, permitiendo visualizar cómo se distribuye esa diversidad a través del tiempo y del continente.

Vimos que en varios componentes un individuo de Lagoa Santa es diferente de los demas del mismo sitio.
Primero podemos verificar que se trate siempre del mismo
```r
### get projections for sumidouro individuals on PCs 2,4,5 and 7
sumi<-projections[grepl("Brazil_Sumidouro_",projections$pop),c("id","pop",paste("V",c(2,4,5,7),sep=""))]
print(sumi)
```

Vemos entonces que es Sumidouro5 (~10400 AP) que se diferencia de los otros. Sumidouro 5, un individuo del Holoceno temprano proveniente de la cueva de Sumidouro (Lagoa Santa, Brasil), destaca por varias razones biológicas y arqueológicas. En primer lugar, su morfología craneal fue reconocida desde mucho tiempo como particularmente distintiva: comparaciones multivariantes indicaron afinidades con poblaciones africanas y austro-melanesias, en lugar de con los indigenas  americanos tardíos o asiáticos [(Neves et al. 2007)](bpb-us-w2.wpmucdn.com). <br>
Más recientemente, datos genómicos han revelado en Sumidouro 5 una señal de ascendencia mitocondrial D4h3a y un perfil autosómico que no se corresponde completamente con otros individuos tempranos de Lagoa Santa, lo que sugiere conexiones genéticas más amplias (incluso con poblaciones amazónicas ancestrales), ([Moreno-Mayar et al. 2018](https://www.science.org/doi/10.1126/science.aav2621), [Ferraz et al. 2023](https://doi.org/10.1038/s41559-023-02114-9)). <br>
Esta combinación de rasgos morfológicos arcaicos y una ascendencia genética poco común indica que Sumidouro 5 pudo pertenecer a un linaje parcialmente distinto a otros pobladores tempranos de Lagoa Santa, aportando evidencia de una mayor diversidad biológica en los primeros habitantes del este de Sudamérica y de posibles rutas de migración complejas durante el Holoceno temprano.

### ACP enfocado en procesos más tardíos de Sudamérica

Ahora vamos a repetir el análisis, esta vez guardando solo los individuos de Central Chile del Holoceno Medio y Tardío, de la Pampa y del Sur de Patagonia.  
El objetivo es focalizarnos en dinámicas más recientes en el Sur de Patagionia.

#### Generar sub-conjunto
Para ello, primero leeremos los datos en formato *eigenstrat* y generaremos un sub-conjunto sin estos individuos.
Como primer paso, vamos a convertir los datos en un objeto `genLight`.

```r
require(dartR.base)
### read file
geno <- read.geno(paste(pref,".geno",sep="")) 
# Read .ind file
ind <- read.table(paste(pref,".ind",sep=""), stringsAsFactors = FALSE,header=F)
names(ind)<-c("id","sex","pop")
# Read .snp file
snp <- read.table(paste(pref,".snp",sep=""), stringsAsFactors = FALSE,header=F)
names(snp)<-c("id","chromosome","cM","position","allele.1","allele.2")
# Build the genlight object
# need to say that 9 are NAs
geno[geno == 9] <- NA
gl <- new("genlight", geno)
#Add individual metadata
indNames(gl) <- ind[,1]                  # individual names
pop(gl) <- as.factor(ind[,3])           # population column
gl@other$sex <- ind[,2]                 #sex

#Add SNP metadata
locNames(gl) <- snp[,1]                  # SNP IDs
chromosome(gl) <- snp[,2]                # chromosome
position(gl) <- snp[,4]  

#Add metrics to the gl object
gl@other$loc.metrics.flags <- list(
    monomorphs = TRUE,
    oneallele = TRUE,
    nosnp = TRUE
)
gl@other$loc.metrics <- data.frame(
    loc.id = locNames(gl),
    chromosome = snp$chromosome,
    cM = snp$cM,
    position = snp$position,
    allele.1 = snp$allele.1,
    allele.2 = snp$allele.2,
    CallRate= NA,
    stringsAsFactors = FALSE
)
gl@other$ind.metrics.flags <- list(
    loc.freq = TRUE,
    het = TRUE,
    callrate = TRUE
)
gl@other$ind.metrics <- data.frame(
    id = indNames(gl),
    Sex = ind$sex,
    Family = ind$pop,
    Call.rate=NA,
    stringsAsFactors = FALSE
)

gl <- gl.recalc.metrics(gl)
```

Ahora vamos a generar la lista de individuos que queremos conservar en el subconjunto de datos.  
Para ello, filtraremos el archivo de información de individuos según los criterios definidos anteriormente.


```r
# make list of  individuals to keep
KeptPops<-meta$pop[ grepl("SouthPatagonia",meta$Region) |  grepl("Pampa",meta$Region) | grepl("CentralChile_MH",meta$Region) | grepl("CentralChile_LH",meta$Region) ] 
toKeepID<-ind$id[  ind$pop %in% KeptPops]
# check from which regions come from the remaining individuals
print(table(meta$Region[ meta$id %in% toKeepID]))
##looks ok --> generate subdataset, removing monorphisms (need to recalculate the metrics for further filtering)
```

Ahora podemos generar el subconjunto de datos.  
Como paso adicional de control de calidad, volveremos a filtrar las variantes monomórficas y aquellas con menos de 2 individuos con valores.  
Además, eliminaremos los individuos que tengan menos de **10.000 SNPs** disponibles, para asegurar que el análisis de ACP tenga suficiente información genética por individuo.

```r
gl<-gl.keep.ind(gl,ind.list =toKeepID,recalc = T,mono.rm = T)
gl<-gl.filter.callrate(gl, threshold=2/nInd(gl),recalc=T)
###check numbers
print(paste("from",nrow(ind),"individuals, we wanted to keep",length(toKeepID),", and the sub-dataset has ",nInd(gl)))
print(paste("from",nrow(snp),"snps",nLoc(gl)," remain in the sub-dataset"))
```

Ahora podemos escribir los ficheros, de la misma forma que vimos en el Módulo 2.
```r
geno_mat <- as.matrix(gl)
# Replace missing data with 9 (per convention)
geno_mat[is.na(geno_mat)] <- 9
###check size
dim(geno_mat)
##write the genotype data. Watchout you have to transpose the matrix
write.table(t(geno_mat), paste(pref,"_subset.geno",sep=""),col.names=F,row.names=F,quote=F,sep="")
##write the individual annotation file
ind_sub<-gl$other$ind.metrics[,c("id","Sex","Family")]
names(ind_sub)<-c("id","Sex","pop")
ind_sub$Sex[is.na(ind_sub$Sex)]<-"U"
write.table(ind_sub, paste(pref,"_subset.ind",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
##write the snp annotation file
snp_sub<-gl$other$loc.metrics[,c("loc.id","chromosome","cM","position","allele.1","allele.2")]
write.table(snp_sub, paste(pref,"_subset.snp",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
```

#### ACP con el sub-conjunto

Con el sub-conjunto ya filtrado y depurado, podemos proceder a realizar el ACP.  
Este análisis permitirá explorar la variación genética enfocándonos en procesos más tardíos de la historia poblacional sudamericana, sin la influencia de los individuos de América del Norte ni de los grupos muy tempranos de Chile y Brasil.

El procedimiento es idéntico al aplicado previamente: generamos el archivo en formato *eigenstrat*, ejecutamos el ACP y luego inspeccionamos tanto los valores propios como la distribución de los individuos en los distintos planos principales.
```r
pc_sub = pca(paste(pref,"_subset.geno",sep=""), scale = TRUE)
```
Los resultados del ACP se guardarán automáticamente en una carpeta llamada `<pref>_subset.pca`,  
donde `<pref>_subset` corresponde al prefijo que usamos para los archivos de entrada.
```r
### we can read the eigenvalues file
eigenvalues_sub<-read.table(paste(pref,"_subset.pca/filteredDataSet_subset.eigenvalues",sep=""),stringsAsFactors=F,header=F)
### let see the porcentage of variance explained per principal components
plot(eigenvalues_sub$V1/sum(eigenvalues_sub$V1)*100, lwd=5, col="mediumseagreen",xlab=("PCs"),ylab="% variance explained",main="For subset")

### we can read the projections of the individuals: one line per individual.
### the individuals are ordered as in the input file (thus reading <pref>.ind.txt we know which column corresponds to which individual)
projections_sub<-read.table(paste(pref,"_subset.pca/",strsplit(pref,split="/")[[1]][2],"_subset.projections",sep=""),stringsAsFactors=F,header=F)
## add info per ind
projections_sub$id=ind_sub$id
projections_sub$pop=ind_sub$pop
###get metafile for plotting
projections_sub<-merge(projections_sub,meta,by=c("id","pop"))

###let's plot the first 10 PCs
pdf(paste(pref,".PCAWithSubSet.pdf",sep=""),height=10)
par(mfrow=c(3,2))
forLeg_sub<-unique(projections_sub[,c("Region","pop","Color","Point")])
forLeg_sub<-forLeg_sub[ order(forLeg_sub$Region,forLeg_sub$pop),]
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=forLeg_sub$Point,pt.bg=forLeg_sub$Color,col=ifelse(forLeg_sub$Point<21,forLeg_sub$Color,"black"),
        legend=paste(forLeg_sub$Region,forLeg_sub$pop),ncol=1,cex=0.5,pt.lwd=0.5)
for(i in seq(1,9,2)){
  plot(projections_sub[,paste("V",i,sep="")],projections_sub[,paste("V",i+1,sep="")],
          pch=projections_sub$Point,
          bg=projections_sub$Color,
          col=ifelse(projections_sub$Point<21,projections_sub$Color,"black"),
          xlab=paste("PC",i),
          ylab=paste("PC",i+1))
}
dev.off()
```
Los gráficos están en un pdf: `<pref>.PCAWithSubSet.pdf`
Add some interpreation

### Limitaciones del ACP con datos de ADN antiguo
El ACP aplicado exclusivamente a individuos antiguos puede producir resultados poco fiables. Esto ocurre porque el ADN antiguo presenta altos niveles de daño, fragmentación y, sobre todo, una gran proporción de datos faltantes. Estos problemas generan desplazamientos artificiales en los componentes principales y pueden distorsionar las relaciones biológicas reales entre individuos o poblaciones.
Existen, sin embargo, soluciones parcialmente efectivas. Una de ellas es usar métodos de proyección:
 - [**lsqproject**](https://github.com/DReichLab/EIG/blob/master/POPGEN/lsqproject.pdf)
 - y/o proyección de individuos antiguos sobre un ACP calculado únicamente con individuos modernos. De este modo, la estructura genética principal se define con datos completos y de alta calidad, y los genomas antiguos se ubican en ese espacio sin influir en la orientación de los ejes.
Sin embargo, la segunda solución suele ser limitada en poblaciones de Ámerica ya que la diversidad genética indígena de Ámerica de las poblaciones modernas suele ser poco representativa de la existente en poblaciones previas a la invasión europea.

Aun así, incluso con estas estrategias, el ACP puede seguir siendo sensible al patrón de datos faltantes típico del ADN antiguo. Por este motivo, una alternativa más robusta consiste en emplear un **MDS (Multidimensional Scaling)** basado en distancias genéticas derivadas del **f<sub>3</sub>-outgroup**, que proporciona una medida menos sesgada de divergencia genética. Este método se presenta como un ejercicio opcional al final del módulo.



## Análisis de estimación de proporciones de componentes genéticos

### Introducción al método

En esta sección utilizaremos el algoritmo **sNMF (Sparse Non-negative Matrix Factorization)**, descrito en  
[Frichot & François, 2015](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12382).  
Este método es conceptualmente similar a **ADMIXTURE** (https://dalexander.github.io/admixture/), pero suele ser más robusto y computacionalmente eficiente, especialmente cuando se trabaja con conjuntos de datos grandes o con genomas incompletos.

El objetivo de sNMF es resumir la variación genética observada en una matriz de genotipos mediante un conjunto de *K* componentes latentes. En términos prácticos, el método realiza los siguientes pasos:

1. **Identifica un número K de componentes genéticos**  
   Agrupa patrones de covariación genética entre individuos, produciendo *K* ejes o patrones básicos de estructura.

2. **Estima las frecuencias alélicas características de cada componente**  
   Esto se almacena en la matriz `.P`, que resume cómo contribuye cada variante genética a cada uno de los *K* componentes.

3. **Estima la contribución de cada componente en cada individuo**  
   Produce la matriz `.Q`, que indica, por ejemplo, si un individuo puede representarse como 40% componente 1, 30% componente 2, etc.

Es importante destacar que estos componentes **no representan poblaciones biológicas discretas ni grupos “puros”**, sino patrones matemáticos que permiten describir la estructura genética de forma compacta en un espacio de dimensión *K*.

> **En resumen:**  
> sNMF, igual que ADMIXTURE, descompone la variación genética en *K* patrones y estima cuánto contribuye cada uno a cada individuo. Es una representación flexible de la estructura genética que no presupone la existencia de poblaciones discretas.

### WEjecutar el algoritmo
Vamos a ejecutar el análisis variando el número de **componentes (K)** de 2 a 6.  
Para cada valor de **K** realizaremos **4 repeticiones independientes** con diferente inicialización.

> Nota: en estudios publicados se suelen hacer entre **10 y 30 repeticiones por K** para asegurar estabilidad y detectar modos alternativos; aquí usamos 4 repeticiones por motivos de tiempo, pero para un análisis definitivo conviene aumentar el número de replicados.

Para cada K compararemos las repeticiones usando el criterio de **cross-entropy** (error de predicción sobre genotipos enmascarados) y escogeremos la réplica con menor cross-entropy como “mejor corrida” para ese K. Más adelante también veremos cómo agrupar réplicas en “modos” cuando haya múltiples soluciones estables.

```r
Kmax=8
admProj = snmf(paste(pref,".geno",sep=""),
                K = 2:Kmax, 
                entropy = TRUE, 
                repetitions = 4,
                project = "new")
```

Los resultados del análisis se guardan automáticamente en una carpeta llamada  
`<pref>.snmf/`.  
Además, el objeto `admProj` que genera la función `snmf()` permite acceder de forma sencilla a todas las salidas: valores de cross-entropy, matrices **Q**, matrices **P**, número de iteraciones y demás información relevante.

### Selección del mejor *K* y de la mejor repetición por *K*

Para cada valor de *K*, el algoritmo genera varias repeticiones independientes.  
La métrica principal para evaluar y comparar estas repeticiones es el **cross-entropy**, que cuantifica qué tan bien el modelo puede predecir genotipos enmascarados.

En términos prácticos:

- El *K* óptimo suele encontrarse en el punto donde el cross-entropy **deja de disminuir de manera marcada** (el típico “punto de saturación”).  
- Dentro de cada *K*, la **mejor repetición** es simplemente aquella con el **cross-entropy más bajo**, es decir, la que mejor ajusta los datos.

Este procedimiento es análogo al de *ADMIXTURE*, donde se selecciona el *K* que minimiza el error de validación cruzada y, para ese *K*, la corrida con el menor error.  
El objeto `admProj` facilita extraer toda esta información y seleccionar de forma sistemática la solución óptima.


Para cada número de componentes (*K*), vamos a evaluar cómo las diferentes iteraciones del algoritmo explican los datos.  
La medida clave es la **entropía cruzada** (*cross-entropy*): cuanto más baja, mejor el ajuste del modelo.

> **¿Qué mide exactamente la entropía cruzada?**  
> El algoritmo *sNMF* busca aproximar la matriz de genotipos  
> **G ≈ Q × P**, donde:
> - **Q** resume, para cada individuo, cuánto contribuye cada uno de los *K* componentes genéticos;
> - **P** contiene las frecuencias alélicas características de esos componentes.
>
> Durante la optimización, *sNMF* compara los genotipos observados con los genotipos “reconstruidos” a partir de **Q** y **P**.  
> La discrepancia entre ambos se resume en un único valor: **la entropía cruzada**.
>
> - Una **entropía alta** indica que el modelo no logra explicar bien los datos con ese *K*  
>   (o que la iteración quedó atrapada en un óptimo local).  
> - Una **entropía baja** indica un mejor ajuste.
>
> Por eso, para cada *K*, la **mejor iteración** es la que tiene la entropía más baja.  
> Y comparar entropías entre diferentes valores de *K* es equivalente al procedimiento de **Cross-Validation en ADMIXTURE**.

```r
### for each K get bet score
listEntro<-list()
for(k in c(2:Kmax)){
  listEntro[[paste("K=",k,sep="")]]<-cross.entropy(admProj, K = k)
}
boxplot(listEntro,main="Cross-Entropopy per K across diferent iterations")
```
Vemos entonces que el modelo con, *K* = 2 es el que mejor se ajusta a los datos, independetemiente de las iteraciones. En general, es un patrón que se obsera cuando se analizan solo datos de ADN antiguo, por el tema de valores faltantes. Volveremos sobre este problema más adelante.

Podemos también seleccionar la iteración para cada <em>K</em> que explica mejor los datos según la entropía. 

```r
listBest<-list()
for(k in c(2:Kmax)){
  listBest[[paste("K=",k,sep="")]]<-which(listEntro[[paste("K=",k,sep="")]] == min(listEntro[[paste("K=",k,sep="")]]))
}
print(unlist(listBest))

```

### Visualización de los resultados

Vamos a graficar las estimaciones de los componentes genéticos.  

> Generar gráficos claros y comparables para este tipo de análisis suele ser trabajoso porque:
> 1. Es necesario **ordenar los individuos** de manera consistente con la estructura biológica esperada, para facilitar la interpretación visual.
> 2. Hay que **mantener colores coherentes entre distintos valores de *K***, ya que los componentes pueden cambiar de posición entre corridas.

Como primer paso, vamos a reorganizar los individuos según su región, utilizando la información contenida en el objeto `meta`.

```r
### preparamos el orden de los individuos en el gráfico
ind<-read.table(paste(pref,".ind",sep=""),stringsAsFactors=F,header=F)
names(ind)<-c("id","sex","pop")
metaPlot<-meta
row.names(metaPlot)<-metaPlot$id
metaPlot<-metaPlot[ind$id,]
metaPlot$orderInInd<-c(1:nrow(metaPlot))
metaPlot$orderPlot<-NA
metaTMP<-metaPlot[0,]
order=0
for(region in c("Beringia_EH","SouthernNorthAmerica_EH",
                 "Brazil_EH","Pampa_MH",
                 "CentralAndes_EH","CentralAndes_MH",
                 "CentralChile_EH","CentralChile_MH","CentralChile_LH",
                 "SouthPatagonia_MH","SouthPatagoniaM_LH","SouthPatagoniaB_LH","SouthPatagoniaT_LH")){
      
      tmp<-metaPlot[metaPlot$Region==region,]
      tmp<-tmp[ order(tmp$pop),]
      tmp$orderPlot<-c(1:nrow(tmp))+order
      metaTMP<-rbind(metaTMP,tmp)
      order=order+nrow(tmp)
      
}
metaPlot<-metaTMP
remove(metaTMP)

## rearrange labels (remove country and date)
change<-function(string){
    tmp<-strsplit(string,split="_")[[1]]
    tmp<-tmp[!grepl("BP",tmp)]
    #return(paste(tmp[c(2:length(tmp))],collapse="\n"))
    return(tmp[length(tmp)])
}
metaPlot$Label<-sapply(metaPlot$pop,change,USE.NAMES=F)

## keep track of individuals sequenced with shotgun 
metaPlot$SG<-grepl(".SG",metaPlot$pop)
```

Ahora vamos a leer los resultados para cada valor de *K* y reordenar la matriz **Q** de acuerdo con el orden deseado de individuos (almacenado en la columna `orderPlot` de la tabla que acabamos de generar).  
Luego graficaremos los resultados.

Para facilitar la interpretación, añadiremos también un punto negro bajo cada barra para identificar a los individuos secuenciados mediante *shotgun*.

```r
pdf(paste(pref,".SNMF.MixedColors.pdf",sep=""),height=10)
listColors<-c("darkgreen","cadetblue","darkorange4","red3","blue1","rosybrown","lightblue3","slateblue1")
par(mfrow=c(Kmax+1,1),mar=c(0.5,3,0.5,0.5))

posX=seq(1.5,nrow(metaPlot)-0.5,length.out=nrow(metaPlot))
plot(0,0,"n",axes=F,ann=F,xlim=c(1,nrow(ind)))
text(x=posX,y=rep(-1,nrow(metaPlot)),
        srt=90,adj=c(0,0.5),col=metaPlot$Color,
        labels=metaPlot$Label,cex=0.7)

for(k in c(Kmax:2)){
  best = listBest[[paste("K=",k,sep="")]]
  # display the Q-matrix
  Q.matrix <- t(as.matrix(Q(admProj, K = k, run = best)))
  Q.matrix<-Q.matrix[,metaPlot$orderInInd]
  bp<-barplot(Q.matrix, 
        border = NA, 
        ylab=paste("K =",k),
        axes=F,
        line=-1,
        space = 0, 
        col = listColors[c(1:k)], 
        )
        
}

plot(0,0,"n",axes=F,ann=F,xlim=c(1,nrow(ind)))
text(x=posX,y=rep(-1,nrow(metaPlot)),
        srt=90,adj=c(0,0.5),col=metaPlot$Color,
        labels=metaPlot$id,cex=0.7)
points(x=posX,y=rep(1,nrow(metaPlot)),
        pch=ifelse(metaPlot$SG,16,NA))
dev.off()
```

Estos gráficos muestran, para cada individuo, cómo se descompone su variación genética en *K* componentes inferidos por el modelo.  
Cada barra corresponde a un individuo y los colores representan la contribución relativa de cada componente.  
Así, patrones compartidos de colores entre individuos o regiones reflejan semejanzas genéticas que el algoritmo resumió en esos componentes.

> Los gráficos producidos por *sNMF* (o *Admixture*) **no mantienen colores coherentes entre distintos valores de *K***.  
> Para comparar fácilmente los resultados, necesitamos **alinear los componentes entre *Ks*** y aplicar una **misma paleta de colores fija**.  
> En este curso seguiremos la siguiente estrategia:
>
> - Elegimos un *K* de referencia (generalmente el *K* máximo) y asignamos un color fijo a cada uno de sus componentes, basándonos en la interpretación biológica preliminar.  
> - Para cada otro *K*, alineamos sus columnas con las del *K* de referencia utilizando el **algoritmo Húngaro**, que encuentra la combinación que maximiza la similitud entre componentes.

Empezamos asignando un color fijo a cada componente del *Kmax*. 
Este paso busca automatizar el proceso de mantener colores coherentes entre distintos valores de *K*.  
Sin embargo, es importante tener en cuenta que los resultados de *sNMF* pueden variar entre corridas, por lo que es posible que el código de abajo necesite pequeños ajustes manuales para lograr una correspondencia visualmente coherente en todos los *K* en las corridas realizadas en cada computadora.
```r
library(clue)

listColors<-c()
###read Q for Kmax
best = listBest[[paste("K=",Kmax,sep="")]]
# display the Q-matrix
Q.Kmax <- t(as.matrix(Q(admProj, K = Kmax, run = best)))
Q.Kmax<-data.frame(Q.Kmax[,metaPlot$orderInInd],stringsAsFactors=F)
names(Q.Kmax)<-metaPlot$id  

##assign "rosybrown" for the component maxized in Beringia
colBer<-which(Q.Kmax[,"USR1"]==max(Q.Kmax[,"USR1"]))
listColors[colBer]<-"rosybrown"

##assign "cadetblue" for the component maxized in SouthPatagonia Terrestre (e.g. I12364)
colTer<-which(Q.Kmax[,"I12364"]==max(Q.Kmax[,"I12364"]))
listColors[colTer]<-"cadetblue"

##assign "blue1" for the component maxized in SouthPatagonia Martime(e.g. I12942)
colMar<-which(Q.Kmax[,"I12942"]==max(Q.Kmax[,"I12942"]))
listColors[colMar]<-"blue1"

##assign "slateblue3" for the component maxized in other SouthPatagonia Terrestre(e.g. IPK13)
colMar<-which(Q.Kmax[,"IPK13"]==max(Q.Kmax[,"IPK13"]))
listColors[colMar]<-"slateblue3"

##assign "red3" for the component maxized in Central Chile (e.g. I1754)
colChil<-which(Q.Kmax[,"I1754"]==max(Q.Kmax[,"I1754"]))
listColors[colChil]<-"red3"

##assign "darkorange4" for the component maxized in CentralAndes (e.g. I0038)
colAnd<-which(Q.Kmax[,"I0038"]==max(Q.Kmax[,"I0038"]))
listColors[colAnd]<-"darkorange4"

##assign "darkgreen" for the component maxized in Brazil (e.g. Sumidouro7)
colBra<-which(Q.Kmax[,"Sumidouro7"]==max(Q.Kmax[,"Sumidouro7"]))
listColors[colBra]<-"darkgreen"

##assign "lightblue3" for the component maxized in MH southPatagonia (e.g. A460)
colMHPat<-which(Q.Kmax[,"A460"]==max(Q.Kmax[,"A460"]))
listColors[colMHPat]<-"lightblue3"

```

Ahora vamos a alinear los *clusters*. Para ello he creado la función `align_Q_down()` (archivo `Modulo3/Align_Q_down.R`) que, al pasar de *K* a *K−1*, realiza los siguientes pasos:

1. Calcula la **distancia euclidiana media** entre las proporciones estimadas del componente *i* en el modelo *K* y las del componente *j* en el modelo *K−1*.  
2. Con esas distancias construye una **matriz de costos** que cuantifica la disimilitud entre cada par de componentes (*i, j*).  
3. Aplica el **algoritmo húngaro** sobre la matriz de costos para obtener la asignación óptima de filas (componentes de *K−1*) a columnas (componentes de *K*).  
4. Devuelve las correspondencias encontradas, que usaremos para reordenar las matrices **Q** y mantener coherencia de colores entre los distintos valores de *K*.

```r

## in the file Modulo3/Align_Q_down.R, I defined a function trying to fit the best colors for K-1 to the ones used for K
source("Modulo3/Align_Q_down.R")
### now apply and plot
pdf(paste(pref,".SNMF.OrderedColors.pdf",sep=""),height=10)
Qref <- t(as.matrix(Q(admProj, K = Kmax, run = best)))
Qref<- Qref[,metaPlot$orderInInd]
par(mfrow=c(Kmax+1,1),mar=c(0.5,3,0.5,0.5))

posX=seq(1.5,nrow(metaPlot)-0.5,length.out=nrow(metaPlot))
plot(0,0,"n",axes=F,ann=F,xlim=c(1,nrow(ind)))
text(x=posX,y=rep(-1,nrow(metaPlot)),
        srt=90,adj=c(0,0.5),col=metaPlot$Color,
        labels=metaPlot$Label,cex=0.7)
bp<-barplot(Qref, 
        border = NA, 
        ylab=paste("K =",Kmax),
        line=-1,
        axes=F,
        space = 0, 
        col = listColors, 
        )

listColorsSub<-listColors
for(k in c((Kmax-1):2)){

    Qk <- t(as.matrix(Q(admProj, K = k, run = best)))
    Qk<-Qk[,metaPlot$orderInInd]
    returnAlign<-align_Q_down(Qref,Qk)
    ColAligned<-returnAlign[[2]]
    Qk<-returnAlign[[1]]
    listColorsSub<-listColorsSub[ sort(ColAligned) ]
    bp<-barplot(Qk, 
        border = NA, 
        ylab=paste("K =",k),
        line=-1,
        axes=F,
        space = 0, 
        col = listColorsSub, 
        )
      
    Qref<-Qk
    print(listColorsSub)

}

plot(0,0,"n",axes=F,ann=F,xlim=c(1,nrow(ind)))
text(x=posX,y=rep(-1,nrow(metaPlot)),
        srt=90,adj=c(0,0.5),col=metaPlot$Color,
        labels=metaPlot$id,cex=0.7)
points(x=posX,y=rep(1,nrow(metaPlot)),
        pch=ifelse(metaPlot$SG,16,NA))

dev.off()
```

Vemos entonces que así resulta más fácil interpretar las estructuras inferidas.  
Sin embargo, los colores no siempre corresponden perfectamente: para un cierto *K*, algunos individuos pueden mostrar una proporción cercana a 1 para un componente, pero esa señal puede fragmentarse cuando aumentamos *K*.  
Esto refleja que cada valor de *K* puede capturar **niveles diferentes de estructura poblacional**, y por lo tanto no existe una correspondencia totalmente rígida entre componentes de distintos *K*.








### Visualización de los diferentes modos por K

> Cuando usamos algoritmos como <em>sNMF</em>, para un mismo número de clusters (K), solemos realizar varias corridas independientes.
> Aunque todas las corridas usan el mismo <em>K</em>, el algoritmo puede llegar a soluciones ligeramente distintas porque:
> - empieza con diferentes condiciones iniciales,
> - la función que optimiza puede tener mínimos locales,
> - los datos pueden sostener más de una partición válida. <br>
> Estas diferentes soluciones estables se llaman **“modos”**. Un **modo** es un patrón consistente de agrupamiento que aparece de forma repetida en varias corridas, es decir una solución que el algoritmo encuentra varias veces porque los datos permiten interpretarse de más de una manera.

> Cada modo representa una interpretación alternativa de la estructura genética/variación estadística en la muestra.
> Visualizar los modos permite:
> - Detectar soluciones mayoritarias: el modo que aparece en la mayoría de corridas refleja la partición más estable o más apoyada por los datos.
> - Identificar soluciones minoritarias pero biológicamente interesantes: modos menos frecuentes pueden capturar patrones sutiles u homogéneos en subgrupos.
> - Distinguir variación real de ruido técnico: si hay muchos modos muy distintos y ninguno predominante, puede indicar que el <em>K</em> es demasiado grande, ya que los datos no soportan una partición clara, o que falta información.
> - Evitar seleccionar al azar una corrida “bonita”: cada corrida individual es solo una realización y los modos integran información sobre todas las corridas.

La herramienta quizás más conocida para realizar la identificación de los modos es [<em>PONG</em>](https://github.com/ramachandran-lab/pong), sin embargo es bastante facíl realizar lo mismo en R con los resuldatos de <em>sNMF</em>.

```r
####################
## PONG-like clustering of snmf replicates for all K. 
####################
## Vamos a usar la distancia euclidiana entre diferentes iteraciones pero se pueden usar otras métricas.

if(! require(proxy)){install.packages("proxy");require(proxy)}


# vamos a guardar los ficheros creados en <pref>.snmf_modes
dir.create(paste(pref,".snmf_modes",sep=""))

##################
# Definición de la función de detección de modo para un K
##################

process_K <- function(k,proj) {
    message("Processing K = ", k)

    nruns <- proj$Kproject[[as.character(k)]]$nrun
    message("Number of runs detected: ", nruns)

    # ----------------------------------------------------
    # Extract Q matrices from all runs
    # ----------------------------------------------------
    Qlist <- lapply(1:nruns, function(r) Q(proj, K = K, run = r))

    # ----------------------------------------------------
    # Align Q matrices (fix label switching)
    # ----------------------------------------------------
    aligned <- LEA::align_Q_matrices(Qlist)

    # Convert the aligned list into a single matrix
    # Each row = run, flattened Q-matrix
    flat <- lapply(aligned, function(Q) as.vector(t(Q)))
    flat <- do.call(rbind, flat)

    # ----------------------------------------------------
    # Compute distances among runs
    # ----------------------------------------------------
    D <- proxy::dist(flat, method = distance_method)

    # ----------------------------------------------------
    # Cluster runs into modes (hierarchical clustering)
    # ----------------------------------------------------
    hc <- hclust(D, method = "average")

    # Automatic mode detection:
    # rule: cut at 10% of maximum height (adjust if needed)
    threshold <- 0.10 * max(hc$height)
    modes <- cutree(hc, h = threshold)

    nmodes <- length(unique(modes))
    message("Identified modes for K=", K, ": ", nmodes)

    # ----------------------------------------------------
    # Compute average Q matrix for each mode
    # ----------------------------------------------------
    mode_Q <- lapply(unique(modes), function(m) {
        selected <- aligned[modes == m]
        Reduce("+", selected) / length(selected)
    })

    # ----------------------------------------------------
    # Save results
    # ----------------------------------------------------
    outdir <- paste0("snmf_modes/K", K)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    # Save clustering of runs
    write.table(
        data.frame(Run = 1:nruns, Mode = modes),
        file = file.path(outdir, "run_modes.txt"),
        row.names = FALSE, quote = FALSE
    )

    # Save mode-average Q matrices
    for (i in seq_along(mode_Q)) {
        write.table(
            mode_Q[[i]],
            file = file.path(outdir, paste0("mode_", i, "_Qmatrix.txt")),
            row.names = FALSE, col.names = FALSE, quote = FALSE
        )
    }

    # ----------------------------------------------------
    # Plot modes (simple barplot for each mode)
    # ----------------------------------------------------
    for (i in seq_along(mode_Q)) {

        df <- data.frame(
            Individual = 1:nrow(mode_Q[[i]]),
            mode_Q[[i]]
        )
        df_long <- reshape2::melt(df, id.vars = "Individual")

        p <- ggplot(df_long, aes(x = Individual, y = value, fill = variable)) +
            geom_bar(stat = "identity", width = 1) +
            theme_minimal() +
            ggtitle(paste("K =", K, " | Mode", i)) +
            ylab("Ancestry proportion")

        ggsave(file.path(outdir, paste0("mode_", i, "_plot.png")),
               p, width = 10, height = 4)
    }

    # Return summary
    list(K = K, nruns = nruns, nmodes = nmodes)
}


# --------------------------
# PROCESS ALL K VALUES
# --------------------------

summary_list <- lapply(Kmin:Kmax, process_K)
summary_table <- do.call(rbind, lapply(summary_list, as.data.frame))

write.table(summary_table,
            file = "snmf_modes/summary_modes_all_K.txt",
            row.names = FALSE, quote = FALSE)

message("\nDONE! Results stored in snmf_modes/\n")




```

### Discusión
#### Sobre el concepto de "Ancestría"
En la literatura se suele hablar de “ancestrías” para describir los componentes identificados con *ADMIXTURE* o *sNMF* (por ejemplo: “este individuo presenta un 40% de ancestría XXX”). Aunque esta forma de comunicar los resultados es práctica, ha sido criticada porque no refleja con exactitud lo que realiza el algoritmo ([Coop, 2022)[https://gcbias.org/wp-content/uploads/2022/07/genetic_similarity_and_genetic_ancestry_groups_current.pdf]. Además, el término “ancestría” puede sugerir la existencia de poblaciones genéticamente “puras”, evocando conceptos asociados a la idea de raza ([Kampourakis & Peterson, 2023)[https://doi.org/10.1093/genetics/iyad002]).
Si bien en publicaciones especializadas se continúa utilizando este término por conveniencia, es recomendable evitarlo en trabajos de divulgación científica, donde puede inducir interpretaciones erróneas o simplificaciones problemáticas.

### Sobre la interpretación de los resultados de estimación de proporciones de componentes genéticos
[Lawson et al. 2018](https://www.nature.com/articles/s41467-018-05257-7) advierieron sobre el riesgo de sobre interpretar los resultados de métodos como *ADMIXTURE* o *sNMF*.
Brevemente, advierten que diferentes hábitos para interpretar los resultados de este tipo de análisis pueden llevar a una mal interpretación de la historia evolutiva que se pretende explicar.
<img width="926" height="576" alt="image" src="https://github.com/user-attachments/assets/68c684df-96f6-476e-b196-5690404b05e6" />

A esto, se añade el problema que diferentes escenarios pueden llevar a observar resultados similares, como lo muestran en su figura 2 abajo mostrada.
<img width="703" height="699" alt="image" src="https://github.com/user-attachments/assets/7ce467e4-12b6-4785-a802-7383f74d3f0e" />




## BONUS: <em>f<sub>3</sub></em>-outgroup

Let's see 
```r
require(admixtools)
```



