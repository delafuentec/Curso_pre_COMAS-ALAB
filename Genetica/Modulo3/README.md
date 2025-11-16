# Analisis de diversidad y estructura genetica

En este modulo, vamos a realizar 2 analisis tipicos de genetica de poblaciones humanas:
1. Un Analisis en Componentes Principales (ACP)
2. Una estimación de coeficientes de estructura genética

> Estos analisis se suelen hacer con <em>smartpca</em> del software [<em>EIGENSOFT</em>](https://github.com/DReichLab/EIG) y <em>ADMIXTURE</em>[https://dalexander.github.io/admixture/], pero su implementación en el paquete <em>LEA</em> de **R** son muy similares.
Además el algoritmo para las estimaciones de coeficientes de estructura genética es más robusto y rápido que <em>ADMIXTURE</em>.

Vamos a leer los datos filtrados que generamos en el modulo anterior.

```r
require(LEA)
### you maybe need to change the directory below
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB/Genetica/")
pref="Modulo3/filteredDataSet"

```

## ACP
En este analísis, se reduce el número de dimensiones en grandes conjuntos de datos a componentes principales que conservan la mayor parte de la información original.
### Primer ACP a escala continental
Vamos a realizar un primer analisis incorporando todos los individuos antiguos que conforman nuestra conjunto de datos filtrado.
Podemos explorar los individuos analizados leyendo el fichero **ind** y cruzandolo con los metadatos
```r
ind<-read.table(paste(pref,".ind",sep=""),stringsAsFactors=F,header=F)
names(ind)<-c("id","sex","pop")
meta=read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.metadataPerind.txt",stringsAsFactors = F,header=T)
### generate table of the number of individual per region
table(meta$Region[ meta$id %in% ind$id])
#### Correr el ACP
```r
pc = pca(paste(pref,".geno",sep=""), scale = TRUE)
```
Los resultados se crean en una carpeta <pref>.pca. 
```r
### we can read the eigenvalues file
eigenvalues<-read.table(paste(pref,".pca/filteredDataSet.eigenvalues",sep=""),stringsAsFactors=F,header=F)
### let see the porcentage of variance explained per principal components
plot(eigenvalues$V1/sum(eigenvalues$V1)*100, lwd=5, col="mediumseagreen",xlab=("PCs"),ylab="% variance explained")
```

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
#pdf(paste(pref,".PCAWithAll.pdf",sep=""),height=10)
par(mfrow=c(3,2))
forLeg<-unique(projections[,c("Region","pop","Color","Point")])
forLeg<-forLeg[ order(forLeg$Region,forLeg$pop),]
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=forLeg$Point,pt.bg=forLeg$Color,col=ifelse(forLeg$Point<21,forLeg$Color,"black"),
        legend=paste(forLeg$Region,forLeg$pop),ncol=1,cex=0.4,pt.lwd=0.5)
for(i in seq(1,9,2)){
  plot(projections[,paste("V",i,sep="")],projections[,paste("V",i+1,sep="")],
          pch=projections$Point,
          bg=projections$Color,
          col=ifelse(projections$Point<21,projections$Color,"black"),
          xlab=paste("PC",i),
          ylab=paste("PC",i+1))
}
#dev.off()
```
Se generaron los graficos de PC2 vs PC1, PC4 vs PC3, ..., P10 vs PC9.


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

Vemos entonces que es Sumidouro5 (~10400BP) que se diferencia de los otros. Sumidouro 5, un individuo del Holoceno temprano proveniente de la cueva de Sumidouro (Lagoa Santa, Brasil), destaca por varias razones biológicas y arqueológicas. En primer lugar, su morfología craneal fue reconocida desde mucho tiempo como particularmente distintiva: comparaciones multivariantes indicaron afinidades con poblaciones africanas y austro-melanesias, en lugar de con los indigenas  americanos tardíos o asiáticos [(Neves et al. 2007)](bpb-us-w2.wpmucdn.com). <br>
Más recientemente, datos genómicos han revelado en Sumidouro 5 una señal de ascendencia mitocondrial D4h3a y un perfil autosómico que no se corresponde completamente con otros individuos tempranos de Lagoa Santa, lo que sugiere conexiones genéticas más amplias (incluso con poblaciones amazónicas ancestrales), ([Moreno-Mayar et al. 2018](https://www.science.org/doi/10.1126/science.aav2621),[Ferraz et al. 2023](https://doi.org/10.1038/s41559-023-02114-9)). <br>
Esta combinación de rasgos morfológicos arcaicos y una ascendencia genética poco común indica que Sumidouro 5 pudo pertenecer a un linaje parcialmente distinto a otros pobladores tempranos de Lagoa Santa, aportando evidencia de una mayor diversidad biológica en los primeros habitantes del este de Sudamérica y de posibles rutas de migración complejas durante el Holoceno temprano.

### ACP enfocandose en procesos mas tardios de Sudmerica
Ahora vamos a rehacer el analisis pero sacando los individuos de Norte America, y los grupos tempranos de Chile (Los Rieles ~12000BP) y de Lagoa Santa (Brazil).
Para eso vamos a leer los datos al formato *eigenstrat* para generar un sub-conjunto sin estos individuos

#### Generar sub-conjunto
Primero vamos a convertir los datos en un objeto `genLight`.
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

Ahora, vamos a generar la lista de individuos que queremos dejar en el subconjunto de datos. 

```r
# get list of inds to remove 
toRemoveID<-meta$id [ meta$Region %in% c("Beringia_EH","CentralChile_EH","SouthernNorthAmerica_EH") | 
		      grepl("Brazil_Sumidouro",meta$pop)]
# make list of remaining individuals
toKeepID<-ind$id[ ! ind$id %in% toRemoveID]
# check from which regions come from the remaining individuals
print(table(meta$Region[ meta$id %in% toKeepID]))
# check the remaining groups from Brazil_EH
print(table(meta$pop[ meta$Region == "Brazil_EH" & meta$id %in% toKeepID]))

##looks ok --> generate subdataset, removing monorphisms (need to recalculate the metrics for further filtering)
```

Ahora podemos generar el subconjunto de datos. Por las dudas hay que volver a filtrar las variantes monomorficas y con alta tasa de valores faltantes.
Y sacaremos los individuos con menos de 10000 SNPs
```r
gl<-gl.keep.ind(gl,ind.list =toKeepID,recalc = T,mono.rm = T)
gl<-gl.filter.callrate(gl, threshold=2/nInd(gl),recalc=T)

### 


###check numbers
print(paste("from",nrow(ind),"individuals, we wanted to keep",length(toKeepID),", and the sub-dataset has ",nInd(gl)))
print(paste("from",nrow(snp),"snps",nLoc(gl)," remain in the sub-dataset"))


```
Ahora podemos escribir los datos del sub-conjunto , tal como lo hicimos en el Modulo2. Los vamos a escribir en ficheros llamados `<pref>_subset.{geno,snp,ind}`.

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

Ahora podemos hacer el ACP en este sub-conjunto
```r
pc_sub = pca(paste(pref,"_subset.geno",sep=""), scale = TRUE)
```
Los resultados se crean en una carpeta <pref>.pca. 
```r
### we can read the eigenvalues file
eigenvalues_sub<-read.table(paste(pref,"_subset.pca/filteredDataSet_subset.eigenvalues",sep=""),stringsAsFactors=F,header=F)
### let see the porcentage of variance explained per principal components
plot(eigenvalues_sub$V1/sum(eigenvalues_sub$V1)*100, lwd=5, col="red",xlab=("PCs"),ylab="% variance explained")
```

```r

### we can read the projections of the individuals: one line per individual.
### the individuals are ordered as in the input file (thus reading <pref>.ind.txt we know which column corresponds to which individual)
projections_sub<-read.table(paste(pref,"_subset.pca/",strsplit(pref,split="/")[[1]][2],"_subset.projections",sep=""),stringsAsFactors=F,header=F)
## add info per ind
projections_sub$id=ind_sub$id
projections_sub$pop=ind_sub$pop
###get metafile for plotting
projections_sub<-merge(projections_sub,meta,by=c("id","pop"))

###let's plot the first 10 PCs
#pdf(paste(pref,".PCAWithSubSet.pdf",sep=""),height=10)
par(mfrow=c(3,2))
forLeg_sub<-unique(projections_sub[,c("Region","pop","Color","Point")])
forLeg_sub<-forLeg_sub[ order(forLeg_sub$Region,forLeg_sub$pop),]
plot(0,0,"n",axes=F,ann=F)
legend("center",pch=forLeg_sub$Point,pt.bg=forLeg_sub$Color,col=ifelse(forLeg_sub$Point<21,forLeg_sub$Color,"black"),
        legend=paste(forLeg_sub$Region,forLeg_sub$pop),ncol=1,cex=0.4,pt.lwd=0.5)
for(i in seq(1,9,2)){
  plot(projections_sub[,paste("V",i,sep="")],projections_sub[,paste("V",i+1,sep="")],
          pch=projections_sub$Point,
          bg=projections_sub$Color,
          col=ifelse(projections_sub$Point<21,projections_sub$Color,"black"),
          xlab=paste("PC",i),
          ylab=paste("PC",i+1))
}
#dev.off()
```

Add some interpreation

> speak of mssing data and PCA


## Estimaciones de coeficientes de estructura genética
### Introducción al método
Vamos a usar el algoritmo  sNMF (Sparse Non-negative Matrix Factorization) descrito en [Frichot & François, 2015](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12382).
Similar a <em>ADMIXTURE</em>[https://dalexander.github.io/admixture/], <em>sNMF</em> es un método que toma una  matriz de genotipos y trata de resumir la variación genética en un conjunto pequeño de componentes.
En términos prácticos, sNMF:
1. Identifica K componentes genéticos:
    patrones de variación que aparecen repetidamente en el conjunto de individuos.)
2. Calcula para cada individuo cuánto contribuye cada componente:
    Genera una matriz que indica, por ejemplo, si un individuo puede describirse como 40% del componente 1, 30% del componente 2, etc.
3. Estima las frecuencias alélicas típicas de cada componente:
    Es decir, cómo se ve genéticamente cada uno de esos patrones.

No asume poblaciones discretas ni grupos biológicos "reales", simplemente resume los datos genéticos en un número K de patrones que ayudan a describir la estructura multidimensional.

> En resumen:
Tanto <em>sNMF</em> como <em>ADMIXTURE</em>  descomponen los datos genéticos en K patrones de variación y estima cuánto contribuye cada patrón a cada individuo.
Es una forma de representar estructuras complejas de variación sin asumir poblaciones discretas.

### Correr el algoritmo
Vamos a realizar el analisis por un numero de **poblaciones ancestral** (o cluster) **K** variando de 2 a 6. Para cada **K** haremos 4 repeticiones independientes.
> Se suelen hacer en los estudios entre 10 y 30 repeticiones por **K**.


```r
Kmax=6
admProj = snmf(paste(pref,".geno",sep=""),
                K = 2:Kmax, 
                entropy = TRUE, 
                repetitions = 4,
                project = "new")
```
Los resultados se escriben en una carpeta `<pref>.snmf/` y el objeto creado (`admProj` permite acceder facilmente a los mismos).

### Selección del mejor <em>K</em> y de la mejor iteración por <em>K</em>

Para cada número de poblacionaes ancestrales (clusters), <em>K</em>, vamos a ver cómo las diferentes iteraciones explican los datos (tiene menor entropía crzuada).
> Cada corrida del algoritmo produce un valor llamado cross-entropy o entropy score. <br>
> Este valor mide qué tan bien el modelo con un número dado de clusters (K) explica los datos genéticos observados.<br>
> Cuanto más baja sea, mejor ajusta el modelo a los datos.<br>
> **¿Qué significa exactamente?**  El algoritmo snmf intenta factorizar la matriz de genotipos **G ≈ Q × F***, con :
> - una matriz ***Q***, que representa proporciones de ancestría por individuo; y
> - una matriz ***F***, que representa frecuencias alélicas por población ancestral. <br>
> Durante el cálculo, <em>sNMF</em> evalúa qué tan bien estas matrices reconstruyen los genotipos observados.
> Esa diferencia entre lo observado y lo esperado se resume en un único valor: la entropía cruzada.
> Si la entropía es alta, significa que el modelo no puede explicar bien la variación genética con ese número de clusters (K), o esa corrida del algoritmo se quedó atrapada en un óptimo local.
> Si la entropía es baja, significa que el modelo explica mejor la estructura genética.
> Por eso, la mejor corrida para un K dado es la que tiene el menor valor de entropía.
> **La entropía cruzada es equivalente al escore de Cross-Validation en Admixture**

```r
### for each K get bet score
listEntro<-list()
for(k in c(2:Kmax)){
  listEntro[[paste("K=",k,sep="")]]<-cross.entropy(admProj, K = k)
}
boxplot(listEntro,main="Cross-Entropopy per K across diferent iteracions")
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
Vamos a graficar los resultados. 
> Es siempre engoroso generar manualmente gráficos lindos de estimaciones de coeficientes de estructura genética: 
> 1. hay que ordenar los individuos según el sentido biologíco esperado para poder visualizar mejor lo que puede explicar la estructura genética capturada
> 2. hay que hacer corresponder los colores de cada componente a diferentes <em>K</em>.

Vamos a reorganizar primero los individuos por región tal como guardado en el objeto `meta`.

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

Ahora vamos a leer los resultados para cada <em>K</em> y reordenar la matriz Q para que los individuos estén en el orden deseado (tal como guardado en la columna `orderPlot` de la tabla que acabamos de generar), y graficar.
Añadimos un punto negro (rótulo abajo) para los individuos secuenciados con shotgun.
```r
listColors<-c("darkgreen","cadetblue","darkorange4","red3","blue1","rosybrown")
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
        line=-1,
        axes=F,
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

```

> Los gráficos de *sNMF* (o *Admixture*) NO usan colores coherentes entre diferentes valores de *K*. Para obtener colores consistentes, debemos alinear los clusters entre *Ks* y usar una misma paleta fija. <br>
> Brevemente, vamos a seguir la siguiente estrategia
> - Elegimos un *K* de referencia (mejor que sea Kmax), y asignamos un color a cada componente en base a lo que observamos previamente.
> - Para cada otro *K), alineamos sus columnas con las del K de referencia usando el algoritmo Húngaro (maximize correlation).

Empezamos con la asignación de un color a cada componente para el *Kmax* . Es un intento de automización, pero es posible que el código aabjo requiera ediciones porque los resultados que se obtienen en cada corrida pueden diferir.
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
##assign "red3" for the component maxized in Central Chile (e.g. I1754)
colChil<-which(Q.Kmax[,"I1754"]==max(Q.Kmax[,"I1754"]))
listColors[colChil]<-"red3"
##assign "darkorange4" for the component maxized in CentralAndes (e.g. I0038)
colAnd<-which(Q.Kmax[,"I0038"]==max(Q.Kmax[,"I0038"]))
listColors[colAnd]<-"darkorange4"
##assign "darkgreen" for the component maxized in Brazil (e.g. Sumidouro7)
colBra<-which(Q.Kmax[,"Sumidouro7"]==max(Q.Kmax[,"Sumidouro7"]))
listColors[colBra]<-"darkgreen"
```

Ahora vamos a tratar de alinear los clusteres. Para esto, escribí una función `align_Q_down()` guardada en el ficher `Modulo3/Align_Q_down.R` que hace lo siguiente:
Yendo de *K* a *K−1*:
1. Calcula el promedio de las distancias euclidianas entre las proporciones estimadas para el componente i del modelo con K in el componente j del modelo con K-1.
2. Esto genera  una matriz de confusión para buscar los clusteres del modelo con K-1 grupos correspondientes a clusteres con modelo con K grupos.
3. Aplica un algoritmo hungaro para asiñar las mejores corerspondencias en base a esta matriz de confusión
4. Retorna las correspondencias

```r

## in the file Modulo3/Align_Q_down.R, I defined a function trying to fit the best colors for K-1 to the ones used for K
source("Modulo3/Align_Q_down.R")
### now apply
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


```

Vemos entonces que es más facíl interpretar!
Sin embargo, a veces los colores no corresponden del todo, porque a un un K dado, unos individuos tienen proporciones cerca del 1 para un componente, y ya no para K superiores. Esto demuestra que cada K puede capturar diferentes estructuras. Ver por ejemplo el caso de los Arquipelagos occidentales con K = 5 y 4 vs con otros K.

AÑADIR INTERPRETACION






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
#### Sobre concepto de "Amcestría"
> En la literatura se suele hablar de “ancestrías” para describir estos componentes (por ejemplo: “este individuo presenta un 40% de ancestría XXX”). Aunque esta forma de comunicar los resultados es práctica, ha sido criticada porque no refleja con exactitud lo que realiza el algoritmo ([Coop, 2022)[https://gcbias.org/wp-content/uploads/2022/07/genetic_similarity_and_genetic_ancestry_groups_current.pdf]. Además, el término “ancestría” puede sugerir la existencia de poblaciones genéticamente “puras”, evocando conceptos asociados a la idea de raza ([Kampourakis & Peterson, 2023)[https://doi.org/10.1093/genetics/iyad002]).
> Si bien en publicaciones especializadas se continúa utilizando este término por conveniencia, es recomendable evitarlo en trabajos de divulgación científica, donde puede inducir interpretaciones erróneas o simplificaciones problemáticas.

### Sobre la interpretación de los resultados de estimaciones de XXX
> [Lawson et al. 2018](https://www.nature.com/articles/s41467-018-05257-7) advierieron sobre el riesgo de sobre interpretar los resultados de métodos como *ADMIXTURE* o *sNMF*.
> Brevemente, diferentes escenarios pueden llevar a observar resultados similares, como lo muestran en su figura 2.



## BONUS: <em>f<sub>3</sub></em>-outgroup

Let's see 
```r
require(admixtools)
```



