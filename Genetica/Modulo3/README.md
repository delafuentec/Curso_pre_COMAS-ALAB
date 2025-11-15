# Analisis de diversidad y estructura genetica

En este modulo, vamos a realizar 2 analisis tipicos de genetica de poblaciones humanas:
1. un Analisis en Componentes Principales (ACP)
2. Un analisis de Admixture 

> Estos analisis se suelen hacer con <em>smartpca</em> del software [<em>EIGENSOFT</em>](https://github.com/DReichLab/EIG) y <em>ADMIXTURE</em>[https://dalexander.github.io/admixture/], pero su implementación en el paquete <em>LEA</em> de **R** son muy similares.
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
plot(eigenvalues$V1/sum(eigenvalues$V1)*100, lwd=5, col="red",xlab=("PCs"),ylab="% variance explained")
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
pdf(paste(pref,".PCAWithSubSet.pdf",sep=""),height=10)
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
dev.off()
```


## f<sub>3</sub>~-outgroup

## Admixture


Vamos a realizar el analisis por un numero de **poblaciones ancestral** (o cluster) **K** variando de 2 a 10. Para cada **K** haremos 4 repeticiones independientes.
> Se suelen hacer en los estudios entre 10 y 30 repeticiones por **K**.



K```r
admProj = snmf(paste(pref,".geno",sep=""),
                K = 2:10, 
                entropy = TRUE, 
                repetitions = 4,
                project = "new")

```
