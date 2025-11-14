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
En este analísis, reduce el número de dimensiones en grandes conjuntos de datos a componentes principales que conservan la mayor parte de la información original.
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
ind<-read.table(paste(pref,".ind",sep=""),stringsAsFactors=F,header=F)
projections$id=ind$V1
projections$pop=ind$V3
###get metafile for plotting
meta=read.table("PatagoniaDataSetWithOutgroups_ALAB2025/Ancient.metadataPerind.txt",stringsAsFactors = F,header=T)
projections<-merge(projections,meta,by=c("id","pop"))

###let;s plot
pdf(paste(pref,".PCAWithAll.pdf",sep=""),height=10)
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
dev.off()
### the
### let see the porcentage of variance explained per principal components
plot(eigenvalues$V1/sum(eigenvalues$V1)*100, lwd=5, col="red",xlab=("PCs"),ylab="% variance explained")
```






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