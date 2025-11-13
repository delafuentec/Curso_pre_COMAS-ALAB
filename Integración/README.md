AGREGAR INTRODUCCION




## Analísis genéticos

Utilizaremos el paquete de R admixtools para estimar el outgroup-ƒ3 entre todos los posibles pares de grupos en la base de datos. Al igual que en el módulo 4 de Genética:

```r
# Cargar modulos 
library(admixtools) 
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(heatmap3)
library(viridis)
library(gridExtra)
library(ggrepel)
library(stringr)
library(patchwork)

# Especificar directorio de trabajo (cambiar según corresponda)
setwd("/pasteur/helix/projects/Hotpaleo/pierre/Projects/Cours/ALAB_2025/Curso_pre_COMAS-ALAB/Integración/")

# Definir input
prefix <- "AmericaByGroups/Genodata_forCompMorphology"

# Definir poblaciones a analizar. Se utilizaran todos los pares posibles de grupos
popA = c('PATAGONIA_MARITIMO',
         'PUERTO_RICO',
         'BRASIL_TARDIO',
         'BRASIL_TEMPRANO',
         'CUBA',
         'CENTRO_PERU',
         'NORTE_PERU',
         'CALIFORNIA',
         'PATAGONIA_TERRESTRE',
         'PERU_TempranoHoloceno',
         'SUR_PERU',
         'PERU_MedioHoloceno',
         'PERICUES',
         'PAMPA',
         'ALASKA')
popB = c('PATAGONIA_MARITIMO',
         'PUERTO_RICO',
         'BRASIL_TARDIO',
         'BRASIL_TEMPRANO',
         'CUBA',
         'CENTRO_PERU',
         'NORTE_PERU',
         'CALIFORNIA',
         'PATAGONIA_TERRESTRE',
         'PERU_TempranoHoloceno',
         'SUR_PERU',
         'PERU_MedioHoloceno',
         'PERICUES',
         'PAMPA',
         'ALASKA')
outg = c("Mbuti")

# Correr f3
f3res <- f3(prefix, outg, popA, popB,
            outgroupmode = TRUE)
# Filtrar resultados en donde pop2 es la misma que pop3
f3res = f3res %>% filter(pop2 != pop3)

```


### Estimación de distancias genéticas y cálculo de MDS

```r
# Construir matriz de distancias genéticas
distGenMDS <- acast(f3res, pop2 ~ pop3, value.var = "est")
diag(distGenMDS) <- 1
distGenMDS <- 1 - distGenMDS #### / max(mdf, na.rm = TRUE)
### we remove ALASKA (we only will use it as outgroup for trees)
distGenMDS<-distGenMDS[row.names(distGenMDS) != "ALASKA" , colnames(distGenMDS) != "ALASKA"]

# Hacer análisis de escalamiento multidimensional (MDS)
# let's do it on 10 dimensions
k=5
GenMDS <- as.data.frame(cmdscale(distGenMDS, k = k))
GenMDS$groupId <- rownames(distGenMDS)
names(GenMDS)[c(1:k)]<-paste0("Dim",c(1:k))
# Agregaremos información de región para facilitar visualización
metadata = data.frame(groupId = c('PATAGONIA_MARITIMO',
                                  'PUERTO_RICO',
                                  'BRASIL_TARDIO',
                                  'BRASIL_TEMPRANO',
                                  'CUBA',
                                  'CENTRO_PERU',
                                  'NORTE_PERU',
                                  'CALIFORNIA',
                                  'PATAGONIA_TERRESTRE',
                                  'PERU_TempranoHoloceno',
                                  'SUR_PERU',
                                  'PERU_MedioHoloceno',
                                  'PERICUES',
                                  'PAMPA'),
                      region = c('Patagonia',
                                 'Caribe',
                                 'Brasil',
                                 'Brasil',
                                 'Caribe',
                                 'Andes Centro-Sur',
                                 'Andes Centro-Sur',
                                 'Norteamerica',
                                 'Patagonia',
                                 'Andes Centro-Sur',
                                 'Andes Centro-Sur',
                                 'Andes Centro-Sur',
                                 'Mesoamerica',
                                 'Pampa'))

colorVec<-c("orange3","mediumseagreen","darkgrey","deeppink", "rosybrown2","steelblue2",'darkblue')
GenMDS$region = metadata$region[match(GenMDS$groupId, metadata$groupId)]

# Grafico de MDS
# Graficaremos todos los posibles pares de dimensiones
dims <- paste0("Dim",c(1:k))
pairs <- combn(dims, 2, simplify = FALSE)

plots <- lapply(pairs, function(p) {
  x <- p[1]
  y <- p[2]
  ggplot(GenMDS, aes_string(x = x, y = y, color = "region",shape="groupId")) +
    geom_point(size = 3) +
    theme_bw() +
    scale_shape_manual(values = c(1:15)) +
    scale_color_manual(values = colorVec) +
    labs(x = x, y = y, shape = "groupId") +
    ggtitle(paste("Genomic\n",x, " vs ", y,sep=""))+
    theme(legend.position = "none", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
})

# Combinar plots
wrap_plots(plots, ncol = 3)
#ggsave("mds_america.png", width = 15, height = 15)
```

Make legend
```r
metadata$groupId<-as.factor(metadata$groupId)
metadata$Point[ order(metadata$groupId)]<-c(1:14)
metadata$region<-as.factor(metadata$region)
regN=0
for(region in unique(sort(metadata$region))){
  regN=regN+1
  metadata$Color[metadata$region ==region]=colorVec[regN]
}

plot(0,0,"n",axes=F,ann=F)
metadata<-metadata[order(metadata$region,metadata$groupId),]
legend("center",pch=metadata$Point,col=metadata$Color,legend=metadata$groupId,ncol=2)

metadata<-rbind(metadata,
                cbind("groupId"="ALASKA","region"="NORTE_AMERICA","Point"=16,"Color"="rosybrown2"))
row.names(metadata)<-metadata$groupId
```


### Arbol de Neighbor-joining con datos genéticos

```r

require(ape)
# Construir matriz de distancias genéticas
distGenNJ <- acast(f3res, pop2 ~ pop3, value.var = "est")
distGenNJ <- 1/distGenNJ
diag(distGenNJ) <- 0
distGenNJ<-distGenNJ/max(distGenNJ)
distGenNJ<-as.dist(distGenNJ)
GenNJ<-bionj(distGenNJ)
GenNJ<-root(GenNJ,outgroup = "ALASKA",resolve.root = F)
if(sum(! GenNJ$tip.label %in% row.names(metadata) )>0){
  print(paste(MorfNJ$tip.label[! GenNJ$tip.label %in% row.names(metadata)]))
  stop("pb groups not in your metadata")
}
write.tree(GenNJ,file="NJwithGenomicData.nwk",append = F,tree.names = F)
GenNJplot<-read.tree("NJwithGenomicData.nwk")
dicoPLOT<-metadata[GenNJplot$tip.label,]
plotStats<-plot(GenNJplot,tip.color = dicoPLOT$Color,cex = 1,use.edge.length = F,show.tip.label=T,label.offset=5,align.tip.label=T,main="NJ with genomic data")
points(x=rep(nrow(dicoPLOT)+2.5,nrow(dicoPLOT)),y=c(1:nrow(dicoPLOT)),
       cex=1,pch=as.numeric(dicoPLOT$Point),col=dicoPLOT$Color)



### we remove ALASKA (we only will use it as outgroup for trees)
distGenNJ_noOut <- acast(f3res, pop2 ~ pop3, value.var = "est")
distGenNJ_noOut <- 1/distGenNJ_noOut
diag(distGenNJ_noOut) <- 0
distGenNJ_noOut<-distGenNJ_noOut[row.names(distGenNJ_noOut) != "ALASKA" , colnames(distGenNJ_noOut) != "ALASKA"]
distGenNJ_noOut<-as.dist(distGenNJ_noOut)

GenNJ_noOut<-bionj(distGenNJ_noOut)
plot(GenNJ_noOut,use.edge.length = F,show.tip.label=T,align.tip.label=T)#,tip.color = dicoPLOT$Color,label.offset=10,cex = 0.3,)

```

## Analísis morfólogicos

### Formateo distancias morfólogicas y calculo MDS
```r

distMorf<-read.csv("AmericaByGroups/MorphologicalDistances.tsv",stringsAsFactors=F,header=T)
row.names(distMorf)<-distMorf$X
distMorf<-distMorf[,-1]
distMorf<-distMorf/max(distMorf)
distMorfMDS<-distMorf[row.names(distMorf) != "ALASKA" , colnames(distMorf) != "ALASKA"]

# Hacer análisis de escalamiento multidimensional (MDS)
MorfMDS <- as.data.frame(cmdscale(distMorfMDS, k = k))
MorfMDS$groupId <- rownames(distMorfMDS)

MorfMDS$region = metadata$region[match(MorfMDS$groupId, metadata$groupId)]

row.names(MorfMDS)<-MorfMDS$groupId
MorfMDS<-MorfMDS[ GenMDS$groupId,]
names(MorfMDS)[c(1:k)]<-paste0("Dim",c(1:k))
# Grafico de MDS
# Graficaremos todos los posibles pares de dimensiones
dims <- paste0("Dim",c(1:k)
pairs <- combn(dims, 2, simplify = FALSE)

plots <- lapply(pairs, function(p) {
  x <- p[1]
  y <- p[2]
  ggplot(MorfMDS, aes_string(x = x, y = y, color = "region",shape="groupId")) +
    geom_point(size = 3) +
    theme_bw() +
    scale_shape_manual(values = c(1:15)) +
    scale_color_manual(values = colorVec) +
    labs(x = x, y = y, shape = "groupId") +
    ggtitle(paste("Morphologic\n",x, " vs ", y,sep=""))+
    theme(legend.position = "none", 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())
})

# Combinar plots
wrap_plots(plots, ncol = 3)
```


### Arbol de Neighbor-joining desde distancias morfólogicas

```r
# Construir matriz de distancias genéticas
distMorfNJ <- as.dist(distMorf)

MorfNJ<-bionj(distMorfNJ)
MorfNJ<-root(MorfNJ,outgroup = "ALASKA",resolve.root = F)
if(sum(! MorfNJ$tip.label %in% row.names(metadata) )>0){
  print(paste(MorfNJ$tip.label[! MorfNJ$tip.label %in% row.names(metadata)]))
  stop("pb groups not in your metadata")
}
write.tree(MorfNJ,file="NJwithMorphologicalData.nwk",append = F,tree.names = F)
MorfNJplot<-read.tree("NJwithMorphologicalData.nwk")
dicoPLOT<-metadata[MorfNJplot$tip.label,]
plotStats<-plot(MorfNJplot,tip.color = dicoPLOT$Color,cex = 1,use.edge.length = F,show.tip.label=T,label.offset=5,align.tip.label=T,main="NJ with morphological data")
points(x=rep(nrow(dicoPLOT)+2.5,nrow(dicoPLOT)),y=c(1:nrow(dicoPLOT)),
       cex=1,pch=as.numeric(dicoPLOT$Point),col=dicoPLOT$Color)






### we remove ALASKA (we only will use it as outgroup for trees)
distMorfNJ_noOut <- distMorf[row.names(distMorf) != "ALASKA" , colnames(distMorf) != "ALASKA"]
distMorfNJ_noOut <- as.dist(distMorfNJ_noOut)

MorfNJ_noOut<-bionj(distMorfNJ_noOut)

plot(MorfNJ_noOut,use.edge.length = F,show.tip.label=T,align.tip.label=T)#,tip.color = dicoPLOT$Color,label.offset=10,cex = 0.3,)


par(mfrow=c(3,1))
for(m in c("ward.D","ward.D2","complete")){
  plot(hclust(distMorfNJ,method=m),main=paste("morfologia con ",m))
}
```



## Comparaciones al fin

Esto es un ensayo. Son análisis explorativos que tratan de responder a la necesidad de dialogo entre las disciplinas. Más allá de articular conclusiones desde diferentes evidencias, pretendemos iniciar una discusión epistemiologíca para articular los datos de diferente índole en un análisis interdisciplinar.
### MDS analysis with procrustes

```r
require(vegan)


### find best match among dimensions (Dim1 from GenMDS corresponds better to DimXX from MorfMDS, and so on)

PROC <- procrustes(GenMDS[,c(1:k)],MorfMDS[,c(1:k)])
Genomic <- as.data.frame(PROC$X)
Genomic$dataset <- "Genomic"
Genomic$groupId <- rownames(PROC$X)

Morphologic <- as.data.frame(PROC$Yrot)
Morphologic$dataset <- "Morphologic"
Morphologic$groupId <- rownames(PROC$Yrot)
names(Morphologic)[c(1:k)]<-paste0("Dim",c(1:k))

# Combine
plot_df <- rbind(Genomic, Morphologic)
plot_df$region = metadata$region[match(plot_df$groupId, metadata$groupId)]
names(plot_df)[c(1:k)]<-paste0("Dim",c(1:k))
# Optional: add a region/color variable (must be aligned with ids)
# plot_df$region <- your_region_vector[match(plot_df$id, names(your_region_vector))]

plot_df$groupId=as.factor(plot_df$groupId)

dims <- paste0("Dim",c(1:5))
pairs <- combn(dims, 2, simplify = FALSE)


# Plot with displacement lines
plots <- lapply(pairs, function(p) {
  x <- p[1]
  y <- p[2]
  ggplot(plot_df, aes_string(x = x, y = y, color = "region",shape="groupId")) +
  geom_path(aes(group = groupId), linewidth = 0.4, alpha = 0.7,arrow=arrow(type = "closed",length=unit(0.1,"inches"))) +   # displacement lines
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = c(1:15)) +
  scale_color_manual(values = colorVec) +
  labs(x = x, y = y, shape = "groupId") +
  ggtitle(paste("Genomic --> Morphologic\n",x, " vs ", y,sep=""))+
  theme(legend.position = "none", 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())

})
wrap_plots(plots, ncol = 3)

```
### if we want to make sure you reorder the 

```r
library(clue)

# Compute correlation matrix between dimensions
corrMDSdim <- abs(cor(GenMDS[, c(1:k)], MorfMDS[, c(1:k)]))
# Convert to a cost matrix by subtracting from the maximum value
cost <- max(corrMDSdim) - corrMDSdim

# Solve the Linear Sum Assignment Problem (Hungarian algorithm)
corres <- solve_LSAP(cost)

# Reorder MorfMDS accordingly
MorfMDS_reordered <- MorfMDS[, as.vector(corres)]
names(MorfMDS_reordered)<-paste0("Dim",c(1:k))
PROCreordered <- procrustes(GenMDS[,c(1:k)],MorfMDS_reordered[,c(1:k)])
Genomic <- as.data.frame(PROCreordered$X)
Genomic$dataset <- "Genomic"
Genomic$groupId <- rownames(PROCreordered$X)

Morphologic <- as.data.frame(PROCreordered$Yrot)
Morphologic$dataset <- "ReorderedMorphologic"
Morphologic$groupId <- rownames(PROCreordered$Yrot)
names(Morphologic)[c(1:k)]<-paste0("Dim",c(1:k))

# Combine
plot_df_reordered <- rbind(Genomic, Morphologic)
plot_df_reordered$region = metadata$region[match(plot_df_reordered$groupId, metadata$groupId)]
names(plot_df_reordered)[c(1:k)]<-paste0("Dim",c(1:k))

plot_df_reordered$groupId=as.factor(plot_df_reordered$groupId)

dims <- paste0("Dim",c(1:k))
pairs <- combn(dims, 2, simplify = FALSE)


# Plot with displacement lines
plots <- lapply(pairs, function(p) {
  x <- p[1]
  y <- p[2]
  ggplot(plot_df_reordered, aes_string(x = x, y = y, color = "region",shape="groupId")) +
  geom_path(aes(group = groupId), linewidth = 0.4, alpha = 0.7,arrow=arrow(type = "closed",length=unit(0.1,"inches"))) +   # displacement lines
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = c(1:15)) +
  scale_color_manual(values = colorVec) +
  labs(x = x, y = y, shape = "groupId") +
  ggtitle(paste("Genomic --> Morphologic reordered\n",x, " vs ", y,sep=""))+
  theme(legend.position = "none", 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())

})
wrap_plots(plots, ncol = 3)
```


### Neighbour-joining comparison with tanglegram

```r
require(dendextend)
require(ape)
library(phytools)

GenNJchanged<-force.ultrametric(root(GenNJ,outgroup = "ALASKA",resolve.root = T))
MorfNJchanged<-force.ultrametric(root(MorfNJ,outgroup = "ALASKA",resolve.root = T))
tanglegram(GenNJchanged,MorfNJchanged,
            main_left="Genomic",
            main_right="Morphological",
            sort=FALSE,
            highlight_branches_lwd = FALSE,
            margin_inner=12,rank_branches=F,
            match_order_by_labels=F
)

