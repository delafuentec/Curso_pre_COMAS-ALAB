# Integración datos genéticos y morfométricos

En esta sección exploraremos la integración entre resultados de análisis genéticos y morfométricos. En primer lugar, estimaremos distancias desde los datos genéticos y morfométricos, para luego:
- hacer analísis descriptivos con ambos y aplicar metódos de comparación de los resultados
- comparar estas distancias y evaluar su posible asociación.

## Análisis genéticos

### Cálculo de las similitudes genéticas entre grupos.
Utilizaremos el paquete de R `admixtools` para estimar el outgroup-ƒ3 entre todos los posibles pares de grupos en la base de datos. Al igual que en el módulo 4 de Genética:

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
setwd("/Documentos/ALAB2025/Integración/")

# Definir input
prefix <- "AmericaByGroups/Genodata_forCompMorphology"

# Definir poblaciones a analizar. Se utilizaran todos los pares posibles de grupos
popA = c('PATAGONIA_MARITIMO',
         'BRASIL_TARDIO',
         'BRASIL_TEMPRANO',
         'CENTRO_PERU',
         'NORTE_PERU',
         'CHANNEL_ISLAND',
         'PATAGONIA_TERRESTRE',
         'PERU_TempranoHoloceno',
         'SUR_PERU',
         'PERU_MedioHoloceno',
         'BAJA_CALIFORNIA_SUR',
         'PAMPA',
         'ALASKA')
popB = c('PATAGONIA_MARITIMO',
         'BRASIL_TARDIO',
         'BRASIL_TEMPRANO',
         'CENTRO_PERU',
         'NORTE_PERU',
         'CHANNEL_ISLAND',
         'PATAGONIA_TERRESTRE',
         'PERU_TempranoHoloceno',
         'SUR_PERU',
         'PERU_MedioHoloceno',
         'BAJA_CALIFORNIA_SUR',
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
distGenMDS <- 1 - distGenMDS

### Sacaremos Alaska (se utilizará como outgroup en NJ)
distGenMDS<-distGenMDS[row.names(distGenMDS) != "ALASKA" , colnames(distGenMDS) != "ALASKA"]

# Hacer análisis de escalamiento multidimensional (MDS)
k=10
GenMDS <- as.data.frame(cmdscale(distGenMDS, k = k))
GenMDS$groupId <- rownames(distGenMDS)
names(GenMDS)[c(1:k)]<-paste0("Dim",c(1:k))

# Agregaremos información de región para facilitar visualización
metadata = data.frame(groupId = c('PATAGONIA_MARITIMO',
                                  'BRASIL_TARDIO',
                                  'BRASIL_TEMPRANO',
                                  'CENTRO_PERU',
                                  'NORTE_PERU',
                                  'CHANNEL_ISLAND',
                                  'PATAGONIA_TERRESTRE',
                                  'PERU_TempranoHoloceno',
                                  'SUR_PERU',
                                  'PERU_MedioHoloceno',
                                  'BAJA_CALIFORNIA_SUR',
                                  'PAMPA'),
                      region = c('Patagonia',
                                 'Brasil',
                                 'Brasil',
                                 'Andes Centro-Sur',
                                 'Andes Centro-Sur',
                                 'Norteamerica',
                                 'Patagonia',
                                 'Andes Centro-Sur',
                                 'Andes Centro-Sur',
                                 'Andes Centro-Sur',
                                 'Mesoamerica',
                                 'Pampa'),
                      region_color = c("steelblue2","mediumseagreen","mediumseagreen","orange3","orange3",
                                       "deeppink","steelblue2","orange3","orange3","orange3",
                                       "darkgrey","rosybrown2"),
                      group_shape = c(1:12))

GenMDS$region = metadata$region[match(GenMDS$groupId, metadata$groupId)]
GenMDS$region_color = metadata$region_color[match(GenMDS$groupId, metadata$groupId)]
GenMDS$group_shape = metadata$group_shape[match(GenMDS$groupId, metadata$groupId)]

GenMDS$region  <- factor(GenMDS$region, levels = c("Norteamerica","Mesoamerica","Andes Centro-Sur", "Brasil","Pampa","Patagonia"))

legend_colors <- setNames(metadata$region_color, metadata$groupId)
legend_shapes <- setNames(metadata$group_shape, metadata$groupId)

# Grafico de MDS
g1 = ggplot(GenMDS, aes(x = Dim1, y = Dim2, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim1', y = 'Dim2') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))   

g2 = ggplot(GenMDS, aes(x = Dim3, y = Dim4, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim3', y = 'Dim4') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))   


g3 = ggplot(GenMDS, aes(x = Dim5, y = Dim6, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim5', y = 'Dim6') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))       


g1 + g2 + g3 + 
  plot_annotation(
    title = "Genetica",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
ggsave("mds_america.png", width = 15, height = 3)

```

### Árbol de Neighbor-joining con distancias genéticas

```r

require(ape)
# Construir matriz de distancias genéticas
distGenNJ <- acast(f3res, pop2 ~ pop3, value.var = "est")
distGenNJ <- 1/distGenNJ
diag(distGenNJ) <- 0
distGenNJ<-distGenNJ/max(distGenNJ)
distGenNJ<-as.dist(distGenNJ)
GenNJ<-bionj(distGenNJ)
GenNJ<-root(GenNJ,outgroup = "ALASKA",resolve.root = T)

if(sum(! GenNJ$tip.label %in% row.names(metadata) )>0){
  print(paste(GenNJ$tip.label[! GenNJ$tip.label %in% row.names(metadata)]))
  stop("pb groups not in your metadata")
}

write.tree(GenNJ,file="NJwithGenomicData.nwk",append = F,tree.names = F)
GenNJplot<-read.tree("NJwithGenomicData.nwk")
dicoPLOT <- metadata[GenNJplot$tip.label, ]

dicoPLOT$Color <- metadata$region_color[ match(GenNJplot$tip.label, metadata$groupId) ]
dicoPLOT$pch_base <- as.numeric(metadata$group_shape[ match(GenNJplot$tip.label,
                                                            metadata$groupId)])

# Hacer figura
pdf("NJ_with_Genomic_Data.pdf", width = 10, height = 10)
plotStats <- plot(
  GenNJplot,
  tip.color = dicoPLOT$Color,
  cex = 1,
  use.edge.length = FALSE,
  show.tip.label = TRUE,
  label.offset = 5,
  align.tip.label = TRUE,
  main = "NJ with genomic data")

points(
  x = rep(nrow(dicoPLOT) + 2.5, nrow(dicoPLOT)),
  y = 1:nrow(dicoPLOT),
  cex = 1.2,
  pch = dicoPLOT$pch_base,
  col = dicoPLOT$Color)
dev.off()

```

## Análisis morfólogicos

### Cálculo de las distancias morfólogicas entre grupos.

```r
datosMorf<-read.table("AmericaByGroups/MorphologicMeansPerGroup.csv",stringsAsFactors=F,header=T,sep="\t")
row.names(datosMorf) <- datosMorf[,1]
datosMorf <- datosMorf[,-1]
datosMorf <- as.matrix(datosMorf)

## función de normalización
normalize_lines <- function(v){ 
  return(v / sqrt(sum(v^2)))
}

datosMorf_norm <- t(apply(datosMorf, 1, normalize_lines))

##Distancias Procrustes
n <- nrow(datosMorf_norm)
distMorf <- matrix(0, n, n)
row.names(distMorf) <- colnames(distMorf) <- row.names(datosMorf_norm)

for(i in 1:n){
  for(j in 1:n){
    if(i < j){
      distMorf[i,j] <- distMorf[j,i] <- sqrt(sum((datosMorf_norm[i,] - datosMorf_norm[j,])^2))
    }
  }
}

### Guardar matriz
write.table(distMorf, "AmericaByGroups/TablasDistancias/MorphologicalDistances.tsv",col.names=T,row.names=T,sep="\t")
```

### MDS con distancias morfólogicas
```r

distMorf<-distMorf/max(distMorf)
distMorfMDS<-distMorf[row.names(distMorf) != "ALASKA" , colnames(distMorf) != "ALASKA"]

# Hacer análisis de escalamiento multidimensional (MDS)
k=10
MorfMDS <- as.data.frame(cmdscale(distMorfMDS, k = k))
MorfMDS$groupId <- rownames(distMorfMDS)

MorfMDS$region = metadata$region[match(MorfMDS$groupId, metadata$groupId)]

row.names(MorfMDS)<-MorfMDS$groupId
MorfMDS<-MorfMDS[ GenMDS$groupId,]
names(MorfMDS)[c(1:k)]<-paste0("Dim",c(1:k))

# Grafico de MDS
g1 = ggplot(MorfMDS, aes(x = Dim1, y = Dim2, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim1', y = 'Dim2') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))   

g2 = ggplot(MorfMDS, aes(x = Dim3, y = Dim4, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim3', y = 'Dim4') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))   


g3 = ggplot(MorfMDS, aes(x = Dim5, y = Dim6, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim5', y = 'Dim6') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))       


g1 + g2 + g3 + 
  plot_annotation(
    title = "Morfología",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
ggsave("mds_america_morfologia.png", width = 15, height = 3)

```


### Árbol de Neighbor-joining desde distancias morfólogicas

```r
# Construir matriz de distancias genéticas
distMorfNJ <- as.dist(distMorf)

MorfNJ<-bionj(distMorfNJ)
MorfNJ<-root(MorfNJ,outgroup = "ALASKA",resolve.root = T)
if(sum(! MorfNJ$tip.label %in% row.names(metadata) )>0){
  print(paste(MorfNJ$tip.label[! MorfNJ$tip.label %in% row.names(metadata)]))
  stop("pb groups not in your metadata")
}
write.tree(MorfNJ,file="NJwithMorphologicalData.nwk",append = F,tree.names = F)

MorfNJ<-read.tree("NJwithMorphologicalData.nwk")
dicoPLOT <- metadata[MorfNJ$tip.label, ]

dicoPLOT$Color <- metadata$region_color[ match(MorfNJ$tip.label, metadata$groupId) ]
dicoPLOT$pch_base <- as.numeric(metadata$group_shape[ match(MorfNJ$tip.label,
                                                            metadata$groupId)])

# Hacer figura
pdf("NJ_with_Morpho_Data.pdf", width = 10, height = 10)
plotStats <- plot(
  MorfNJ,
  tip.color = dicoPLOT$Color,
  cex = 1,
  use.edge.length = FALSE,
  show.tip.label = TRUE,
  label.offset = 5,
  align.tip.label = TRUE,
  main = "NJ with genomic data")

points(
  x = rep(nrow(dicoPLOT) + 2.5, nrow(dicoPLOT)),
  y = 1:nrow(dicoPLOT),
  cex = 1.2,
  pch = dicoPLOT$pch_base,
  col = dicoPLOT$Color)
dev.off()

```

## Comparaciones de datos morfológicos y genéticos

En esta sección realizaremos un análisis de carácter exploratorio, cuyo objetivo es tratar de responder a la necesidad de dialogo entre las disciplinas. Más allá de articular conclusiones desde diferentes evidencias, pretendemos iniciar una discusión epistemiológica para articular los datos de diferente índole en un análisis interdisciplinar.

### MDS analysis with procrustes

Anteriormente realizamos un MDS de datos genéticos y otro de datos morfológicos. En primera instancia, encontraremos la mejor correspondencia entre dimensiones de los distintos MDS (es decir, Dim1 de GenMDS tiene una mejor correspondencia con DimXX de MorfMDS, etc).

```r

#install.package('vegan') # si fuese necesario
require(vegan)

PROC <- procrustes(GenMDS[,c(1:k)],MorfMDS[,c(1:k)])
Genomic <- as.data.frame(PROC$X)
Genomic$dataset <- "Genomic"
Genomic$groupId <- rownames(PROC$X)

Morphologic <- as.data.frame(PROC$Yrot)
Morphologic$dataset <- "Morphologic"
Morphologic$groupId <- rownames(PROC$Yrot)
names(Morphologic)[c(1:k)]<-paste0("Dim",c(1:k))

# Combinar
plot_df <- rbind(Genomic, Morphologic)
plot_df$region = metadata$region[match(plot_df$groupId, metadata$groupId)]
names(plot_df)[c(1:k)]<-paste0("Dim",c(1:k))

plot_df$groupId=as.factor(plot_df$groupId)

# Graficar
g1 = ggplot(plot_df, aes(x = Dim1, y = Dim2, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  geom_path(aes(group = groupId, color = groupId), linewidth = 0.4, alpha = 0.7,
            arrow=arrow(type = "closed",length=unit(0.1,"inches")),show.legend = FALSE) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim1', y = 'Dim2') +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))    

g2 = ggplot(plot_df, aes(x = Dim3, y = Dim4, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  geom_path(aes(group = groupId, color = groupId), linewidth = 0.4, alpha = 0.7,
            arrow=arrow(type = "closed",length=unit(0.1,"inches")),show.legend = FALSE) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim3', y = 'Dim4') +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))   


g3 = ggplot(plot_df, aes(x = Dim5, y = Dim6, color = groupId, shape=groupId)) +
  geom_point(size = 3) +
  geom_path(aes(group = groupId, color = groupId), linewidth = 0.4, alpha = 0.7,
            arrow=arrow(type = "closed",length=unit(0.1,"inches")),show.legend = FALSE) +
  theme_bw() +
  scale_shape_manual(values = legend_shapes, name = 'Grupo') +
  scale_color_manual(values = legend_colors, name = 'Grupo') +
  labs(x = 'Dim5', y = 'Dim6') +
  theme(legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing = unit(0.3, "lines"))     


g1 + g2 + g3 + 
  plot_annotation(
    title = "Genomic --> Morphologic",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )
ggsave("mds_america_genomics-morfologia.png", width = 15, height = 5)

```

### Neighbour-joining comparison with tanglegram

```r
require(dendextend)
require(ape)
library(phytools)
library(colorspace)

GenNJchanged<-force.ultrametric(root(GenNJ,outgroup = "ALASKA",resolve.root = T))
MorfNJchanged<-force.ultrametric(root(MorfNJ,outgroup = "ALASKA",resolve.root = T))

#Graficar y guardar PDF
pdf("NJ_Genomics_vs_Morpho.pdf", width = 10, height = 10)
par(oma = c(0, 0, 4, 0)) 
plot(cophylo(GenNJchanged, MorfNJchanged), lwd=2, fsize=0.8, col="black", dcol="black", lty=c(1,1))
mtext(at = -0.25, side=3, text = "Genomic", xpd = TRUE)
mtext(at = 0.25, side=3, text = "Morphology", xpd = TRUE)
dev.off()

```

### Poner a prueba la asociación entre matrices de distancias morfológicas y genéticas.

#### Test de Mantel
Este método es aceptable para un análisis exploratorio, pero no recomendable como prueba principal, ya que sufre de varias limitaciones:
  - Alta inflación de falsos positivos (problemas de autocorrelación espacial en los datos)
  - No distingue bien entre correlación directa y correlación inducida por una variable tercera (geografía, estructura)
  - Baja potencia estadística.

```r
require(vegan)

distGenMDS<-distGenMDS[ colnames(distMorfMDS),colnames(distMorfMDS)]
mantel_res <- mantel(as.dist(distMorfMDS), as.dist(distGenMDS), method = "pearson", permutations = 99,)

print(mantel_res)

```

### Evaluar si los patrones de variación morfología están asociada a la variación genética

Para esto vamos a realizar un dbRDA.
> El análisis de redundancia (RDA) es una ordenación restringida que utiliza distancias euclidianas. <br>
  > El análisis de redundancia basado en la distancia (dbRDA por su sigla en inglés; Legendre y Anderson, 1999) es una forma más general de realizar el mismo tipo de ordenación restringida.<br>
  > Su generalidad se debe a que puede aplicarse a cualquier matriz de distancia o disimilitud. <br>
  > La conexión entre el RDA y el dbRDA es simplemente el paso de expresar una matriz de disimilitud no euclidiana en un espacio euclidiano (lo que hace un PCoA). 
> Por lo tanto el procedimiento es:
  > 1. Aplicar un PCoA a la matriz de distancia para explicar (en nuestro caso distancias morfológicas) y mantener todas las coordenadas principales. <br>
  > 2. Realizar un RDA sobre las coordenadas principales. Es decir: <br>
  >  2.a. Las coordenadas principales (de la variable para explicar) se someten a una regresión lineal múltiple multivariada (variables explicativas siendo los ejes de la variable explicativa, en nuestro caso del PCoA a partir de la distancias genéticas). <br>
>  2.b. Los valores ajustados y los residuos se someten a PCoA separadas. En este procedimiento se maximiza el ajuste entre la matriz de distancias para explicar que contiene las coordenadas PCoA (entonces, para distancias morfológicas) y la variable explicativa (es decir, del PCoA para distancias genéticas).


Para hacer este análisis hay que hacer un PCoA (muy similar al MDS), con la variable explicativa (distancias genéticas). Guardaremos solo las dimensiones que, acumuladas, explican hasta el 70% de la varianza.

```r
print("PCoA on genetic distances")
pco_gen <- pcoa(as.dist(distGenMDS))
gen_scores <- pco_gen$vectors
Cumul_eig_gen <- pco_gen$values$Cumul_eig
keep_gen<-which(Cumul_eig_gen<0.4)
gen_scores_k <- gen_scores[, keep_gen]


print("now dbRDA")
dbrda_Morf <- capscale(as.dist(distMorfMDS) ~ gen_scores_k, add = TRUE)

print("test significance")
anova_all_Morf<-anova(dbrda_Morf, permutations = 99)
print(anova_all_Morf)
```

> En un análisis dbRDA, si hay una sola variable explicativa, con `anova()`, se evalúa **si el predictor (genética) explica una proporción significativa de la variación morfológica.**.
> Si el valor *p* es menor de 0.05: 
  - El modelo completo explica la morfología significativamente mejor que una asociación aleatoria.  

Ahora bien, la distribución geografía y la temporalidad  pueden ser un factores confundentes.
Aunque las observaciones para un grupo (p.ej. Brasil Tardío) no se obtuvieron para individuos de los mismos sitios, es decir con las mismas coordenadas geográficas ni temporalidad, para el propósito de este ejercicio vamos a considerar que así fue. 
procedimiento:
  - La temporalidad se codifica de 1 a 3 entre los grandes períodos (Holoceno Temprano, Medio y Tardío). Ojo, esta aproximación va a aplanar mucho las distancias temporales entre grupos con tan solo 3 valores posibles 0 (mismo periodo), 1 (Holoceno medio con Holoceno tardío o temprano), 2 (Holoceno temprano vs tardío).
  - Asignamos para las coordenadas geográficas, las de un punto arbitrario en la región considerada.

> Cuando se hace un dbRDA con varias variables explicativas, el procedimiento es igual a lo descrito previamente, solo que en el segundo paso se hace la regresión multivariada en todos los ejes de los PCoA de las diferentes variables explicativas.<br>
  > Es decir, en nuestra caso, se hará la regresión multivariada de los ejes del PCoA de las distancias morfológicas contra los ejes de los PCaA para las distancias genéticas, temporales y geográficas.  

Vamos a construir matrices de distancias geograficas y temporales:

```r
### Leer metadata
meta<-read.table("AmericaByGroups/SummaryTable_PublishedData.tsv",stringsAsFactors = F,header=T,sep="\t")
## info for temporality and geography are repeated across Genetic and Morphological data
meta<-meta[ meta$Tipo.de.dato == "Morfológico",]
## we don't analyze Alaska
meta<-meta[ meta$Región != "ALASKA",]


### Asignamos un valor numérico por rango temporal:
meta$TEMP<-as.numeric(ifelse(meta$Temporalidad=="LH",3,
                             ifelse(meta$Temporalidad=="MH",2,
                                    ifelse(meta$Temporalidad=="EH",1,"UPS???"))))
distGEO<-matrix(0,nrow(meta),nrow(meta))
row.names(distGEO)<-colnames(distGEO)<-meta$Región
distTEMP<-distGEO

for(i in c(1:nrow(distGEO))){
  for(j in c(1:nrow(distGEO))){
    distTEMP[i,j]=abs(meta$TEMP[i]-meta$TEMP[j])
    distGEO[i,j]=sqrt((meta$Latitud[i]-meta$Latitud[j])^2 + 
                        (meta$Longitud[i]-meta$Longitud[j])^2)
  }
}

distGEO<-distGEO/max(distGEO)
distTEMP<-distTEMP/max(distTEMP)
```

Realizaremos el PCoA para estos dos variables (temporalidad / geografia):

```r
print("PCoA on geograhic distances")
pco_geo <- pcoa(as.dist(distGEO))
pco_geo
##Only 2 dimensions catch the whole variance, with the first ~ 90%, let's keep both
geo_scores <- pco_geo$vectors

print("PCoA on temporal distances")
pco_temp <- pcoa(as.dist(distTEMP))
pco_temp
##Only 1 dimension catch the whole variance, with the first ~ 90%, let's keep both
temp_scores <- pco_temp$vectors
```

Ahora podemos volver a hacer el dbRDA incorporando estas dos variables explicativas:
```r
print("dbRDA with Genetics+Geography+Temporality")
dbrda_Morf_covar <- capscale(as.dist(distMorfMDS) ~ gen_scores_k + temp_scores + geo_scores, add = TRUE)

print("test significance")
anova_all_Morf_covar<-anova(dbrda_Morf_covar, permutations = 99)
print(anova_all_Morf_covar)
```

Y probar si este modelo con las variables temporales y geográficas es mejor que el modelo sin:
```r
print("compare models")
anova(dbrda_Morf, dbrda_Morf_covar, permutations = 999)
```

Ahora vamos a graficar las 2 primeras dimensiones del dbRDA en dos situaciones: 1) cuando la variable explicativa en el PCoA basado en distancias genéticas solamente, y 2) cuando también consideramos las distancias temporales y geográficas.

```r
plot(dbrda_Morf, main="Morphology ~ Genetics")
plot(dbrda_Morf_covar, main="Morphology ~ Genetics + Temporality + Geography")
```



