AGREGAR INTRODUCCION




## Generación de matriz de distancias genéticas

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
setwd("~/Dropbox/2025/CONFERENCES/ALAB2025/Workshop/Analysis/Genodata_forCompMorphology/")

# Definir input
prefix <- "Genodata_forCompMorphology"

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
         'PERICUES')
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
         'PERICUES')
outg = c("Mbuti")

# Correr f3
f3res <- f3(prefix, outg, popA, popB,
            outgroupmode = TRUE)
# Filtrar resultados en donde pop2 es la misma que pop3
f3res = f3res %>% filter(pop2 != pop3)

```

### Estimación de distancias y cálculo de MDS

```r
# Construir matriz de distancias genéticas
mdf <- acast(f3res, pop2 ~ pop3, value.var = "est")
diag(mdf) <- 1
mdf <- 1 - mdf / max(mdf, na.rm = TRUE)

# Hacer análisis de escalamiento multidimensional (MDS)
mdf2 <- as.data.frame(cmdscale(mdf, k = 5))
mdf2$groupId <- rownames(mdf2)

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
                                  'PERICUES'),
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
                                 'Mesoamerica'))

mdf2$region = metadata$region[match(mdf2$groupId, metadata$groupId)]

# Grafico de MDS
$ Graficaremos todos los posibles pares de dimensiones
dims <- paste0("V", 1:5)
pairs <- combn(dims, 2, simplify = FALSE)

plots <- lapply(pairs, function(p) {
  x <- p[1]
  y <- p[2]
  ggplot(mdf2, aes_string(x = x, y = y, color = "region")) +
    geom_point() +
    geom_text_repel(aes(label = groupId), color = "black", size = 2.5) +
    theme_bw() +
    scale_color_manual(values = c("orange2", "steelblue1", "pink3","#F0E442","#009E73",'#D55E00')) +
    ggtitle(paste(x, "vs", y)) +
    theme(legend.position = "none", plot.title = element_text(size = 10))
})

# Combinar plots
wrap_plots(plots, ncol = 3)
ggsave("mds_america.png", width = 15, height = 15)
```


### Estimación de distancias y calculo de árbol de Neighbor-joining

[convertir out-f3 a distancia]

[generar matriz]

[calcular NJ]

[Figura]


## Generación de matriz de distancias morfológicas
