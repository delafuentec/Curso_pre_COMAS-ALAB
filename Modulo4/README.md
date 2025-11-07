# Análisis de genomas antiguos: estadísticos F
En este practico utilizaremos el paquete de R admixtools2 para evaluar la afinidad genética y posible eventos de mestizaje en individuos antiguos de Patagonia. Para ello se utilizarán los estadísticos outgroup-ƒ3 y ƒ4

##       Instalación de paquetes necesarios (sólo una vez) en Rstudio

### 1A. Instalar paquete admixtools2

`install.packages("devtools")` # correr si "devtools" no está instalado 

`devtools::install_github("uqrmaie1/admixtools")`

### 2A: En caso de problemas, instalar manualmente los siguientes paquetes:
`install.packages("Rcpp")`

`install.packages("tidyverse")`

`install.packages("igraph")`

`install.packages("plotly")`

### 2B: Volver a intentar
`devtools::install_github("uqrmaie1/admixtools")`

### 3A: Otra opción para isntallar admixtools2
`install.packages("remotes")`

`remotes::install_github("uqrmaie1/admixtools")`


##       Análisis

### Cargar modulos
`library(admixtools)`

`library(tidyverse)`

`library(ggplot2)`

`library(dplyr)`

### Especificar directorio de trabajo (cambiar según corresponda)
`setwd("~/Dropbox/2025/CONFERENCES/ALAB2025/Workshop/Analysis/3.Practico_Fstatistics/")`

### Leer archivo con metadatos
Este archivo contiene información sobre los individuos a analizar

`metadf = read.csv("finalSet_metadata.csv", 
                  sep = ",", header = T,na.strings=c(""," ","NA"))`

### Prefijo de archivos eigenstrat
`prefix <- "finalSet" `

### Caso 1: Evaluar deriva genética compartida entre representante de población ancestral en Beringia e individuos antiguos de América utilizando outgroup-ƒ3

Primero tenemos que definir las poblaciones a utilizar:

`popA = c('USA_Ancient_Beringian.SG')`

`popB = c('USA_Anzick.SG',
         'USA_Nevada_SpiritCave_11000BP.SG',
         'Brazil_Sumidouro_10100BP.SG',
         'Brazil_LapaDoSanto_9600BP',
         'Peru_Lauricocha_8600BP',
         'Chile_WesternArchipelago_800BP.SG',
         'Chile_WesternArchipelago_1200BP.SG',
         'Chile_BeagleChannel_800BP.SG',
         'Chile_LosRieles_12000BP',
         'Chile_Conchali_700BP',
         'Chile_WesternArchipelago_Ayayema_5100BP.SG',
         'Chile_PuntaSantaAna_7300BP.SG',
         'Argentina_NorthTierradelFuego_500BP',
         'Argentina_BeagleChannel_1500BP',
         'Argentina_MitrePeninsula_400BP',
         'Argentina_NorthTierradelFuego_LaArcillosa2_5800BP',
         'Argentina_BeagleChannel_100BP',
         'Chile_SouthernContinent_400BP',
         'Chile_NorthTierradelFuego_100BP',
         'Argentina_ArroyoSeco2_7700BP',
         'Chile_StraitOfMagellan_100BP.SG',
         'Argentina_NorthTierradelFuego_100BP.SG',
         'Argentina_BeagleChannel_100BP.SG')`

`outg = c("Mbuti")`        

Análisis:

`f3res <- f3(prefix, outg, popA, popB,
            outgroupmode = TRUE)`

Explorar los resultados en Rstudio:

Podemos simplemente evaluar los resultados mirando la tabla con resultados o graficarlos. A continuación dos ejemplos:

1. Visualización clásica:

`ggplot(f3res, aes(x = reorder(pop3, est), y = est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (est - se*3), ymax = (est + se*3)),
                width = 0.2) +
  coord_flip() +
  labs(x = "",
       y = expression(f[3]),
       title = "f3(Y, AncientBeringia; Mbuti)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))`


<img width="2400" height="1800" alt="outf3_v1" src="https://github.com/user-attachments/assets/c6692d67-312f-4fff-979c-3cd6d2974e10" />


2. Mapa de calor

  a) Agregar coordenadas geográficas:

  `f3res$lat = metadf$latitude[match(f3res$pop3, metadf$popId)]`
  
  `f3res$lon = metadf$longitude[match(f3res$pop3, metadf$popId)]`

  b) Graficar:

`library(tidyverse)`

`library(ggthemes)`

`world <- map_data("world")`

`latlimits <- c(-56,56)` bottom, top

`longlimits <- c(-140,-35)` left, right

`
g1 = ggplot() +
  geom_polygon(data = world,aes(x = long, y = lat, group = group), 
               fill = NA, color = 'black', linewidth = 0.2) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_quickmap(xlim = longlimits, ylim = latlimits) +
  xlab('') + ylab('') +
  theme_map() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        axis.ticks = element_blank())
`


`
g1 + 
  geom_point(data = f3res, aes(x = lon, y = lat, fill = est, stroke=1/5), 
                        shape = 21, color = "black", size = 3) +
  scale_fill_distiller(palette = "RdYlBu",name = "Outgroup-f3")
`


<img width="1459" height="1357" alt="outf3_v2" src="https://github.com/user-attachments/assets/ecaf33d8-f782-4077-b726-6840ab74d141" />


#### Preguntas y ejercicios adicionales
¿Es posible identificar una tendencia?

¿Hay diferencias significativas entre las afinidades genéticas de Ancient Beringia y los otros individuos antiguos de América?

Repetir con Anzick es posición de popA (y sacar Anzick de popB)


