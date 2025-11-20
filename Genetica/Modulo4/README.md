# Análisis de genomas antiguos: estadísticos F
En este practico utilizaremos el paquete de R admixtools2 para evaluar la afinidad genética y posible eventos de mestizaje en individuos antiguos de América, particularmente de Patagonia. Para ello se utilizarán los estadísticos outgroup-ƒ3 y ƒ4. Archivos para este módulo puede descargarse aquí:


https://www.dropbox.com/scl/fi/jiehanqgpciljyu1g8rqh/Input_modulo4.zip?rlkey=j8r08c3064sx0wsz52zz2clq6&dl=0



##       Instalación de paquetes necesarios (sólo una vez) en Rstudio

### 1A. Instalar paquete admixtools2
```r
install.packages("devtools") # correr si "devtools" no está instalado 
devtools::install_github("uqrmaie1/admixtools")
```

### 2A: En caso de problemas, instalar manualmente los siguientes paquetes:
```r
install.packages("Rcpp")
install.packages("tidyverse")
install.packages("igraph")
install.packages("plotly") 
```

### 2B: Volver a intentar
```r
devtools::install_github("uqrmaie1/admixtools")
```
### 3A: Otra opción para instalar admixtools2
```r
install.packages("remotes")
remotes::install_github("uqrmaie1/admixtools")
```

#       Análisis: outgroup-ƒ3

### Cargar modulos
```r
library(admixtools)
library(tidyverse)
library(ggplot2)
library(dplyr)
```

### Especificar directorio de trabajo (cambiar según corresponda)
```r
setwd("~/Documentos/Modulo3")
```

### Leer archivo con metadatos
Este archivo contiene información sobre los individuos a analizar

```r
metadf = read.csv("finalSet_metadata.csv", 
                  sep = ",", header = T,na.strings=c(""," ","NA"))
```

### Prefijo de archivos eigenstrat
```r
prefix <- "finalSet"
```

### Caso 1: Evaluar deriva genética compartida entre representante de población ancestral en Beringia e individuos antiguos de América utilizando outgroup-ƒ3

Primero tenemos que definir las poblaciones a utilizar:

```r
popA = c('USA_Ancient_Beringian.SG')
popB = c('USA_Anzick.SG',
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
         'Argentina_BeagleChannel_100BP.SG')
outg = c("Mbuti")        
```

Análisis:
```r
f3res <- f3(prefix, outg, popA, popB,
            outgroupmode = TRUE)
```

Explorar los resultados en Rstudio:

Podemos simplemente evaluar los resultados mirando la tabla con resultados o graficarlos. A continuación dos ejemplos:

1. Visualización clásica:

```r
ggplot(f3res, aes(x = reorder(pop3, est), y = est)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = (est - se*3), ymax = (est + se*3)),
                width = 0.2) +
  coord_flip() +
  labs(x = "",
       y = expression(f[3]),
       title = "f3(Y, AncientBeringia; Mbuti)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
```

<img width="1200" height="900" alt="outf3_v1" src="https://github.com/user-attachments/assets/c6692d67-312f-4fff-979c-3cd6d2974e10" />


2. Mapa de calor

  a) Agregar coordenadas geográficas:

```r
f3res$lat = metadf$latitude[match(f3res$pop3, metadf$popId)]
f3res$lon = metadf$longitude[match(f3res$pop3, metadf$popId)]
```

  b) Graficar:

```r
library(tidyverse)
library(ggthemes)

world <- map_data("world")

latlimits <- c(-56,56) # bottom, top 
longlimits <- c(-140,-35) # left, right

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


g1 + 
  geom_point(data = f3res, aes(x = lon, y = lat, fill = est, stroke=1/5), 
                        shape = 21, color = "black", size = 3) +
  scale_fill_distiller(palette = "RdYlBu",name = "Outgroup-f3")
```


<img width="1000" height="900" alt="outf3_v2" src="https://github.com/user-attachments/assets/ecaf33d8-f782-4077-b726-6840ab74d141" />


## Preguntas y ejercicios adicionales
1. ¿Es posible identificar una tendencia?
2. ¿Hay diferencias significativas entre las afinidades genéticas de Ancient Beringia y los otros individuos antiguos de América?

### Caso 2: Repetir con Anzick es posición de popA (y sacar Anzick de popB).


### Caso 3: Repetir usando sólo poblaciones de Patagonia, distinguiendo entre: individuos del Holoceno Medio (en popA) versus Tardío (en popB)

```r
popA = c('Chile_WesternArchipelago_Ayayema_5100BP.SG','Chile_PuntaSantaAna_7300BP.SG','Argentina_NorthTierradelFuego_LaArcillosa2_5800BP')
popB = c('Chile_WesternArchipelago_800BP.SG',
         'Chile_WesternArchipelago_1200BP.SG',
         'Chile_BeagleChannel_800BP.SG',
         'Argentina_NorthTierradelFuego_500BP',
         'Argentina_BeagleChannel_1500BP',
         'Argentina_MitrePeninsula_400BP',
         'Argentina_BeagleChannel_100BP',
         'Chile_SouthernContinent_400BP',
         'Chile_NorthTierradelFuego_100BP',
         'Chile_StraitOfMagellan_100BP.SG',
         'Argentina_NorthTierradelFuego_100BP.SG',
         'Argentina_BeagleChannel_100BP.SG')
outg = c("Mbuti")      


# Luego de realizar el análisis, se puede graficar de la siguiente manera:


ggplot(f3res, aes(x = reorder(pop3, est), y = est, color = pop2)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = (est - se*3), ymax = (est + se*3)),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("orange2", "steelblue1", "pink3")) +
  coord_flip() +
  labs(x = "Target population",
       y = expression(f[3]),
       color = "Pop2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")
```

¿Qué grupos tardíos tienen mayor afinidad con individuos tempranos? ¿Se observan diferencias significativas?


#       Análisis: d-statistics o f4

### Cargar modulos
```r
library(admixtools)
library(tidyverse)
library(ggplot2)
library(dplyr)

# Especificar directorio de trabajo (cambiar según corresponda)
setwd("~/Dropbox/2025/CONFERENCES/ALAB2025/Workshop/Analysis/3.Practico_Fstatistics/")

# Leer archivo con metadatos. Este archivo contiene información sobre los individuos a analizar

metadf = read.csv("finalSet_metadata.csv", 
                  sep = ",", header = T,na.strings=c(""," ","NA"))

# Prefijo de archivos eigenstrat
prefix <- "finalSet" 
```

### Caso 1: Linajes tempranos
Usando outgroup-ƒ3 evaluamos, por separado, la afinidad genética entre individuos antiguos de América con un linaje temprano de Beringia (USA_Ancient_Beringian.SG) y un linaje temprano de América representado por Anzick, siendo este último el representante mas temprano de un linaje característico de poblaciones en Sudamérica. A continuación, evaluaremos las afinidades genéticas de individuos antiguos de América contra estos dos individuos tempranos utilizando el estadístico f4. ¿Qué expectativas tenemos respecto a estos resultados? 

A recordar: \
D ~ 0, pop1 y pop2 están simétricamente relacionadas a pop3 \
D >> 0, pop1 comparte mas alelos con pop3 que con pop2 \
D << 0, pop2 comparte mas alelos con pop3 que con pop1 \
No sólo es importante el valor del estadístico. El valor de z nos indicará si este valor es significativamente distinto de 0. \

Tenemos entonces que definir la posición de 4 poblaciones. Podemos visualizarlo de la siguiente manera:

<img width="800" height="400" alt="tree_plot" src="https://github.com/user-attachments/assets/d9a5f618-e006-46a4-a2c5-090bfaa388a8" />


En la posición de pop1 usaremos a USA_Ancient_Beringian.SG, en pop2 a USA_Anzick.SG, pop3 será una lista de individuos antiguos de América, mientras que pop4 es un outgroup (Mbuti). 

```r
popC = c('USA_Nevada_SpiritCave_11000BP.SG',
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
         'Argentina_BeagleChannel_100BP.SG')


dstat1 = qpdstat(prefix, pop1="USA_Ancient_Beringian.SG", pop2="USA_Anzick.SG", pop3=popC, pop4="Mbuti")
```

¿Están pop1 y pop2 igualmente relacionados a los individuos de pop3? ¿Se ajusta este resultado a las expectativas?

### Caso 2: Linajes tempranos de Sudamérica
Usaremos ahora los individuos con fechados mas tempranos (~10.000) de América y evaluaremos si presentan afinidades genéticas diferenciales con un set de individuos antiguos de América. En particular, compararemos Anzick con otros individuos tempranos. ¿Qué expectativas tenemos respecto a estos resultados? 

```r
# Definir grupos en pop2 
popB = c('USA_Nevada_SpiritCave_11000BP.SG',
         'Brazil_Sumidouro_10100BP.SG',
         'Brazil_LapaDoSanto_9600BP',
         'Chile_LosRieles_12000BP')

# Definir grupos en pop3
popC = c('Peru_Lauricocha_8600BP',
         'Chile_WesternArchipelago_800BP.SG',
         'Chile_WesternArchipelago_1200BP.SG',
         'Chile_BeagleChannel_800BP.SG',
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
         'Argentina_BeagleChannel_100BP.SG')

# Correr análisis
dstat2 = qpdstat(prefix, pop1="USA_Anzick.SG", pop2=popB, pop3=popC, pop4="Mbuti")

# Figura con resultados:
ggplot(dstat2, aes(x = pop3, y = est, color = pop2)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = (est - se*3), ymax = (est + se*3)),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("orange2", "steelblue1", "pink3","yellow2")) +
  coord_flip() +
  labs(x = "Target population",
       y = " <- Target",
       color = "Pop2") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")
```

¿Qué podemos concluir de este resultado?


### Caso 3: Patagonia
Ahora realizaremos un análisis específico de Patagonia, evaluando si existen diferencias genéticas significativas entre individuos del Holoceno Medio y del Holoceno Tardío en la región. Específicamente, evaluaremos si algún individuo del Holoceno Medio tiene mayor afinidad genética hacia algún individuo/grupo del Holoceno Tardío. 

Los individuos del Holoceno Medio son: \
`'Chile_WesternArchipelago_Ayayema_5100BP.SG', 'Chile_PuntaSantaAna_7300BP.SG', 'Argentina_NorthTierradelFuego_LaArcillosa2_5800BP'`

Los individuos/grupos del Holoceno tardío son: \
`'Chile_WesternArchipelago_800BP.SG',
         'Chile_WesternArchipelago_1200BP.SG',
         'Chile_BeagleChannel_800BP.SG',
         'Argentina_NorthTierradelFuego_500BP',
         'Argentina_BeagleChannel_1500BP',
         'Argentina_MitrePeninsula_400BP',
         'Argentina_BeagleChannel_100BP',
         'Chile_SouthernContinent_400BP',
         'Chile_NorthTierradelFuego_100BP',
         'Chile_StraitOfMagellan_100BP.SG',
         'Argentina_NorthTierradelFuego_100BP.SG',
         'Argentina_BeagleChannel_100BP.SG'`


¿Qué posición ocupará cada uno de estos individuos en el análisis? \
Grafique sus resultados

Podemos también realizar preguntas mas específicas en la región, tales como la relación entre individuos cazadores-recolectores marinos y terrestres, o quienes habitaron la región de los archipiélagos occidentales versus Tierra del Fuego. Proponga configuraciones apropiadas para realizar estos análisis. 








