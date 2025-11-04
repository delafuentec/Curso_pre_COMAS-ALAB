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

# Cargar modulos
`library(admixtools)`

`library(tidyverse)`

`library(ggplot2)`

`library(dplyr)`


