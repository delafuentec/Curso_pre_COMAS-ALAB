# Datos de genotipos: Control de Calidad y Parentesco

## Objetivos
En esta sesión, vamos a:
- repasar la organización de datos de genotipos.
- realizar unos pasos de control de calidad.
- estimar el grado de parentesco entre individuos.
- generar datos listos para llevar a cabo análisis de genética de poblaciones.

## Los datos
Vamos a analizar los datos de genotipos de individuos antiguos y modernos del Cono Sur de Ámerica, y particularmente del Sur de Patagonia.
Existen muchos formatos para este tipo de dato. Vamos a enfocarnos en dos de los más comunes: [</em>plink binary</em>](https://www.cog-genomics.org/plink2/input) y [<em>eigenstrat</em>](https://reich.hms.harvard.edu/software/InputFileFormats).<br>
Ambos formatos, se constituyen de 3 ficheros: 
- uno con la matriz de <strong>genotipos</strong> que puede ser binario, es decir no lisible como texto (<strong>.bed</strong> para <em>plink</em> o <strong>.geno</strong> para <em>eigenstrat</em>).
- uno con la información de las <strong>variantes</strong> analizadas (<strong>.bim</strong> para <em>plink</em> y <strong>.snp</strong> para <em>eigenstrat</em> ).
- uno con la información de los <strong>individuos</strong> analizados (<strong>.fam</strong> para <em>plink</em> y <strong>.ind</strong> para <em>eigenstrat</em>).

Noten que tienen formatos asociados muy similares, y para los cuales la matriz de genotipos es un fichero de texto.

### Fichero de anotación de las variantes
Los ficheros <strong>.bim</strong> y <strong>.snp</strong> contienen una linea por variante y 6 columnas con la misma información. La única diferencia es el orden de las 2 primeras columnas:
- ID de la variante: columna 1 para .snp y 2 para .bim
- cromosoma: columna 2 para .snp y 1 para.bim
- posición genética: columna 3. Se trata de la posición de la variante en función de los eventos de recombinación pasados. Para este curso, no es relevante.
- posición física: columna 4. Se trata del par de base en cual se encuentra la variante en la secuencia del cromosoma.
- Alelo 1: columna 5. Suele ser el alelo de menor frecuencia en la muestra analizada
- Alelo 2: columna 6. Ser ser el alelo de mayor frecuencia en la muestra analizada.

Noten entonces que ambos formatos son compatibles solamente para variantes bialelicas.




Existen varios paquetes de R que permiten leer y procesar estos formatos, pero vamos a usar <huge>XXX</huge>.






