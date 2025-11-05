# Datos de genotipos: Control de Calidad y Parentesco

## Objetivos
En esta sesión, vamos a:
- repasar la organización de datos de genotipos
- realizar unos pasos de control de calidad
- estimar el grado de parentesco entre individuos
- generar datos listos para llevar a cabo análisis de genética de poblaciones.

## Los datos
<p>Vamos a analizar los datos de genotipos de individuos antiguos y modernos del Cono Sur de Ámerica, y particularmente del Sur de Patagonia.</p>
<p>Existen muchos formatos para este tipo de dato. Vamos a enfocarnos en dos de los más comunes: [*plink binary*]{https://www.cog-genomics.org/plink2/input} y [*eigenstrat*]{https://reich.hms.harvard.edu/software/InputFileFormats}.<br>
Ambos formatos, se constituyen de 3 ficheros: 
- uno con la información de las **variantes** analizadas (**.bim** para *plink* y **.snp** para *eigenstrat* )
- uno con la información de los **individuos** analizados (**.fam** para *plink* y **.ind** para *eigenstrat*)
- uno con la matriz de **genotipos** que puede ser binario, es decir no lisible como texto (**.bed** para *plink* o **.geno** para *eigenstrat*).
Noten que ambos formatos *plink* y *eigenstrat* tienen formatos similares para los cuales la matriz de genotipos es un fichero de texto.



