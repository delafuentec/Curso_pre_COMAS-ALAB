# Datos de genotipos: Control de Calidad y Parentesco

## Objetivos
En esta sesión, vamos a:
- repasar la organización de datos de genotipos.
- realizar unos pasos de control de calidad.
- estimar el grado de parentesco entre individuos.
- generar datos listos para llevar a cabo análisis de genética de poblaciones.

## Los datos
<p>Vamos a analizar los datos de genotipos de individuos antiguos y modernos del Cono Sur de Ámerica, y particularmente del Sur de Patagonia.</p>
Existen muchos formatos para este tipo de dato. Vamos a enfocarnos en dos de los más comunes: [<em>plink binary</em>]{https://www.cog-genomics.org/plink2/input} y [<em>eigenstrat</em>]{https://reich.hms.harvard.edu/software/InputFileFormats}.<br>
Ambos formatos, se constituyen de 3 ficheros: 
- uno con la matriz de <strong>genotipos</strong> que puede ser binario, es decir no lisible como texto (<strong>.bed</strong> para <em>plink</em> o <strong>.geno</strong> para <em>eigenstrat</em>).
- uno con la información de las <strong>variantes</strong> analizadas (<strong>.bim</strong> para <em>plink</em> y <strong>.snp</strong> para <em>eigenstrat</em> ).
- uno con la información de los <strong>individuos</strong> analizados (<strong>.fam</strong> para <em>plink</em> y <strong>.ind</strong> para <em>eigenstrat</em>).

Noten que ambos formatos <em>plink</em> y <em>eigenstrat</em> tienen formatos similares para los cuales la matriz de genotipos es un fichero de texto.



