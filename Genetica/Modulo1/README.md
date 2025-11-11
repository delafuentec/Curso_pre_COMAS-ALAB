# Tutorial de Galaxy para filtrar y alinear datos de NGS antiguos

El objetivo de este tutorial es realizar el procedimiento básico de datos de secuenciación NGS de individuos antiguos, incluyendo filtros de calidad, alineamiento, estimación de porcentaje de ADN endógeno y llamado de variantes (haploide). \
Se utilizarán datos de tres muestras, las cuales tienen distintos valores de ADN endógeno. Para el propósito de este tutorial, los datos se limitaron a ADN mitocondrial.

El procesamiento de los datos incluirá los siguientes pasos:

### 1. Subir archivos a galaxy:
####  a) Referencia: archivo fasta
`NC_012920_rCRS.fasta`

Descargado en: https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1

####  b) Archivos de muestras: 3 muestras, formato fastq.gz

`sample1_E1_L1_P1.fastq.gz`\
`sample2_E1_L1_P1.fastq.gz` \
`sample3_E1_L1_P1.fastq.gz`

Estos archivos están disponibles aquí: 
https://www.dropbox.com/scl/fi/5dl822pp8mv81sher4d44/Input_modulo1.zip?rlkey=33p8nprde64c05qq9703z925w&dl=0


### 2. FastQC: control de calidad y filtros iniciales  
Este paso nos permite evaluar la calidad de la secuenciación, número de reads (fragmentos o lecturas), largo de reads, presencia de adaptadores, etc.

### 3. Cutadapt
Se utilizará la herramienta Cutadapt para:

  -Eliminar secuencias de adaptadores de reads \
  -Eliminar Ns \
  -Eliminar fragmentos menores a 25pb \

##### Secuencias de adaptadores (Illumina)

Read 1 Adapter
  +Insert 3' (End) Adapters: Adapter1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
  +Insert 5' (End) Adapters:
  Adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT


### 4. Alineamiento
Se utilizará `BWA aln` para realizar el alineamiento de los reads al genoma de referencia. Luego, se utilizará `samtools` para filtrar reads según calidad y finalmente `Picard MarkDuplicates` para eliminar duplicados de PCR.

### 5. Llamado de variantes
Al trabajar sólo ADN mitocondrial, el llamado de variantes de puede realizar utilizando un procedimiento simple con `bcftools mpileup` seguido de `bcftools call`.

El resultado es un archivo en formato VCF. Adicionalmente se generará un archivo fasta de consenso utilizando `bcftools consensus`

### 6. Bonus: Identificación de haplogrupo mitocondrial
El archivo VCF generado en paso 5 puede ser analizado utilizando la plataforma haplogrep (https://haplogrep.i-med.ac.at/), lo que nos permite hacer una primera caracterización de haplogrupos mitocondriales. 

## Actividades

Se trabajará en la plataforma gratuita Galaxy: https://galaxy-main.usegalaxy.org/

1. Crearse una cuenta. 
2. Subir archivos fasta y fastq. 
3. Acceder al workflow para analizar datos.
4. Correr workflow para las 3 muestras compartidas.

¿Cuantos reads tiene cada muestra? \
¿Tienen adaptadores? \
¿Cuantos reads alinearon a la referencia? \
¿Cuantos reads alineados únicos tienen? \
¿Cuál es el porcentaje de ADN endógeno en las muestras? \
¿Cuál es el linaje de mtDNA en cada muestra? 

5. Descargar archivos VCF, resultado del llamado de variantes, y subir a la plaforma haplogrep para la caracterización de haplogrupos mitocondriales. 

6. Descargar dos archivos fastq.gz de la web y procesar.

Ejemplos:\
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/008/ERR4352398/ERR4352398.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR435/000/ERR4352390/ERR4352390.fastq.gz


¿Tienen adaptadores? \
¿Cuantos reads tienen? \
¿Cuantos reads alinearon a la referencia? \
¿Cuantos reads alineados únicos tienen? \
¿Cuál es el porcentaje de ADN endógeno en las muestras? \
¿Cuál es el linaje de mtDNA en cada muestra? 


