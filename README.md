# Broad-et-al.-2023

# Data Availability
Gene expression sequencing data files are available on The University of Queensland Research Data Manager repository. The data files generated at each step of this process are also available.

# Bioinformatics
Gene expression sequencing was conducted by the Beijing Genomics Institute (BGI) using the DNBseq platform. Pooled samples were sequenced randomly with respect to Gravitropic/Agravitropic phenotype, control and rotation group and time group (Samples were taken at 5 time points). 

## Quality filtering
Basic quality control of the raw reads for all 60 samples was conducted using FASTQC (SOURCE) on default parameters. Trimmomatic (SOURCE) was used to identify the first ten base pairs of each read as low quality and trim them. A report for all individuals was compiled from the FASTQC outputs using MultiQC (SOURCE). 

## Assembly and Mapping
We assembled and mapped the trimmed reads onto a S. lautus reference transcriptome (SOURCE Wilkinson et al., 2021) using HISAT2 (SOURCE). 

## Pricniple Component Analysis

## Differential Gene Expression

## WGCNA


