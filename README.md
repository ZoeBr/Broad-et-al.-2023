# Broad-et-al.-2023

# Data Availability
Gene expression sequencing data files are available on The University of Queensland Research Data Manager repository. The data files generated at each step of this process are also available.

# Bioinformatics
Gene expression sequencing was conducted by the Beijing Genomics Institute (BGI) using the DNBseq platform. Pooled samples were sequenced randomly with respect to Gravitropic/Agravitropic phenotype, control and rotation group and time group (Samples were taken at 5 time points). 

## Quality filtering
Basic quality control of the raw reads for all 60 samples was conducted using FASTQC (SOURCE) on default parameters. Trimmomatic (SOURCE) was used to identify the first ten base pairs of each read as low quality and trim them. A report for all individuals was compiled from the FASTQC outputs using MultiQC (SOURCE). 

## Assembly and Mapping
We assembled and mapped the trimmed reads onto a S. lautus reference transcriptome (SOURCE Wilkinson et al., 2021) using HISAT2 (SOURCE). The SAM files from HISAT2 were then converted to BAM files using Samtools (SOURCE) and the reads were then processed using TPM calculator (Source) to normalise read counts. 

The aligned reads were then translated into GFF3 format using GMAP/GSNAP (SOURCE).

## Data Filtering
We eliminated the 176,287 transcripts with missing read counts across 60 pools. We then set a read count threshold, filtering out transcripts with fewer than 10 reads in any pool to ensure robust differential gene expression identification. To manage potential alternative splicing of genes, we identified transcript sequences with a high degree of overlap (>95% sequence intersection) and clustered them together. Following these data filtering steps, we were left with 269,209 gene sequences that we could confidently use for further data analysis.

### Gene Expression Quantification and Analysis
Changes in transcript abundance were evaluated using a log 2-fold change analysis. We then used Cuffdiff2 (SOURCE) to identify differentially expressed genes. We used a false discovery rate (FDR) of 0.05% to account for multiple testing in combination with a 2-fold change cutoff.  

We further identified DE genes using the limma trend and Voom functions (https://-bioconductor.org/packages/release/bioc/html/limma.html). Limma is an R-Bioconductor software package that provides an integrated solution for differential expression analyses of data from gene expression experiments. Using limma Trend and Voom, we fitted 2 models to the data to allow us to identify genes that are either differentially expressed between two groups across all time points, or alternatively between two groups at a single specified time point. These models were: (1) log CPM ~ 0 + group + time + Batch and (2) log CPM ~ 0 + group*time + Batch, or 0 + group + time + group:time + Batch.  Here "group" represents the four phenotype/treatment combinations: agravitropic control, agravitropic rotated, gravitropic control, gravitropic rotated. We treated all variables, including time, as discrete factors.

### Gene Set and Gene Ontology Enrichment
We performed ortholog mapping against the Arabidopsis thaliana (Arabidopsis) genome. We used ShinyGO v0.77 (http://bioinformatics.sdstate.edu/go/) for GO term enrichment and visualisation, setting a FDR threshold of 10% (Ge et al., 2019). This software creates networks of GO terms, where nodes represent the enriched GO terms, and node size corresponds to the number of genes within that term. The edges (lines between nodes) indicate connections between GO terms, reflecting the multifunctional nature of many genes. The weight of the edge demonstrates the number of genes involved in both connected GO terms. We proceeded to select the top 10 GO terms for further in-depth analysis.

## Weighted Gene Co-expression Network Analysis
We employed Weighted Gene Co-expression Network Analysis (WGCNA) (http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA) to identify modules of co-expressed genes associated with gravitropism (Langfelder and Horvath, 2008). For this, we used the gravitropic control vs gravitropic rotated differentially expressed (DE) gene subset, as changes in the expression of these genes are expected to indicate a functional response to gravitropism. The aim was to observe alterations in co-expression patterns transitioning from functional responses (typical of gravitropic plants) to non-functional or absent responses (as seen in agravitropic and un-rotated control plants, respectively). To this end, we constructed a co-expression network on the selected genes for each of the 4 groups. WGCNA first constructs a similarity co-expression matrix using Pearson’s correlations, called the correlation matrix. Then an adjacency matrix is created. Next, the topological overlap matrix (TOM) is constructed from the adjacency matrix to evaluate the shared connections in the network, offering a robust method to discern biologically meaningful clusters of co-expressed genes (Langfelder and Horvath, 2008). Co-expression analysis was based on variance stabilised expression using the deseq2 library (varianceStabilizingTransformation function). Genes with >0.99 correlation across all the groups were merged, reducing count from 517 to 428 genes. We constructed unsigned networks with power 12, that is, the adjacency score between two genes was |corr|^12, where corr is the Pearson correlation of the variance stabilised expression across the 15 data points in the given group. Using Cytoscape version 3.8.2 (http://cytoscape.org), we visualised the resulting networks by exporting the edge and node results from WGCNA through the ExportNetworkToCytoscape function. Cytoscape, a powerful open-source software platform widely used in computational biology, visualises complex networks (Su et al., 2014). In our study, Cytoscape facilitated the creation and interpretation of intricate gene interaction networks, showcasing the relationships among differentially expressed genes, gene ontology-enriched genes, and those detected via WGCNA.




