# Differential Gene Expression and Correlation Network Analysis of RNA-Seq Data From
# Upland Chorus Frog (*P. feriarum*) brains 

This repository contains scripts (R and Java) and input data asssociated to the manuscript
published in BMC Bioinformatics titled: Neurogenomic Divergence During Speciation by
Reinforcement of Mating Behaviors in Chorus Frogs (*Pseudacris*). The manuscript can be
accessed using DOI: XXXXX. Java code was contributed by Alan Lemmon and Python code by
Oscar Ospina. Code in other languages contributed by Alan Lemmon and Oscar Ospina.

The scripts and input files were used to detect differentially expressed genes in the
brain associated with synaptic trasnmission, as well as overall gene expression patterns.
An R script used to generate weighted correlation networks (WGCNA) is included. Code to
extract raw counts from BAM files is also provided.

## Pre-processing
The `rawReadCounter` folder contains scripts and an example file to obtain raw counts from 
BAM files (brain RNA-seq data mapped to *P. feriarum* reference transcriptome).
* `rawReadCounter/inputs/sample_I29357_BAMmappings.txt`: An example of input file (sample I29357,
allopatric female) from which raw counts were obtained. This file was directly extracted from
the headers of the corresponding BAM file.
* `rawReadCounter/code/ExtractRawCountsFromCovBed2*`: Counts reads mapped to a given transcript. Produces 
one file per sample. The counter is wrapped to loop through all samples via de `getRawCounts.sh`
script.
* `rawReadCounter/code/collatesumRawCols.py`: Takes output fles produced by the java raw counter and puts them
together in a .csv table (input for edgeR and most downstream analyses).

## Analysis/Statistical Assessment
The `IdentifyCandidateTranscripts` folder contains the script to find genes relevant to synaptic 
transmission during acoustic processing within the reference transcriptome generated in this study. 
* `IdentifyCandidateTranscripts/code/ExtractCandidatesFromAnnotation3.java`: This Java script 
uses a list of pre-selected terms (Additional File 11 of the supplemental materials) and the 
reference transcriptome annotations (DOI: 10.5281/zenodo.4709988). The list of genes related to 
synaptic transmission tested in this procedure was extracted from neurophysiological studies.

The `CandidateDEAnalysis` folder includes a script and input data to test for differential
expression of synaptic transmission genes as selected by `ExtractCandidatesFromAnnotation3.java`.
* `CandidateDEAnalysis/code/DE_genes_script_CandidateGenes.Rmd`: An RMarkdown document to run
the R script detecting differentially expressed genes via `edgeR`. The output includes
likelihood-ratio tests results, p-values, and FDR corrections.
* `CandidateDEAnalysis/inputs/CandidateTranscriptsFinal.xlsx`: The list of selected synaptic
transmission genes to test.
* `CandidateDEAnalysis/inputs/samples.csv`: Metadata of RNA-seq samples used in this study. The
same file is used for DE analysis on the whole transcriptome.
* `CandidateDEAnalysis/inputs/sumRawCounts.csv`: Raw counts obtained from BAM files of mapped
RNA-seq reads of *P. feriarum* brains. The same file is used for DE analysis on the whole 
transcriptome.

The `wholeTranscriptomeDEGenes` folder contains the code and inputs to test for differentially
expressed genes on the whole brain RNA-seq counts.
* `wholeTranscriptomeDEGenes/code/DE_genes_script.Rmd`: This RMarkdown document contains the
R code to test for differentially expressed genes in the whole brain RNA-seq transcriptome.
It uses the same raw counts and sample metadata of `CandidateDEAnalysis`.
* `wholeTranscriptomeDEGenes/code/heatmap_mtx.R`: An R function that outputs a matrix for
generation of heatmaps.
* `wholeTranscriptomeDEGenes/inputs/sprotAnnots_Trinotate_UniqueGenes_SYMBOLS.txt`: A file
containing Trinity contig IDs and their Uniprot symbols for annotation of DE genes.

The `wgcna_RNAseq` folder contains the script for generation of a weighted correlation 
network analysis (WGNCA).
* `wgcna_RNASeq/code/wgcna_script.R`: This R script generates a weighted correlation network
of the most variable genes among the samples. It takes as input the raw count file and sample
metadata in the `CandidateDEAnalysis` folder. The script also performs gene set enrichment
analysis via `g:Profiler`.

The `RandomizationTests` folder contains code and input data to generate randomization tests
for significant enrichment of synaptic activity-related genes within the WGCNA.
* `RandomizationTests/code/RandomizationTestsMethods.R`: R code to perform randomization tests on 
the WGCNA modules to test for enrichment in synaptic activity genes.
* `RandomizationTests/inputs/CPM_normalized_counts.csv.zip`: The counts-per-million normalized
matrix, resulting from processing of the raw counts of brain RNA-seq reads.
* `RandomizationTests/inputs/ModuleSummary.txt`: Summary data from the WGCNA.

