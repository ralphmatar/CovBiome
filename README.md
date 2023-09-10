# CovBiome
Pipeline for an RNA-seq exploratory analysis of the microbiome of COVID-19 patients collected from nasopharyngeal swabs

Welcome! This repository contains the scripts for my Master's thesis:

*Creation of a library of interfering sequences in SARS-Cov-2 nasal swabs, to improve molecular diagnostics and provide insight about the nasal microbiome*

In this thesis, I present a data description of 954 COVID-19 patients,
sequenced using Next Generation Sequencing (NGS) total RNA-seq from paired-end nasophayngeal swabs (NP).
This repositorz contains the scripts I used for the analysis of this data.

The main pipeline is a snakemake workflow *COMA.snakemake*, a second pipeline *performance.snakemake* is used
for benchmarking during the comparison of using 2 separate genomes or a concatenated one. The remaing python
and R scripts are for either generation of matrices or visualization.

I hope that this data and tools will be useful for researchers initiating in the bioinformatics analysis.
