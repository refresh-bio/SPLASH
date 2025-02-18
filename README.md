# SPLASH 2
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/splash/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/SPLASH/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/splash.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/splash)
<!-- [![Docker pulls](https://img.shields.io/github/downloads/refresh-bio/splash/total?label=downloads&logo=docker)](https://github.com/refresh-bio/SPLASH/pkgs/container/splash) --> <!-- DOES NOT WORK, shows jsut release download -->

## Introduction
SPLASH is an unsupervised and reference-free unifying framework to discover regulated sequence variation through statistical analysis of k-mer composition in both DNA and RNA sequence. 
It leverages our observation that detecting sample-regulated sequence variation, such as alternative splicing, RNA editing, gene fusions, V(D)J, transposable element mobilization, allele-specific splicing, genetic variation in a population, and many other regulated events can be unified in theory and in practice.
SPLASH analyzes the k-mer composition of raw sequencing reads to identify constant sequences (anchors) that are followed by sample-specific target variation and provides valid p-values (Chaung et al. 2023). 
SPLASH is reference-free, sidestepping the computational challenges associated with alignment, enabling fast discovery and statistical precision.

The first version of [SPLASH](https://github.com/salzman-lab/nomad/) implementedf in Python proved its usefulness. Here we provide SPLASH2, a new and improved implementation in C++ and Python (Kokot et al. 2024).
This new version is much more efficient and allows for the analysis of datasets >1TB size in hours on a workstation or even a laptop. 

We have also extended the SPLASH framework to barcoded single-cell and spatial analysis, called sc-SPLASH (Dehghannasiri et al. 2024), enabling the detection of regulated sequence variation at single-cell resolution in high-throughput single-cell (10x) and spatial (Visium) transcriptomics. sc-SPLASH is integrated into the SPLASH2 pipeline and can be invoked by setting the input parameter `technology = 10x` for 10x scRNA-Seq analysis or `technology = visium` for Visium spatial analysis.

## How does it work

A key concept of SPLASH is the analysis of composition of pairs of substrings *anchor*&ndash;*target* across many samples.
The substrings can be adjacent in reads or can be separated by a *gap*.

The image below presents the SPLASH pipeline on a high-level.
![image](https://github.com/refresh-bio/SPLASH/assets/9378882/8210fee0-c877-4374-9938-e3c01ea69e76)

<!-- ![image](https://user-images.githubusercontent.com/9378882/225988504-70266e4d-37e0-4c85-8c95-e47ad208cda9.png) -->

<!-- ![image](https://user-images.githubusercontent.com/9378882/224449978-309a4708-0fa1-4cb8-8483-a32e36ec2d58.png) -->

## Compactors

Compactors is a new statistical approach to local seed-based assembly. It comes as a part of SPLASH package and was particularly suited to assemble regions divere across across samples (see figure below). However, it can be used as an independent assembler on any types of seeds provided by the user.


![compactors-idea-v2](https://github.com/user-attachments/assets/49ee9aaf-54b8-4383-80be-3d225862e8bf)


## Installation, usage, example

Please visit our [Wiki page](../../wiki).

## References
Marek Kokot*, Roozbeh Dehghannasiri*, Tavor Baharav, Julia Salzman, and Sebastian Deorowicz.
[Scalable and unsupervised discovery from raw sequencing reads using SPLASH2](https://www.nature.com/articles/s41587-024-02381-2), Nature Biotechnology (2024)

Roozbeh Dehghannasiri*, Marek Kokot*, Sebastian Deorowicz, and Julia Salzman. [sc-SPLASH provides ultra-efficient reference-free discovery in barcoded single-cell sequencing](https://doi.org/10.1101/2024.12.24.630263), bioRxiv (2024)

Kaitlin Chaung*, Tavor Baharav*, George Henderson, Ivan Zheludev, Peter Wang, and Julia Salzman. [SPLASH: A statistical, reference-free genomic algorithm unifies biological discovery](https://www.cell.com/cell/fulltext/S0092-8674(23)01179-0), Cell (2023)
 
Tavor Baharav, David Tse, and Julia Salzman. 
[OASIS: An interpretable, finite-sample valid alternative to Pearson’s X2 for scientific discovery](https://www.pnas.org/doi/10.1073/pnas.2304671121), PNAS (2024)

George Henderson, Adam Gudys, Tavor Baharav, Punit Sundaramurthy, Marek Kokot, Peter L. Wang, Sebastian Deorowicz, Allison F. Carey, and Julia Salzman.
[Ultra-efficient, unified discovery from microbial sequencing with SPLASH and precise statistical assembly](https://www.biorxiv.org/content/10.1101/2024.01.18.576133v1.full)
bioRxiv 2024.01.18.576133 (2024)
