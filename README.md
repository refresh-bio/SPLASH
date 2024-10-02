# SPLASH 2
## Introduction
SPLASH is an unsupervised and reference-free unifying framework to discover regulated sequence variation through statistical analysis of k-mer composition in both DNA and RNA sequence. 
SPLASH leverages our observation that detecting sample-regulated sequence variation, such as alternative splicing, RNA editing, gene fusions, V(D)J, transposable element mobilization, allele-specific splicing, genetic variation in a population, and many other regulated events can be unified–in theory and in practice.
This is achieved with a simple model, SPLASH, that analyzes k-mer composition of raw sequencing reads (Chaung et al. 2022). 
SPLASH finds constant sequences (anchors) that are followed by a set of sequences (targets) with sample-specific target variation and provides valid p-values. 
SPLASH is reference-free, sidestepping the computational challenges associated with alignment and making it significantly faster and more efficient than alignment, and enabling discovery and statistical precision not currently available, even from pseudo-alignment.

The first version of [SPLASH](https://github.com/salzman-lab/nomad/) pipeline proved its usefulness.
It was implemented mainly in Python with the use of NextFlow.
Here we provide a new and improved implementation based in C++ and Python (Kokot et al. 2024).
This new version is much more efficient and allows for the analysis of datasets >1TB size in hours on a workstation or even a laptop.

## How does it work

A key concept of SPLASH is the analysis of composition of pairs of substrings *anchor*&ndash;*target* across many samples.
The substrings can be adjacent in reads or can be separated by a *gap*.

The image below presents the SPLASH pipeline on a high-level.
![image](https://github.com/refresh-bio/SPLASH/assets/9378882/8210fee0-c877-4374-9938-e3c01ea69e76)

<!-- ![image](https://user-images.githubusercontent.com/9378882/225988504-70266e4d-37e0-4c85-8c95-e47ad208cda9.png) -->

<!-- ![image](https://user-images.githubusercontent.com/9378882/224449978-309a4708-0fa1-4cb8-8483-a32e36ec2d58.png) -->

## Compactors

Compactors is a new statistical approach to local seed-based assembly. It comes as a part of SPLASH package and was particularly suited to assemble regions divere across across samples (see figure below). However, it can be used as an independent assembler on any types of seeds provided by the user.

![compactors-idea-v2](https://github.com/user-attachments/assets/0d9e4fac-949f-452c-baff-0dd842979899)


## Installation, usage, example

Please visit our [Wiki page](https://github.com/refresh-bio/SPLASH/wiki).

## References
Marek Kokot, Roozbeh Dehghannasiri, Tavor Baharav, Julia Salzman, and Sebastian Deorowicz.
[Scalable and unsupervised discovery from raw sequencing reads using SPLASH2](https://www.nature.com/articles/s41587-024-02381-2), Nature Biotechnology (2024), https://doi.org/10.1038/s41587-024-02381-2
 
Kaitlin Chaung, Tavor Baharav,  Ivan Zheludev, Julia Salzman. [A statistical, reference-free algorithm subsumes myriad problems in genome science and enables novel discovery](https://doi.org/10.1101/2022.06.24.497555), bioRxiv (2022)
 
Tavor Baharav, David Tse, and Julia Salzman. 
[An Interpretable, Finite Sample Valid Alternative to Pearson’s X2 for Scientific Discovery](https://www.biorxiv.org/content/10.1101/2023.03.16.533008), bioRxiv (2023)


George Henderson, Adam Gudys, Tavor Baharav, Punit Sundaramurthy, Marek Kokot, Peter L. Wang, Sebastian Deorowicz, Allison F. Carey, Julia Salzman.
[Ultra-efficient, unified discovery from microbial sequencing with SPLASH and precise statistical assembly](https://www.biorxiv.org/content/10.1101/2024.01.18.576133v1.full)
bioRxiv 2024.01.18.576133 (2024)
