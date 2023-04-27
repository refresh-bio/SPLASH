# NOMAD 2.0
## Introduction
NOMAD is an unsupervised and reference-free unifying framework to discover regulated sequence variation through statistical analysis of k-mer composition in both DNA and RNA sequence. 
NOMAD leverages our observation that detecting sample-regulated sequence variation, such as alternative splicing, RNA editing, gene fusions, V(D)J, transposable element mobilization, allele-specific splicing, genetic variation in a population, and many other regulated events can be unified–in theory and in practice.
This is achieved with a simple model, NOMAD, that analyzes k-mer composition of raw sequencing reads (Chaung et al. 2022). 
NOMAD finds constant sequences (anchors) that are followed by a set of sequences (targets) with sample-specific target variation and provides valid p-values. 
NOMAD is reference-free, sidestepping the computational challenges associated with alignment and making it significantly faster and more efficient than alignment, and enabling discovery and statistical precision not currently available, even from pseudo-alignment.

The first version of [NOMAD](https://github.com/salzman-lab/nomad/) pipeline proved its usefullness.
It was implemented mainly in Python with the use of NextFlow.
Here we provide a new and improved implementation based in C++ and Python (Kokot et al. 2023).
This new version is much more efficient and allows for the analysis of datasets >1TB size in hours on a workstation or even a laptop.

## How does it work

A key concept of NOMAD is the analysis of composition of pairs of substrings *anchor*&ndash;*target* across many samples.
The substrings can be adjacent in reads or can be separated by a *gap*.

The image below presents the NOMAD pipeline on a high-level.

![image](https://user-images.githubusercontent.com/9378882/225988504-70266e4d-37e0-4c85-8c95-e47ad208cda9.png)

<!-- ![image](https://user-images.githubusercontent.com/9378882/224449978-309a4708-0fa1-4cb8-8483-a32e36ec2d58.png) -->

## Installation
### Precompiled binaries
The easiest way to get nomad is to use [precompiled release](https://github.com/refresh-bio/NOMAD/releases).
To get the version 2.0.0 and run the example is is sufficient to do:
```
curl -L https://github.com/refresh-bio/NOMAD/releases/download/v2.0.0/nomad-2.0.0.linux.x64.tar.gz | tar xz
cd example
./run-example.sh
```
### Compile from sources
NOMAD is implemented as a number of applications written in the C++ programming language and a Python wrapper to run the whole pipeline.
Currently the software may be used only under Linux. 
A compiler supporting C++17 is needed to compile the code.
Use following snippet to install NOMAD.
```
git clone https://github.com/refresh-bio/nomad
cd nomad
make -j
sudo make install
```
#### Running the example
To verify the installation on small example one may perform:
```
cd example
./download.py #download examplary data
nomad input.txt #run the pipeline with default parameters
```
The result consists of two TSV files, namely,
 1. `result.after_correction.all_anchors.tsv` 
 2. `result.after_correction.scores.tsv`
The first file contain all unfiltered anchors found by the pipeline.
The second file contains only anchors whose corrected p-value is below 0.05.
## Understanding the output
There are following columns in the resulting tsv files
|                 Column            |     Meaning                            |     Notes                                                                              |
| ----------------------------------| -------------------------------------- | -------------------------------------------------------------------------------------- |
| anchor                            | anchor                                 |                                                                                        |
| pval_opt                          | p-value from alternating maximization  | number of iterations and other parameters can be configured with switches. Optimization formulation and statistical exposition in (Baharav et al. 2023).             |
| effect_size_bin                   | measure of anchor's effect size        | bounded between [0,1], indicates how well the data is divided between 2 groups of columns. 0 indicates no different between groups, 1 indicates two groups of columns with disjoint row supports (Baharav et al. 2023).                    |
| pval_asymp_base                   | asymptotically valid p-value based on alternating maximization                  | compute optimizing c and f used for pval_opt, and evaluate approximate p-value given by asymptotic normality (Baharav et al. 2023). |
| pval_base                         | p-value based on random partitionings      | Compute base p-value (using num_rand_cf random c,f), Bonferroni correction to yield valid p-value. |
| M                                 | anchor count                           | number of anchor occurences in the input (= the sum of elements in contingency table)  |
| anch_uniqTargs                    | number of unique targets for anchor    | (=number of rows in contingency table)                                                  |
| number_nonzero_samples            | number of samples containing anchor    | (=number of columns in contingency table)                                              |
| target_entropy                    | measure of diversity of target distribution         | entropy of empirical target distribution. Compute the empirical target distribution $(p_1,...,p_T)$ over the $T$ observed targets by summing counts across all samples and normalizing, and output the entropy of this empirical target distribution $H = -\sum_i p_i \log_2 (p_i)$.            |
| avg_hamming_distance_max_target   | average hamming distance to most abundant target | for each observed target (count), compute the hamming distance between it and the most abundant target (highest counts), average this over all targets. Useful in identifying / filtering SNPs (this measure will be low). |
| avg_hamming_distance_all_pairs    | average hamming distance between all pairs of targets | For the M counts in the matrix, compute the hamming distance between all pairs of targets, output the average.                                                                                       |
| avg_edit_distance_max_target      | average Levenshtein distance to most abundant target | Same as avg_hamming_distance_max_target, but Levenshtein distance as opposed to hamming distance. If an anchor has high avg_hamming_distance_max_target but low avg_edit_distance_max_target, this is indicative that the discrepancy between the top targets is due to an insertion or deletion. If the two quantities are similar, then this discrepancy is most likely due to a SNP.                                                                                        |
| avg_edit_distance_all_pairs       | average Levenshtein distance between all pairs of targets | Same as avg_hamming_distance_all_pairs but with Levenshtein distance as opposed to hamming distance.                                                                                       |
| pval_opt_corrected  | Benjamini-Yekutieli corrected pval_opt       |   only present in `*.scores.tsv` file                                                             |

## Input format
In the example the `input.txt` file was used. This file defines the set of input samples for the algorithm.
Its format is one sample per line. Each line should contain the name of a sample and (after space) path to the input sample.

**Important note:** if relative path is specified it is relative to the current working directory, not the directory of `input.txt`.

## Additional output
### Most frequent targets per each anchor
By default NOMAD will store 2 most frequent targets per each anchor in the resulting TSV files. This should be sufficient for splicing, but for RNA editing/missmatches 4 may be a better choice. It may be set with `--n_most_freq_targets` switch. If the number of targets for a given anchor is lower than specified value there will be a single `-` for each missing target.
### SATC format
NOMAD stores intermediate and optional output files in SATC format (**S**ample **A**nchor **T**arget **C**ount).
### Sample representation
The unique id is assigned to each sample. The ids are consecutive numbers starting with 0. The first sample from the input file gets id 0, the second one gets 1, and so on. By default NOMAD will store the mapping in `sample_name_to_id.mapping.txt` file, but this may be redefined with `--sample_name_to_id` parameter. This mapping file may be useful to access the data stored in SATC format.
### Output sample, anchor, target, count
#### Textual
When `--dump_sample_anchor_target_count_txt` switch is used there will be an additional output directory (named `result_dump` by default, but the `results` part may be redefined with `--outname_prefix` switch). This directory will contain a number of files (equal to the number of bins, default 128, may be redefined with `--n_bins` switch). The extension of these files is `.satc.dump`. Each line of these files is a tab-separated list of <sample_name> <anchor> <target> <count>. This is the easiest way to be able to reproduce contingency tables used during computation. Each file contains some number of anchors, but it is assured that a specific anchor is present in a single file. Since these text files could be large, it may be proficient to use binary (SATC) files instead.
#### Binary
When `--dump_sample_anchor_target_count_binary` switch is used there will be an additional output directory (named `result_satc` by default, but the `results` part may be redefined with `--outname_prefix` switch). This directory will contain a number of files (equal to the number of bins, default 128, may be redefined with `--n_bins` switch). The extension of these files is `.satc`. These are binary files in SATC format. Their content may be converted to textual representation with `satc_dump` program (part of the NOMAD package). Each file contains some number of anchors, but it is assured that a specific anchor is present in a single file. Since these text files could be large, it may be proficient to use binary (SATC) files instead, especially if one wants to investigate only some of all anchors.
### satc_dump
To convert SATC files into textual representation one may use `satc_dump` program. The simples usage is
```
 satc_dump input.satc output.satc.dump
```
There are also additional parameters that may be useful, namely:
 * `--anchor_list` &mdash; path to text file containing anchors separated by whitespaces, only anchors from this file will be dumped
 * `--sample_names` &mdash; path for decode sample id, each line should contain <sample_name> <sample_id>
 * `--n_bins <int>` &mdash; if set to value different than 0 the input is interpreted as a list of bins (each bin in separate line, first list is bin_0, second line is bin_1, etc. (in case of ill-formed input results will be incorrect)
 * `--separately` &mdash; if set with n_bins != 0 output param will be treated as suffix name and there will be output for each bin
 
 If `--sample_names` is not used in the output there will be sample ids instead of its names. NOMAD by default write mapping to `sample_name_to_id.mapping.txt` file (may be redefined with `--sample_name_to_id` switch of NOMAD.
 
 ### Interpreting results
 `plotGeneration.py` provides basic functionality to visualize contingency tables. To use, run NOMAD with the --dump_sample_anchor_target_count_binary flag. Then, `plotGeneration.py` can be run with inputs:
- `--dsName` &mdash; dataset name to be used in title of plots
- `outFolder` &mdash; path to save files to
- `metadataPath` &mdash; path to .tsv file containing two columns titled `sampleName` and `metadata`
- `satcFolder` &mdash; path to folder containing .satc files
- `pvDfPath` &mdash; path to p-value .tsv file (result.after_correction.scores.tsv)
- `sampleMappingTxt` &mdash; path to sample_name_to_id.mapping.txt file
- `anchorFile` &mdash; file of anchors to be plotted, one anchor per line
- `satc_dump_file` &mdash; path to satc_dump utility file, within NOMAD /bin folder
- `skipSATC` &mdash; optional flag, if .satc files have already been dumped for this set of anchors into the output folder, and the desire is to regenerate plots.

The generated output is 1 file per anchor, each file containing 4 subplots. At bottom left is the raw I x J contingency table, where each column is a sample, and each row represents a target (low abundance targets not shown). Each sample’s target counts are normalized by $n_j$ and plotted. At bottom right, the targets for the contingency table are displayed in a I x k table (corresponding to the rows of the contingency table, visually aligned). Each target of length k is shown, with basepairs colored as in the below colorbar. The target sequence for each row is displayed on the right. At top left,  the sample metadata is shown in a 1 x J table. Each entry corresponds to the column in the contingency table directly below it. The sample metadata is shown in a colorbar below. At top right, the column counts ($n_j$) are shown in a 1 x J table, where the colorbar below provides the scale. Again, columns are sorted as the contingency table.


The script `c_analysis.ipynb` shows how the saved c vectors can be loaded in for further analysis. `--dump_Cjs` must be enabled for this.
 
## Biological interpretation and classification of anchors
To facilitate downstream analysis of anchors, we provide a postprocessing script `NOMAD_extendor_classification.R`, that can be run on the anchors file generataed from the NOMAD run to classify anchors to biologically meaningful events such as alternative splicing, and base pair changes. `NOMAD_extendor_classification.R` needs the following inputs:

- `directory` &mdash; the output directory used for the NOMAD run
- `which_anchors_file` &mdash; flag to decide which anchor file (after correction or all anchors) to use, could be "after_correction" or "all" 
- `effect_size_cutoff` &mdash; the effect size cutoff for significant anchors (default 0.2) 
- `num_samples_cutoff` &mdash; the minimum number of sampels for an anchor to be called (default 20)
- `STAR_executable` &mdash; path to STAR executable file
- `STAR_reference` &mdash; path to STAR index files for the reference genome
- `samtools_executable` &mdash; path to samtools executable file
- `bedtools_executable` &mdash; path to bedtools executable file
- `annotated_splice_juncs` &mdash; path to the file containing annotated splice junctions from reference genome
- `annotated_exon_boundaries` &mdash; path to the file containing annotated exon boundaries from reference genome
- `bowtie2_executable` &mdash; path to bowtie2 executable file
- `bowtie2_univec_index` &mdash; path to bowtie2 index for univec
- `bowtie2_reference` &mdash; path to bowtie2 index for the reference genome
- `paralogs_file` &mdash; path to file containing list of paralogous genes from reference genome
 
The script will generate a file `classified_anchors.tsv` in the same directory used for NOMAD run, containing significant anchors along with their biological classification and alignment information.

## Annotation files for human genome (based on T2T):
The following three files are needed by the classification script for annotation anchors from human datasets and can be downloaded from the following links:
- `Annotated exon coordinates`: https://drive.google.com/file/d/1R-4-ICDAzmIBgQmlOF22nNrCWoSgrmHi/view?usp=share_link
- `Annotated splice junctions`: https://drive.google.com/file/d/1owlOQyP1z4cyFvYcAAA-qQmc-K6jGbs9/view?usp=share_link
- `Paralogous genes`: https://drive.google.com/file/d/1mqGft4tPlx8X3cRYoqQnPeXaonLfSbGa/view?usp=share_link

## Configuration
There is a lot of parameters allowing to customize the pipeline. They can be grouped into several categories. 
The parameters will be displayed when running nomad without parameters (or with `--help`).
 
 ### Base configuration
 * `--outname_prefix` &mdash; prefix of output file names (default: result)
 * `--anchor_len` &mdash; anchor length (default: 27)
 * `--gap_len` &mdash; gap length, if 'auto' it will be inferred from the data, in the opposite case it must be an int (default: 0)
 * `--target_len` &mdash; target length (default: 27)
 * `--anchor_list` &mdash; list of accepted anchors, this is path to plain text file with one anchor per line without any header (default accept all achnors)
 * `--pvals_correction_col_name` &mdash; for which column correction should be applied (default: pval_opt)
### Filters and thresholds:
 * `--poly_ACGT_len` &mdash; filter out all anchors containing poly(ACGT) of length at least <poly_ACGT_len> (0 means no filtering) (default: 8)
 * `--anchor_unique_targets_threshold` &mdash; filter out all anchors for which the number of unique targets is <= anchor_unique_targets_threshold (default: 1)
 * `--anchor_count_threshold` &mdash; filter out all anchors for which the total count <= anchor_count_threshold (default: 50)
 * `--anchor_samples_threshold` &mdash; filter out all anchors for which the number of unique samples is <= anchor_samples_threshold (default: 1)
 * `--anchor_sample_counts_threshold` &mdash; filter out anchor from sample if its count in this sample is <= anchor_sample_counts_threshold (default: 5)
 * `--n_most_freq_targets_for_stats` &mdash; use at most n_most_freq_targets_for_stats for each contignency table, 0 means use all (default: 0)
 * `--fdr_threshold` &mdash; keep anchors having corrected p-val below this value (default: 0.05)
### Additional output configuration:
 * `--dump_Cjs` &mdash; output Cjs (default: False)
 * `--max_pval_opt_for_Cjs` &mdash; dump only Cjs for anchors that have pval_opt <= max_pval_opt_for_Cjs (default: 0.1)
 * `--n_most_freq_targets` &mdash; number of most frequent tragets printed per each anchor in stats mode (default: 2)
 * `--with_effect_size_cts` &mdash; if set effect_size_cts will be computed (default: False)
 * `--sample_name_to_id` &mdash; file name with mapping sample name <-> sammpe id (default: sample_name_to_id.mapping.txt)
 * `--dump_sample_anchor_target_count_txt` &mdash; if set contignency tables will be generated in text format (default: False)
 * `--dump_sample_anchor_target_count_binary` &mdash; if set contignency tables will be generated in binary (SATC) format, to convert to text format later `satc_dump` program may be used, it may take optionally mapping from id to sample_name (--sample_names param) (default: False)
### Tuning statistics computation:
 * `--opt_num_inits` &mdash; the number of altMaximize random initializations (default: 10)
 * `--opt_num_iters` &mdash; the maximum number of iterations in altMaximize (default: 50)
 * `--num_rand_cf` &mdash; the number of random c and f used for pval_base (default: 50)
 * `--opt_train_fraction` &mdash; in calc_stats mode use this fraction to create training data from contingency table (default: 0.25)
### Technical and performance-related:
 * `--bin_path` &mdash; path to a directory where satc, satc_dump, satc_merge, sig_anch, kmc, kmc_tools binaries are (if any not found there nomad will check if installed and use installed) (default: ./)
 * `--tmp_dir` &mdash; path to a directory where temporary files will be stored (default: let nomad decide)
 * `--n_threads_stage_1` &mdash; number of threads for the first stage, too large value is not recomended because of intensive disk access here, but may be profitable if there is a lot of small size samples in the input (default: 4)
 * `--n_threads_stage_1_internal` &mdash; number of threads per each stage 1 thread (default: 8)
 * `--n_threads_stage_2` &mdash; number of threads for the second stage, high value is recommended if possible, single thread will process single bin (default: 32)
 * `--n_bins` &mdash; the data will be split in a number of bins that will be merged later (default: 128)
 * `--kmc_use_RAM_only_mode` &mdash; if set may increase performance but also RAM-usage (default: False)
 * `--kmc_max_mem_GB` &mdash; maximal amount of memory (in GB) KMC will try to not extend (default: 12)
 * `--dont_clean_up` &mdash; if set then intermediate files will not be removed (default: False)
 * `--logs_dir` &mdash; director where run logs of each thread will be stored (default: logs)

## References

Kaitlin Chaung, Tavor Baharav,  Ivan Zheludev, Julia Salzman. [A statistical, reference-free algorithm subsumes myriad problems in genome science and enables novel discovery](https://doi.org/10.1101/2022.06.24.497555), bioRxiv (2022)
 
Tavor Baharav, David Tse, and Julia Salzman. 
[An Interpretable, Finite Sample Valid Alternative to Pearson’s X2 for Scientific Discovery](https://www.biorxiv.org/content/10.1101/2023.03.16.533008), bioRxiv (2023)
 
Marek Kokot, Roozbeh Dehghannasiri, Tavor Baharav, Julia Salzman, and Sebastian Deorowicz.
[NOMAD2 provides ultra-efficient, scalable, and unsupervised discovery on raw sequencing reads] (https://..) 
bioRxiv (2023)
