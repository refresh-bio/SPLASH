# SPLASH 2
## Introduction
SPLASH is an unsupervised and reference-free unifying framework to discover regulated sequence variation through statistical analysis of k-mer composition in both DNA and RNA sequence. 
SPLASH leverages our observation that detecting sample-regulated sequence variation, such as alternative splicing, RNA editing, gene fusions, V(D)J, transposable element mobilization, allele-specific splicing, genetic variation in a population, and many other regulated events can be unified–in theory and in practice.
This is achieved with a simple model, SPLASH, that analyzes k-mer composition of raw sequencing reads (Chaung et al. 2022). 
SPLASH finds constant sequences (anchors) that are followed by a set of sequences (targets) with sample-specific target variation and provides valid p-values. 
SPLASH is reference-free, sidestepping the computational challenges associated with alignment and making it significantly faster and more efficient than alignment, and enabling discovery and statistical precision not currently available, even from pseudo-alignment.

The first version of [SPLASH](https://github.com/salzman-lab/nomad/) pipeline proved its usefulness.
It was implemented mainly in Python with the use of NextFlow.
Here we provide a new and improved implementation based in C++ and Python (Kokot et al. 2023).
This new version is much more efficient and allows for the analysis of datasets >1TB size in hours on a workstation or even a laptop.

## How does it work

A key concept of SPLASH is the analysis of composition of pairs of substrings *anchor*&ndash;*target* across many samples.
The substrings can be adjacent in reads or can be separated by a *gap*.

The image below presents the SPLASH pipeline on a high-level.
![image](https://github.com/refresh-bio/SPLASH/assets/9378882/8210fee0-c877-4374-9938-e3c01ea69e76)

<!-- ![image](https://user-images.githubusercontent.com/9378882/225988504-70266e4d-37e0-4c85-8c95-e47ad208cda9.png) -->

<!-- ![image](https://user-images.githubusercontent.com/9378882/224449978-309a4708-0fa1-4cb8-8483-a32e36ec2d58.png) -->

## Installation
### Precompiled binaries
The easiest way to get SPLASH is to use [precompiled release](https://github.com/refresh-bio/SPLASH/releases).
To get version 2.3.0 and run the example, is is sufficient to do the following:
```
curl -L https://github.com/refresh-bio/SPLASH/releases/download/v2.3.0/splash-2.3.0.linux.x64.tar.gz | tar xz
cd example
./run-example.sh
```

### Docker container
It is possible to run splash using the docker container available through GitHub packages (https://github.com/orgs/refresh-bio/packages/container/package/splash).

To pull the image run:
```
sudo docker pull ghcr.io/refresh-bio/splash:2.3.0 # replace version number if needed
```
Example of how to run splash with docker.
Prerequisites:
 * There is `input.txt` in the current directory
 * All input splash files defined in the `input.txt` are in the same or child directory
```
sudo docker run -v `pwd`:/home/ubuntu ghcr.io/refresh-bio/splash:2.3.0 splash input.txt
```

### Singularity container
Sometimes, `sudo` is unavailable (for example, on HPC).
In such a case docker container may be transformed into a singularity container (https://docs.sylabs.io/guides/latest/user-guide/).

To pull the singularity version of splash use:
```
singularity pull docker://ghcr.io/refresh-bio/splash:2.3.0 # replace version as needed
```
This will result in a `splash_2.3.0.sif` file created in the current directory.
To execute splash using this file run:
```
./splash_2.3.0.sif splash input.txt
```
It is also possible to run without pulling first:
```
singularity run docker://ghcr.io/refresh-bio/splash:2.3.0 splash input.txt
```
It may be necessary to configure bind for singularity (`-B` parameter) depending on your configuration.

### Compile from sources
SPLASH is implemented as several applications written in the C++ programming language and a Python wrapper to run the whole pipeline.
A compiler supporting C++17 is needed to compile the code.
You can use the following snippet to compile SPLASH.
```
git clone --recurse-submodules https://github.com/refresh-bio/splash
cd splash
make -j
```
After this splash may be run as follows:
```
bin/splash # this will print help
```

#### Install the compiled splash
The simplest way of installing splash after compilation is to run:
```
sudo make install
```
If the `PREFIX` environment variable is not defined the above will install splash in `/usr/local/bin` which requires `sudo`.
It is possible to override the install location so that `sudo` is not required.
For example, to install in the user's home directory, one may run:
```
make install PREFIX=~/splash
```
or
```
export PREFIX=~/splash
make install
```

#### Uninstall splash
To uninstall splash one may run:
```
sudo make uninstall
```
The same `PREFIX` should be used like for installation (depending on this `sudo` may be not required).

## Running the example
To verify the installation on a small example, one may perform the following:
```
cd example
./download.py #download exemplary data
splash input.txt #run the pipeline with default parameters
```
The result consists of two TSV files, namely,
 1. `result.after_correction.all_anchors.tsv` 
 2. `result.after_correction.scores.tsv`
The first file contains all unfiltered anchors found by the pipeline.
The second file contains only anchors whose corrected p-value is below 0.05.

## Configuration
There are a lot of parameters allowing to customize the pipeline. They can be grouped into several categories. 
The parameters will be displayed when running splash without parameters (or with `--help`).
 ### Input
 * `input_file` — positional argument, path to the file where input samples are defined, the  format is: per each line 
 `{sample_name} {path}`
 path is a fastq[.gz] file in case of non-10X and txt
file for 10X/Visium where the content of text file is
`{first_file.fastq[.gz]},{second_file.fastq[.gz]}`
per line
 ### Base configuration
 * `--outname_prefix` — prefix of output file names (default: result)
 * `--anchor_len` — anchor length (default: 31)
 * `--gap_len` — gap length, if 'auto' it will be inferred from the data, in the opposite case it must be an int (default: 0)
 * `--target_len` — target length (default: 31)
 * `--anchor_list` — list of accepted anchors, this is path to plain text file with one anchor per line without any header (default accept all anchors)
 * `--pvals_correction_col_name` — for which column correction should be applied (default: pval_opt)
 * `--technology` — Technology used to generate the input data, must be one of 'base', '10x', 'visium' (default: `base`)
 * `--without_compactors` — if used compactors will not be run (default: False)
 * `--compactors_config` — optional json file with compactors configuration, example file content: { "num_threads": 4, "epsilon": 0.001 } (default: )
 * `--lookup_table_config` — optional json file with configuration of lookup_table, if not specified lookup_tables are not used  (default: )
### Filters and thresholds:
 * `--poly_ACGT_len` — filter out all anchors containing poly(ACGT) of length at least <poly_ACGT_len> (0 means no filtering) (default: 8)
 * `--artifacts` — path to artifacts, each anchor containing artifact will be filtered out (default: )
 * `--dont_filter_illumina_adapters` — if used anchors containing Illumina adapters will not be filtered out (default: False)
 * `--anchor_unique_targets_threshold` — filter out all anchors for which the number of unique targets is <= anchor_unique_targets_threshold (default: 1)
 * `--anchor_count_threshold` — filter out all anchors for which the total count <= anchor_count_threshold (default: 50)
 * `--anchor_samples_threshold` — filter out all anchors for which the number of unique samples is <= anchor_samples_threshold (default: 1)
 * `--anchor_sample_counts_threshold` — filter out anchor from sample if its count in this sample is <= anchor_sample_counts_threshold (default: 5)
 * `--n_most_freq_targets_for_stats` — use at most n_most_freq_targets_for_stats for each contingency table, 0 means use all (default: 0)
 * `--n_most_freq_targets_for_dump` — use when dumping satc (txt or binary), resulting file will only contain data for n_most_freq_targets_for_dump targets in each anchor, 0 means use all (default: 0)
 * `--fdr_threshold` — keep anchors having corrected p-val below this value (default: 0.05)
 * `--min_hamming_threshold` — keep only anchors with a pair of targets that differ by >= min_hamming_threshold (default: 0)
 * `--keep_top_n_target_entropy` — select keep_top_n_target_entropy records with highest target entropy (0 means don't select) (default: 10000)
 * `--keep_top_n_effect_size_bin` — select keep_top_n_effect_size_bin records with highest effect size bin (0 means don't select) (default: 20000)
 * `--keep_significant_anchors_satc` — if set there will be additional output file in SATC format with all significant anchors (default: False)
 * `--keep_top_target_entropy_anchors_satc` — if set there will be additional output file in SATC format with top target entropy significant anchors (default: False)
 * `--keep_top_effect_size_bin_anchors_satc` — if set there will be additional output file in SATC format with top effect size bin anchors (default: False)
### Additional output configuration:
 * `--dump_Cjs` — output Cjs (default: False)
 * `--max_pval_opt_for_Cjs` — dump only Cjs for anchors that have pval_opt <= max_pval_opt_for_Cjs (default: 0.1)
 * `--n_most_freq_targets` — number of most frequent targets printed per each anchor (default: 10)
 * `--with_effect_size_cts` — if set effect_size_cts will be computed (default: False)
 * `--with_pval_asymp_opt` — if set pval_asymp_opt will be computed (default: False)
 * `--without_seqence_entropy` — if set sequence entropy for anchor and most freq targets will not be computed (default: False)
 * `--sample_name_to_id` — file name with mapping sample name <-> sample id (default: sample_name_to_id.mapping.txt)
 * `--dump_sample_anchor_target_count_txt` — if set contingency tables will be generated in text format (default: False)
 * `--dump_sample_anchor_target_count_binary` — if set contingency tables will be generated in binary (SATC) format, to convert to text format later `satc_dump` program may be used, it may take optionally mapping from id to sample_name (--sample_names param) (default: False)
 * `--supervised_test_samplesheet` — if used script for finding/visualizing anchors with metadata-dependent variation will be run (forces `--dump_sample_anchor_target_count_binary`) (default: )
 * `--supervised_test_anchor_sample_fraction_cutoff` — the cutoff for the minimum fraction of samples for each anchor (default: 0.4)
 * `--supervised_test_num_anchors` — maximum number of anchors to be tested example (default: 20000)
### Tuning statistics computation:
 * `--opt_num_inits` — the number of altMaximize random initializations (default: 10)
 * `--opt_num_iters` — the maximum number of iterations in altMaximize (default: 50)
 * `--num_rand_cf` — the number of random c and f used for pval_base (default: 50)
 * `--num_splits` — the number of contingency table splits (default: 1)
 * `--opt_train_fraction` — use this fraction to create training data from contingency table (default: 0.25)
 * `--without_alt_max` — if set int alt max and related stats will not be computed (default: False)
 * `--without_sample_spectral_sum` — if set sample spectral sum will not be computed (default: False)
 * `--Cjs_samplesheet` — path to file with predefined Cjs for non-10X supervised mode (default: )
### Technical and performance-related:
 * `--bin_path` — path to a directory where satc, satc_dump, satc_merge, sig_anch, kmc, kmc_tools, and other binaries are (if any not found there splash will check if installed and use installed) (default: bin)
 * `--tmp_dir` — path to a directory where temporary files will be stored (default: let splash decide)
 * `--n_input_generators` — number of input generators running in parallel, used only when input generators are defined (0 means auto adjustment)  (default: 0)
 * `--n_threads_stage_1` — number of threads for the first stage, too large value is not recomended because of intensive disk access here, but may be profitable if there is a lot of small size samples in the input (0 means auto adjustment) (default: 0)
 * `--n_threads_stage_1_internal` — number of threads per each stage 1 thread (0 means auto adjustment) (default: 0)
 * `--n_threads_stage_1_internal_boost` — multiply the value of n_threads_stage_1_internal by this (may increase performance but the total number of running threads may be high) (default: 1)
 * `--n_threads_stage_2` — number of threads for the second stage, high value is recommended if possible, single thread will process single bin  (0 means auto adjustment) (default: 0)
 * `--n_bins` — the data will be split in a number of bins that will be merged later (default: 128)
 * `--kmc_use_RAM_only_mode` — True here may increase performance but also RAM-usage (default: False)
 * `--kmc_max_mem_GB` — maximal amount of memory (in GB) KMC will try to not extend (default: 12)
 * `--dont_clean_up` — if set then intermediate files will not be removed (default: False)
 * `--logs_dir` — director where run logs of each thread will be stored (default: logs)

### 10x/visium processing
 * `--cbc_len` — call barcode length (in case of 10X/Visium data) (default: 16)
 * `--umi_len` — UMI length (in case of 10X/Visium data) (default: 12)
 * `--soft_cbc_umi_len_limit` — allow additional symbols (beyond cbc_len + umi_len in _1.fastq 10X file UMI (default: 0)
 * `--cbc_filtering_thr` — how to filter cbcs, if 0 do the same as umi tools, in the opposite case keep cbcs with freq >= <cbc_filtering_thr> (default: 0)
 * `--cell_type_samplesheet` — path for mapping barcode to cell type, is used Helmert-based supervised mode is turned on (default: )
 * `--export_cbc_logs` — use if need cbc log files (default: False)
 * `--predefined_cbc` — path to file with predefined CBCs (default: )
 * `--export_filtered_input` — use if need filtered FASTQ files  (default: False)
 * `--allow_strange_cbc_umi_reads` — use to prevent the application from crashing when the CBC+UMI read length is outside the acceptable range (either shorter than CBC+UMI or longer than CBC+UMI+soft_cbc_umi_len_limit) (default: False)

### Postprocessing
 * `--postprocessing_item` — path to JSON defining postprocessing, may be defined multiple times, will be executed in the order of provided arguments (default: [])
 * `--exclude_postprocessing_item` — Path to JSON defining postprocessing to exclude from the default or provided postprocessing items (default: None)

## Understanding the output
There are following columns in the resulting tsv files
|                 Column            |     Meaning                            |     Notes                                                                              |
| ----------------------------------| -------------------------------------- | -------------------------------------------------------------------------------------- |
| anchor                            | anchor                                 |                                                                                        |
| pval_opt                          | p-value from alternating maximization  | number of iterations and other parameters can be configured with switches. Optimization formulation and statistical exposition in (Baharav et al. 2023).             |
| pval_sample_spectral_sum          | p-value from spectral sum              | may be disabled with `--without_sample_spectral_sum`                                   |
| effect_size_bin                   | measure of anchor's effect size        | bounded between [0,1], indicates how well the data is divided between 2 groups of columns. 0 indicates no different between groups, 1 indicates two groups of columns with disjoint row supports (Baharav et al. 2023).                    |
| pval_asymp_opt                    | asymptotically valid p-value based on alternating maximization                  | compute optimizing c and f used for pval_opt, and evaluate approximate p-value given by asymptotic normality (Baharav et al. 2023). Only if `--with_pval_asymp_opt` was used. |
| pval_base                         | p-value based on random partitionings      | Compute base p-value (using num_rand_cf random c,f), Bonferroni correction to yield valid p-value. |
| M                                 | anchor count                           | number of anchor occurences in the input (= the sum of elements in contingency table)  |
| anch_uniqTargs                    | number of unique targets for anchor    | (=number of rows in contingency table)                                                  |
| number_nonzero_samples            | number of samples containing anchor    | (=number of columns in contingency table)                                              |
| target_entropy                    | measure of diversity of target distribution         | entropy of empirical target distribution. Compute the empirical target distribution $(p_1,...,p_T)$ over the $T$ observed targets by summing counts across all samples and normalizing, and output the entropy of this empirical target distribution $H = -\sum_i p_i \log_2 (p_i)$.            |
| avg_no_homopolymer_targets        |                                        |                                                                                            |
| avg_hamming_distance_max_target   | average hamming distance to most abundant target | for each observed target (count), compute the hamming distance between it and the most abundant target (highest counts), average this over all targets. Useful in identifying / filtering SNPs (this measure will be low). |
| avg_hamming_distance_all_pairs    | average hamming distance between all pairs of targets | For the M counts in the matrix, compute the hamming distance between all pairs of targets, output the average.                                                                                       |
| avg_edit_distance_max_target      | average Levenshtein distance to most abundant target | Same as avg_hamming_distance_max_target, but Levenshtein distance as opposed to hamming distance. If an anchor has high avg_hamming_distance_max_target but low avg_edit_distance_max_target, this is indicative that the discrepancy between the top targets is due to an insertion or deletion. If the two quantities are similar, then this discrepancy is most likely due to a SNP.                                                                                        |
| avg_edit_distance_all_pairs       | average Levenshtein distance between all pairs of targets | Same as avg_hamming_distance_all_pairs but with Levenshtein distance as opposed to hamming distance.                                                                                       |
| anchor_2mer_seq_entropy           | 2-mer sequence entropy of the anchor                      |                                                                        |
| anchor_3mer_seq_entropy           | 3-mer sequence entropy of the anchor                      |                                                                        |
| most_freq_target_<i>              | `i`-th most frequent target for current anchor            | The number of reported targets may be controlled with `--n_most_freq_targets` |
| cnt_most_freq_target_<i>            | count of the `i`-th most frequent target for current anchor | The number of reported targets may be controlled with `--n_most_freq_targets` |
| most_freq_target_<i>_2mer_seq_entropy| 2-mer sequence entropy of the `i`-th most frequent target for current anchor | The number of reported targets may be controlled with `--n_most_freq_targets` |
| most_freq_target_<i>_3mer_seq_entropy| 3-mer sequence entropy of the `i`-th most frequent target for current anchor | The number of reported targets may be controlled with `--n_most_freq_targets` |
| pval_opt_corrected  | Benjamini-Yekutieli corrected pval_opt       |   only present in `*.scores.tsv` file                                                             |

## Input format
In the example, the `input.txt` file was used. This file defines the set of input samples for the algorithm.
Its format is one sample per line. Each line should contain the name of a sample and (after space) path to the input sample.

**Important note:** if a relative path is specified it is relative to the current working directory, not the directory of `input.txt`.

### Input generators
SPLASH handles natively fasta/fastq (gz or not) input files.
To support other formats (like SAM, BAM, CRAM, SRA) one may first prepare input data in the above format with an appropriate command.
It is also possible to use a built-in input generator mechanism.
Input format is:
```
<sample_name_i> <path_to_reads_i> [generator_i]
```
where the `generator_i` is an optional command that SPLASH will run before main processing.
**Warning:** This command will be run even if the file at `path_to_reads_i` already exists.
For example, to retrieve samples from SRA and run SPLASH one may define input as follows:
```
SRR30680504 SRR30680504_1.fastq fasterq-dump SRR30680504
SRR30680505 SRR30680505_1.fastq fasterq-dump SRR30680505
```
(of course, it may be profitable to configure fasterq-dump to generate fasta instead of fastq, or to use prefetch first, etc.).
The same rule may be applied to convert BAM/SAM/CRAM to fasta/q.
A couple of examples are in `example` directory:
 - `input-with-generators-wget.txt`
 - `input-with-generators-fasterq-dump.txt`

In the case of 10X/Visum data, the idea is similar, but instead of defining generating commands in the main input file, they are defined on a lower level (see an example in `example/input-10X-with-generators-wget.txt` and `example/S*-generators-wget.txt`).

The number of generators running in parallel may be defined with `--n_input_generators` parameter (if not set, this will be auto-adjusted).

## Example Applications
Given the broad applications of SPLASH, in analysis_notebooks folder, we provide the analysis downstream of SPLASH on how to interpret the results for a few major applications in which SPLASH has been applied so far. We should note that SPLASH is quite general and continues to be applied in new genomics problems.  
#### Splicing analysis for RNA-Seq
The notebook (`analysis_notebooks/SPLASH_splicing_analysis_notebook.Rmd`) provides detailed step-by-step instructions on how SPLASH results can be interpreted for an alternative splicing analysis. The reference genome in this example is human but it could be replaced with any other organism with any quality of transcriptome annotation. 

## Additional output
### Most frequent targets per each anchor
By default SPLASH will store 10 most frequent targets per each anchor in the resulting TSV files. 2 should be sufficient for splicing, but for RNA editing/mismatches 4 may be a better choice. It may be set with `--n_most_freq_targets` switch. If the number of targets for a given anchor is lower than specified value there will be a single `-` for each missing target.
### SATC format
SPLASH stores intermediate and optional output files in SATC format (**S**ample **A**nchor **T**arget **C**ount).
### Sample representation
The unique id is assigned to each sample. The ids are consecutive numbers starting with 0. The first sample from the input file gets id 0, the second one gets 1, and so on. By default SPLASH will store the mapping in `sample_name_to_id.mapping.txt` file, but this may be redefined with `--sample_name_to_id` parameter. This mapping file may be useful to access the data stored in SATC format.
### Output sample, anchor, target, count
#### Textual
When `--dump_sample_anchor_target_count_txt` switch is used there will be an additional output directory (named `result_dumps` by default, but the `results` part may be redefined with `--outname_prefix` switch).
This directory will contain a number of files (equal to the number of bins, default 128, may be redefined with `--n_bins` switch).
The extension of these files is `.satc.dump`.
Each line of these files is a tab-separated list of <sample_name> <anchor> <target> <count>.
This is the easiest way to be able to reproduce contingency tables used during computation.
Each file contains some number of anchors, but it is assured that a specific anchor is present in a single file.
Since these text files could be large, it may be proficient to use binary (SATC) files instead, especially if one wants to investigate only some of all anchors.
#### Binary
When `--dump_sample_anchor_target_count_binary` switch is used there will be an additional output directory (named `result_satc` by default, but the `results` part may be redefined with `--outname_prefix` switch).
This directory will contain a number of files (equal to the number of bins, default 128, may be redefined with `--n_bins` switch).
The extension of these files is `.satc`.
These are binary files in SATC format.
Their content may be converted to textual representation with `satc_dump` program (part of the SPLASH package).
Each file contains some number of anchors, but it is assured that a specific anchor is present in a single file.

### satc_dump
To convert SATC files into textual representation one may use `satc_dump` program.
The simples usage is
```
 satc_dump input.satc output.satc.dump
```
There are also additional parameters that may be useful, namely:
 * `--anchor_list` — path to text file containing anchors separated by whitespaces (or tsv file with header containing column named 'anchor'), only anchors from this file will be dumped
 * `--sample_names` — path for decode sample id, each line should contain <sample_name> <sample_id>
 * `--n_bins <int>` — if set to value different than 0 the input is interpreted as a list of bins (each bin in separate line, first list is bin_0, second line is bin_1, etc. (in case of ill-formed input results will be incorrect))
 * `--separately` — if set with n_bins != 0 output param will be treated as suffix name and there will be output for each bin
 * `--binary` — if set the output will be in binary SATC format instead of text (satc_dump my be used as a anchor-based filter)

 If `--sample_names` is not used in the output there will be sample ids instead of its names. SPLASH by default write mapping to `sample_name_to_id.mapping.txt` file (may be redefined with `--sample_name_to_id` switch of SPLASH.

 `satc_dump` may be also used to determine in which bin given anchor may be present (to which bin it is mapped basing on hash value).
 Below is the usage example for such a case.
 ```
 satc_dump --which-bin --n_bins <int> <anchor>
 ```
 
 ### Interpreting results
 `plotGeneration.py` provides basic functionality to visualize contingency tables. To use, run SPLASH with the --dump_sample_anchor_target_count_binary flag. Then, `plotGeneration.py` can be run with inputs:
- `--dsName` — dataset name to be used in title of plots
- `outFolder` — path to save files to
- `metadataPath` — path to .tsv file containing two columns titled `sampleName` and `metadata`
- `satcFolder` — path to folder containing .satc files
- `pvDfPath` — path to p-value .tsv file (result.after_correction.scores.tsv)
- `sampleMappingTxt` — path to sample_name_to_id.mapping.txt file
- `anchorFile` — file of anchors to be plotted, one anchor per line
- `satc_dump_file` — path to satc_dump utility file, within SPLASH /bin folder
- `skipSATC` — optional flag, if .satc files have already been dumped for this set of anchors into the output folder, and the desire is to regenerate plots.

The generated output is 1 file per anchor, each file containing 4 subplots. At bottom left is the raw I x J contingency table, where each column is a sample, and each row represents a target (low abundance targets not shown). Each sample’s target counts are normalized by $n_j$ and plotted. At bottom right, the targets for the contingency table are displayed in a I x k table (corresponding to the rows of the contingency table, visually aligned). Each target of length k is shown, with base pairs colored as in the below colorbar. The target sequence for each row is displayed on the right. At top left,  the sample metadata is shown in a 1 x J table. Each entry corresponds to the column in the contingency table directly below it. The sample metadata is shown in a color bar below. At top right, the column counts ($n_j$) are shown in a 1 x J table, where the colorbar below provides the scale. Again, columns are sorted as the contingency table.


The script `c_analysis.ipynb` shows how the saved c vectors can be loaded in for further analysis. `--dump_Cjs` must be enabled for this.

## Compactors

Compactors is a new statistical approach to local seed-based assembly coming as a part of the SPLASH package. While compactors is run automatically by SPLASH with anchors used as seeds, it can be also employed as a standalone local assembler. Here one can find a more detailed [documentation](./COMPACTORS.md).
 
## Biological interpretation and classification of anchors
To facilitate downstream analysis of anchors, we provide a postprocessing script `SPLASH_extendor_classification.R`, that can be run on the anchors file generated from the SPLASH run to classify anchors to biologically meaningful events such as alternative splicing, and base pair changes. `SPLASH_extendor_classification.R` needs the following inputs:

- `directory` — the output directory used for the SPLASH run
- `which_anchors_file` — flag to decide which anchor file (after correction or all anchors) to use, could be "after_correction" or "all" 
- `effect_size_cutoff` — the effect size cutoff for significant anchors (default 0.2) 
- `num_samples_cutoff` — the minimum number of samples for an anchor to be called (default 20)
- `STAR_executable` — path to STAR executable file
- `samtools_executable` — path to samtools executable file
- `bedtools_executable` — path to bedtools executable file
- `bowtie2_executable` — path to bowtie2 executable file
- `STAR_reference` — path to STAR index files for the reference genome
- `bowtie2_reference` — path to bowtie2 index for the reference genome
- `bowtie2_univec_index` — path to bowtie2 index for univec
- `annotated_splice_juncs` — path to the file containing annotated splice junctions from the reference transcriptome (can be either downloaded or generated from `SPLASH_build.R`)
- `annotated_exon_boundaries` — path to the file containing annotated exon boundaries from the reference transcriptome (can be either downloaded or generated from `SPLASH_build.R`)
- `gene_coords_file` — path to the file containing gene coordinates from the reference transcriptome (can be either downloaded or generated from `SPLASH_build.R`)
- `paralogs_file` — (optional) path to file containing list of paralogous genes from the reference genome
 
The script will generate a file `classified_anchors.tsv` in the same directory used for SPLASH run, containing significant anchors along with their biological classification and alignment information.

## Building index and annotation files needed for running classification script 
For running the classification script for a given reference genome/transcriptome you first need to obtain a fasta file for the reference genome and a gtf file for the transcriptome annotation. You then need to do the following 3 steps (note that all index/annotation files from these 3 steps should be generated from the same fasta and gtf file):
- **STAR index**: you need [STAR](https://github.com/alexdobin/STAR) index for reference genome. You can use default parameters to build the index: `STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index_files --genomeFastaFiles $fasta file$ --sjdbGTFfile $gtf file$`
- **Bowtie index**: you also need to [Bowtie2](https://github.com/BenLangmead/bowtie2) index that can be built using this command: `bowtie2-build $fasta file$ Bowtie_index_files/$Bowtie_index_name$`
- **Three annotation files**: three files are needed for annotated exon boundaries, annotated splice sites, and gene coordinates that should be built by running script `SPLASH_build.R`. You need three inputs for `SPLASH_build.R` script: `$gtf_file$` (absolute path to the gtf file), `$hisat2_directory$` (directory containing HISAT2 codes downloaded from [HISAT2 repository](https://github.com/DaehwanKimLab/hisat2), the script assumes that there are two python scripts at: `$hisat2_directory$/extract_exons.py` and `$hisat2_directory$/extract_splice_sites.py`), `$outfile_name$` (the name used for the annotation files that script will generate). The script can be run using the following command:  
`Rscript SPLASH_build.R $gtf_file$ $hisat2_directory$ $outfile_name$`
If the script is run successfully, it will generate 3 output annotation files in the same directory as the script: `$outfile_name$_known_splice_sites.txt` (for annotated splice sites), `$outfile_name$_exon_coordinates.bed` (for annotated exon boundaries), and `$outfile_name$_genes.bed` (for annotated gene coordinates)

## Downloading pre-built annotation files for human and mouse genomes:
 The human files were built for both [T2T assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) and [GRCh38 assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). The mouse files were built based on [mm39 assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/). The annotation files can be downloaded using the following links:
- **Human (T2T)**:
   - `Annotated exon coordinates`: https://drive.google.com/file/d/1R-4-ICDAzmIBgQmlOF22nNrCWoSgrmHi/view?usp=share_link
   - `Annotated splice junctions`: https://drive.google.com/file/d/1owlOQyP1z4cyFvYcAAA-qQmc-K6jGbs9/view?usp=share_link
   - `Gene coordinates`: https://drive.google.com/file/d/1L0A7iGXEYiOsPQ0QiJayKPybJ79ZDi2F/view?usp=sharing
   - `Paralogous genes`: https://drive.google.com/file/d/1mqGft4tPlx8X3cRYoqQnPeXaonLfSbGa/view?usp=share_link
- **Human (GRCh38)**:
   - `Annotated exon coordinates`: https://drive.google.com/file/d/1oK6OgQnFFVvybBo0EZ5aIyeoZLAtMyZF/view?usp=sharing
   - `Annotated splice junctions`: https://drive.google.com/file/d/1izVHy1m-ddlNgJtFKfWcHdtkc_Y5bHHP/view?usp=sharing
   - `Gene coordinates`: https://drive.google.com/file/d/1REfnl9ZNYcsb-1jSurDHcsL7QFJ00JEp/view?usp=sharing
 - **Mouse (mm39)**:
   - `Annotated exon coordinates`: https://drive.google.com/file/d/1npE0rkxhsDtJk3FeMdfuZwc5Elfuk4bq/view?usp=sharing
   - `Annotated splice junctions`: https://drive.google.com/file/d/1iJhf421nMRDC0uCo_0jh7Nkns8NAieTE/view?usp=sharing
   - `Gene coordinates`: https://drive.google.com/file/d/1V8By-yq7AmgXY-XDhipgjjsamL0ghhJa/view?usp=sharing
   - `Paralogous genes`: https://drive.google.com/file/d/1Uf28bE2XF9Y2w57ARlUGfO5agiFFXx2S/view?usp=sharing

## References
Marek Kokot, Roozbeh Dehghannasiri, Tavor Baharav, Julia Salzman, and Sebastian Deorowicz.
[Scalable and unsupervised discovery from raw sequencing reads using SPLASH2](https://www.nature.com/articles/s41587-024-02381-2), Nature Biotechnology (2024)
 
Kaitlin Chaung, Tavor Baharav,  Ivan Zheludev, Julia Salzman. [A statistical, reference-free algorithm subsumes myriad problems in genome science and enables novel discovery](https://doi.org/10.1101/2022.06.24.497555), bioRxiv (2022)
 
Tavor Baharav, David Tse, and Julia Salzman. 
[An Interpretable, Finite Sample Valid Alternative to Pearson’s X2 for Scientific Discovery](https://www.biorxiv.org/content/10.1101/2023.03.16.533008), bioRxiv (2023)
