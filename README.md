# NOMAD 2.0
## Introduction
NOMAD is an unsupervised and reference-free unifying framework to discover regulated sequence variation through statistical analysis of k-mer compulsion in both DNA and RNA sequence. 
NOMAD leverages our observation that detecting sample-regulated sequence variation, such as alternative splicing, RNA editing, gene fusions, V(D)J, transposable element mobilization, allele-specific splicing, and genetic variation in a population, among and many other regulated events can unified–in theory and in practice.
This is achieved with a simple model that analyzes k-mer composition of raw sequencing reads (Chaung et al. 2022). 
NOMAD finds constant sequences (anchors) that are followed by a set of sequences (targets) with sample-specific target variation and provides valid p-values (Chaung et al. 2022). 
NOMAD is reference-free, sidestepping the computational challenges associated with alignment and making it significantly faster and more efficient than alignment, and enabling discovery and statistical precision not currently available, even from pseudo-alignment.

The first version of [NOMAD](https://github.com/salzman-lab/nomad/) pipeline proved its usefullness.
It was implemented mainly in Python with a use of NextFlow.
Here we provide its new version implemented in the C++ and Pyton programming languages.

**TODO: ADD SOME MORE DESCRIBTION HERE, SOME BIOLOGICAL TEXT WOULD BE GREAT**

## How does it work

**TODO: HOW DETAILED WE WANT TO BE HERE??**

![image](https://user-images.githubusercontent.com/9378882/224449978-309a4708-0fa1-4cb8-8483-a32e36ec2d58.png)
## Installation
NOMAD is implemented as a number of applications written in the C++ programming language and a Python wrapper to run the whole pipeline.
Currently the software may be used only under Linux. 
A compiler supporting C++17 is needed to compile the code.
Use following snippet to install NOMAD.
```
git clone https://github.com/refresh-bio/nomad
cd nomad
make -j32
sudo make install
```
## Running the example
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
| pval_rand_init_alt_max            | p-value form alternating maximization  | number of iterations and other parameters can be configured with switches              |
| effect_size_bin                   | **TODO**                               |                                                                                        |
| pval_asymp_base                   | **TODO**                               |                                                                                        |
| pval_base                         | **TODO**                               |                                                                                        |
| M                                 | anchor count                           | number of anchor occurences in the input (= the sum of elements in contingency table)  |
| anch_uniqTargs                    | number of unique targets for anchor    | (=number of rows in contingency table)                                                 |
| number_nonzero_samples            | number of samples containing anchor    | (=number of columns in contingency table)                                              |
| target_entropy                    | **TODO**                               |                                                                                        |
| avg_hamming_distance_max_target   | **TODO**                               |                                                                                        |
| avg_hamming_distance_all_pairs    | **TODO**                               |                                                                                        |
| avg_edit_distance_max_target      | **TODO**                               |                                                                                        |
| avg_edit_distance_all_pairs       | **TODO**                               |                                                                                        |
| pval_rand_init_alt_max_corrected  | corrected pval_rand_init_alt_max       |   only `*.scores.tsv` file                                                             |

## Input format
In the example the `input.txt` file was used. This file defines the set of input samples for the algorithm.
Its format is one sample per line. Each line should contain the name of a sample and (after space) path to the input sample.

**Important note:** if relative path is specified it is relative to the current working directory, not the directory of `input.txt`.

## Configuration
There is a lot of parameters allowing to customize the pipeline. These parameters are splitted in a number of groups. 
These parameters will be displayed when running nomad without parameters (or with `--help`).
Below the groups with parameters are listed.
### Base configuration
 * `--outname_prefix` - prefix of output file names (default: result)
 * `--anchor_len` - anchor length (default: 27)
 * `--gap_len` - gap length, if 'auto' it will be inferred from the data, in the opposite case it must be an int (default: 0)
 * `--target_len` - target length (default: 27)
 * `--anchor_list` - list of accepted anchors, this is path to plain text file with one anchor per line without any header (default accept all achnors)
 * `--pvals_correction_col_name` - for which column correction should be applied (default: pval_rand_init_alt_max)
### Filters and thresholds:
 * `--poly_ACGT_len` - filter out all anchors containing poly(ACGT) of length at least <poly_ACGT_len> (0 means no filtering) (default: 8)
 * `--anchor_unique_targets_threshold` - filter out all anchors for which the number of unique targets is <= anchor_unique_targets_threshold (default: 1)
 * `--anchor_count_threshold` - filter out all anchors for which the total count <= anchor_count_threshold (default: 50)
 * `--anchor_samples_threshold` - filter out all anchors for which the number of unique samples is <= anchor_samples_threshold (default: 1)
 * `--anchor_sample_counts_threshold` - filter out anchor from sample if its count in this sample is <= anchor_sample_counts_threshold (default: 5)
 * `--n_most_freq_targets_for_stats` - use at most n_most_freq_targets_for_stats for each contignency table, 0 means use all (default: 0)
 * `--fdr_threshold` - keep anchors having corrected p-val below this value (default: 0.05)
### Additional output configuration:
 * `--dump_Cjs` - output Cjs (default: False)
 * `--max_pval_rand_init_alt_max_for_Cjs` - dump only Cjs for anchors that have pval_rand_init_alt_max <= max_pval_rand_init_alt_max_for_Cjs (default: 0.1)
 * `--n_most_freq_targets` - number of most frequent tragets printed per each anchor in stats mode (default: 0)
 * `--with_effect_size_cts` - if set effect_size_cts will be computed (default: False)
### Tuning statistics computation:
 * `--generate_alt_max_cf_no_tires` - the number of altMaximize runs (default: 10)
 * `--altMaximize_iters` - the number of iteration in altMaximize (default: 50)
 * `--num_rand_cf` - the number of rand cf (default: 50)
 * `--train_fraction` - in calc_stats mode use this fraction to create train X from contingency table (default: 0.25)
### Technical and performance-related:
 * `--bin_path` - path to a directory where satc, satc_dump, satc_merge, sig_anch, kmc, kmc_tools binaries are (if any not found there nomad will check if installed and use installed) (default: ./)
 * `--tmp_dir` - path to a directory where temporary files will be stored (default: let nomad decide)
 * `--n_threads_stage_1` - number of threads for the first stage, too large value is not recomended because of intensive disk access here, but may be profitable if there is a lot of small size samples in the input (default: 4)
 * `--n_threads_stage_1_internal` - number of threads per each stage 1 thread (default: 8)
 * `--n_threads_stage_2` - number of threads for the second stage, high value is recommended if possible, single thread will process single bin (default: 32)
 * `--n_bins` - the data will be split in a number of bins that will be merged later (default: 128)
 * `--kmc_use_RAM_only_mode` - if set may increase performance but also RAM-usage (default: False)
 * `--kmc_max_mem_GB` - maximal amount of memory (in GB) KMC will try to not extend (default: 12)
 * `--dont_clean_up` - if set then intermediate files will not be removed (default: False)
 * `--logs_dir` - director where run logs of each thread will be stored (default: logs)

## Sample representation
The unique id is assigned to each sample. The ids are consecutive numbers starting with 0. The first sample from the input file gets id 0, the second one gets 1, and so on.
