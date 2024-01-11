#!/usr/bin/env python3
import os
import sys
import subprocess
from time import gmtime, strftime
import threading
import queue
import tarfile
import shutil

##################################################
#       Config section
##################################################

#path to a directory where satc, satc_dump, satc_merge, sig_anch, kmc, kmc_toolsbinaries are
bin_path=os.path.join("../bin")

#number of threads for the first stage, too large value is not recomended because of intensive disk access here
n_threads_stage_1 = 4

#number of threads for the second stage, high value is recommended if possible, single thread will process single bin
n_threads_stage_2 = 32

#one may define a list of accepted anchors, this is plain text file with one anchor per line without any header
anchor_list = ""

#the data will be split in a number of bins that will be merged
#here other approach is used to have better bin size balance
n_bins = 128

anchor_len=27
target_len=27
#gap_len=23
gap_len="auto"

anchor_unique_targets_threshold=1
anchor_count_threshold=50
anchor_samples_threshold=1
anchor_sample_counts_threshold=5

#filter out all anchors containing poly(ACGT) of length at least <poly_ACGT_len> (0 means no filtering)
poly_ACGT_len=8
#number of most frequent tragets printed per each anchor in stats mode
n_most_freq_targets=2

#in calc_stats mode the number of altMaximize runs
generate_alt_max_cf_no_tires=10

#in calc_stats mode the number of iteration in altMaximize
altMaximize_iters=50

#in calc_stats mode use this fraction to create train X from contingency table
train_fraction=0.25
# True here may increase performance but also RAM-usage
kmc_use_RAM_only_mode=True

# if True statictics (pvals, etc.) will be computed
calculate_stats=True

#if True SVD will not be computed
without_SVD=True
#if True effect_size_cts will also be computed
with_effect_size_cts=False
# if True contignency tables will be generated in text format
dump_sample_anchor_target_count_txt=False

# if True contignency tables will be generated in binary (satc) format
# to convert to text format later satc_dump program may be used, it may take optionally maping from id to sample_name (--sample_names param)
dump_sample_anchor_target_count_binary=True

#correct pvals, calculate_stats must be set to True for this to work (if not the script will set it)
enable_pvals_correction=True
#for which column correction should be applied
pvals_correction_col_name="pval_rand_init_alt_max"
fdr_threshold=0.05

#dump Cjs
dump_Cjs=False

# dump only Cjs for anchors that have pval_rand_init_alt_max <=  max_pval_rand_init_alt_max_for_Cjs
max_pval_rand_init_alt_max_for_Cjs=0.1

# possible options here are
# * splash - text format like in nextflow splash
# * satc - different order of elements per line
# Warning: our approach uses just ints to represent sample, so it is unaware of a sample name
# in this script there is an unique identifier created for each sample, which may be used for decoding
# it is writen to "sample_name_to_id.mapping.txt" file
satc_merge_dump_format="splash"

#prefix of used to generate output file names
# each file name will start with this and it will have .bin{bin_id} suffix
outname_prefix="result"

#path to the file where input samples are defined
# the format is: per each line {sample_name}<space>{path}
# path is a fastq[.gz] file

input_file="input.txt"

#if True then intermediate files will be removed when not needed 
clean_up=True

# download the example data, after frist run should be set to false
download_data=False

#directory where logs will be stored
logs_dir="logs"

##################################################
#       End of config section
##################################################

if download_data:
    cmd = "(cd data; python ./download.py)"
    os.system(cmd)


def run_cmd(cmd):    
    #print(cmd)
    subprocess.run(cmd, shell=True)

_kmc_use_RAM_only_mode_param = "--kmc_use_RAM_only_mode" if kmc_use_RAM_only_mode else ""
_calculate_stats_param = "--calculate_stats" if calculate_stats else ""
_without_SVD_param = "--without_SVD" if without_SVD else ""
_with_effect_size_cts_param = "--with_effect_size_cts" if with_effect_size_cts else ""
_dump_sample_anchor_target_count_txt_param = "--dump_sample_anchor_target_count_txt" if dump_sample_anchor_target_count_txt else ""
_dump_sample_anchor_target_count_binary_param = "--dump_sample_anchor_target_count_binary" if dump_sample_anchor_target_count_binary else ""
_enable_pvals_correction_param = "--enable_pvals_correction" if enable_pvals_correction else ""
_clean_up_param = "--clean_up" if clean_up else ""
_dump_Cjs_param = "--dump_Cjs" if dump_Cjs else ""

_anchor_list_param = f"--anchor_list {anchor_list}" if anchor_list != "" else ""

cmd=f"./splash.py \
    --bin_path {bin_path} \
    --n_threads_stage_1 {n_threads_stage_1} \
    --n_threads_stage_2 {n_threads_stage_2} \
    {_anchor_list_param} \
    --n_bins {n_bins} \
    --anchor_len {anchor_len} \
    --target_len {target_len} \
    --gap_len {gap_len} \
    --poly_ACGT_len {poly_ACGT_len} \
    --anchor_unique_targets_threshold {anchor_unique_targets_threshold} \
    --anchor_count_threshold {anchor_count_threshold} \
    --anchor_samples_threshold {anchor_samples_threshold} \
    --anchor_sample_counts_threshold {anchor_sample_counts_threshold} \
    --n_most_freq_targets {n_most_freq_targets} \
    --generate_alt_max_cf_no_tires {generate_alt_max_cf_no_tires} \
    --altMaximize_iters {altMaximize_iters} \
    --train_fraction {train_fraction} \
    {_kmc_use_RAM_only_mode_param} \
    {_calculate_stats_param} \
    {_without_SVD_param} \
    {_with_effect_size_cts_param} \
    {_dump_sample_anchor_target_count_txt_param} \
    {_dump_sample_anchor_target_count_binary_param} \
    {_enable_pvals_correction_param} \
    --pvals_correction_col_name {pvals_correction_col_name} \
    --fdr_threshold {fdr_threshold} \
    --satc_merge_dump_format {satc_merge_dump_format} \
    --outname_prefix {outname_prefix} \
    {_dump_Cjs_param} \
    --max_pval_rand_init_alt_max_for_Cjs {max_pval_rand_init_alt_max_for_Cjs} \
    {_clean_up_param} \
    --logs_dir {logs_dir} \
    {input_file}"

run_cmd(cmd)
