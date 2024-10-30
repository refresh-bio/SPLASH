#!/usr/bin/env python3
import argparse
import os
import sys
import re
from time import gmtime, strftime, localtime
import subprocess
import threading
import multiprocessing
import queue
import gzip
import shutil
import uuid
from sys import platform
import time
import json

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

SPLASH_VERSION="2.11.2"

parser = argparse.ArgumentParser(
                    prog = "splash",
                    description = "Welcome to SPLASH\nVersion: " + SPLASH_VERSION,
                    #epilog = 'Text at the bottom of help',
                    #formatter_class=SmartFormatter
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                    )

parser.add_argument("input_file", help="path to the file where input samples are defined, the format is: per each line {sample_name}<space>{path}, path is a fastq[.gz] file in case of non-10X and txt file for 10X/Visium where the content of text file is {first_file.fastq[.gz]},{second_file.fastq[.gz]} per line")

def check_anchor_len(len):
    ilen = int(len)
    if ilen < 1 or ilen > 31:
        raise argparse.ArgumentTypeError(f"anchor_len must be in [1;31]")
    return ilen

def check_target_len(len):
    ilen = int(len)
    if ilen < 1 or ilen > 31:
        raise argparse.ArgumentTypeError(f"target_len must be in [1;31]")
    return ilen

EXPERIMENTAL_INPUT_GENERATORS = False

group_base_configuration = parser.add_argument_group('Base configuration')
group_base_configuration.add_argument("--outname_prefix", default="result", type=str, help="prefix of output file names")
group_base_configuration.add_argument("--anchor_len", default=31, type=check_anchor_len, help="anchor length")
group_base_configuration.add_argument("--gap_len", default="0", type=str, help="gap length, if 'auto' it will be inferred from the data, in the opposite case it must be an int")
group_base_configuration.add_argument("--target_len", default=31, type=check_target_len,  help="target length")
group_base_configuration.add_argument("--anchor_list", default="", type=str, help="list of accepted anchors, this is path to plain text file with one anchor per line without any header")
group_base_configuration.add_argument("--pvals_correction_col_name", default="pval_opt", type=str, help="for which column correction should be applied")
group_base_configuration.add_argument("--technology", default="base", type=str, help="Technology used to generate the input data, must be one of 'base', '10x', 'visium'")
group_base_configuration.add_argument("--without_compactors", default=False, action='store_true', help="if used compactors will not be run")
group_base_configuration.add_argument("--compactors_config", default="", type=str, help="optional json file with compactors configuration, example file content: { \"num_threads\": 4, \"epsilon\": 0.001 }")
group_base_configuration.add_argument("--lookup_table_config", default="", type=str, help="optional json file with configuration of lookup_table, if not specified lookup_tables are not used")

group_filters_and_thresholds = parser.add_argument_group('Filters and thresholds')
group_filters_and_thresholds.add_argument("--poly_ACGT_len", default=8, type=int, help="filter out all anchors containing poly(ACGT) of length at least <poly_ACGT_len> (0 means no filtering)")
group_filters_and_thresholds.add_argument("--artifacts", default="", type=str, help="path to artifacts, each anchor containing artifact will be filtered out")
group_filters_and_thresholds.add_argument("--dont_filter_illumina_adapters", default=False, action='store_true', help="if used anchors containing Illumina adapters will not be filtered out")
group_filters_and_thresholds.add_argument("--anchor_unique_targets_threshold", default=1, type=int, help="filter out all anchors for which the number of unique targets is <= anchor_unique_targets_threshold")
group_filters_and_thresholds.add_argument("--anchor_count_threshold", default=50, type=int, help="filter out all anchors for which the total count <= anchor_count_threshold")
group_filters_and_thresholds.add_argument("--anchor_samples_threshold", default=1, type=int, help="filter out all anchors for which the number of unique samples is <= anchor_samples_threshold")
group_filters_and_thresholds.add_argument("--anchor_sample_counts_threshold", default=5, type=int, help="filter out anchor from sample if its count in this sample is <= anchor_sample_counts_threshold")
group_filters_and_thresholds.add_argument("--n_most_freq_targets_for_stats", default=0, type=int, help="use at most n_most_freq_targets_for_stats for each contingency table, 0 means use all")
group_filters_and_thresholds.add_argument("--n_most_freq_targets_for_dump", default=0, type=int, help="use when dumping satc (txt or binary), resulting file will only contain data for n_most_freq_targets_for_dump targets in each anchor, 0 means use all")
group_filters_and_thresholds.add_argument("--fdr_threshold", default=0.05, type=float, help="keep anchors having corrected p-val below this value")
group_filters_and_thresholds.add_argument("--min_hamming_threshold", default=0, type=int, help="keep only anchors with a pair of targets that differ by >= min_hamming_threshold")
group_filters_and_thresholds.add_argument("--keep_top_n_target_entropy", default=10000, type=int, help="select keep_top_n_target_entropy records with highest target entropy (0 means don't select)")
group_filters_and_thresholds.add_argument("--keep_top_n_effect_size_bin", default=20000, type=int, help="select keep_top_n_effect_size_bin records with highest effect size bin (0 means don't select)")
group_filters_and_thresholds.add_argument("--keep_significant_anchors_satc", default=False, action='store_true', help="if set there will be additional output file in SATC format with all significant anchors")
group_filters_and_thresholds.add_argument("--keep_top_target_entropy_anchors_satc", default=False, action='store_true', help="if set there will be additional output file in SATC format with top target entropy significant anchors")
group_filters_and_thresholds.add_argument("--keep_top_effect_size_bin_anchors_satc", default=False, action='store_true', help="if set there will be additional output file in SATC format with top effect size bin anchors")

group_additional_out = parser.add_argument_group('Additional output configuration')
group_additional_out.add_argument("--dump_Cjs", default=False, action='store_true', help="output Cjs")
group_additional_out.add_argument("--max_pval_opt_for_Cjs", default=0.10, type=float, help="dump only Cjs for anchors that have pval_opt <= max_pval_opt_for_Cjs")
group_additional_out.add_argument("--n_most_freq_targets", default=10, type=int, help="number of most frequent targets printed per each anchor")
group_additional_out.add_argument("--with_effect_size_cts", default=False, action='store_true', help="if set effect_size_cts will be computed")
group_additional_out.add_argument("--with_pval_asymp_opt", default=False, action='store_true', help="if set pval_asymp_opt will be computed")
group_additional_out.add_argument("--without_seqence_entropy", default=False, action='store_true', help="if set sequence entropy for anchor and most freq targets will not be computed")
group_additional_out.add_argument("--sample_name_to_id", default="sample_name_to_id.mapping.txt", type=str, help="file name with mapping sample name <-> sample id")
group_additional_out.add_argument("--dump_sample_anchor_target_count_txt", default=False, action='store_true', help="if set contingency tables will be generated in text format")
group_additional_out.add_argument("--dump_sample_anchor_target_count_binary", default=False, action='store_true', help="if set contingency tables will be generated in binary (satc) format, to convert to text format later satc_dump program may be used, it may take optionally mapping from id to sample_name (--sample_names param)")
group_additional_out.add_argument("--satc_merge_dump_format", default="splash", type=str, choices=['splash', 'satc'], help="splash - text format like in nextflow, satc - different order of elements per line")
group_additional_out.add_argument("--supervised_test_samplesheet", default="", help="if used script for finding/visualizing anchors with metadata-dependent variation will be run (forces --dump_sample_anchor_target_count_binary)")
group_additional_out.add_argument("--supervised_test_anchor_sample_fraction_cutoff", default=0.4, type=float, help="the cutoff for the minimum fraction of samples for each anchor")
group_additional_out.add_argument("--supervised_test_num_anchors", default=20000, type=int, help="maximum number of anchors to be tested example")

group_tuning_stats = parser.add_argument_group('Tuning statistics computation')
group_tuning_stats.add_argument("--opt_num_inits", default=10, type=int, help="the number of altMaximize random initializations")
group_tuning_stats.add_argument("--opt_num_iters", default=50, type=int, help="the number of iteration in altMaximize")
group_tuning_stats.add_argument("--num_rand_cf", default=50, type=int, help="the number of random c and f used for pval_base")
group_tuning_stats.add_argument("--num_splits", default=1, type=int, help="the number of contingency table splits")
group_tuning_stats.add_argument("--opt_train_fraction", default=0.25, type=float, help="use this fraction to create train X from contingency table")
group_tuning_stats.add_argument("--without_alt_max", default=False, action='store_true', help="if set int alt max and related stats will not be computed")
group_tuning_stats.add_argument("--without_sample_spectral_sum", default=False, action='store_true', help="if set sample spectral sum will not be computed")
group_tuning_stats.add_argument("--compute_also_old_base_pvals", default=False, action='store_true', help="if set old pvals and effect size will be computed")
group_tuning_stats.add_argument("--Cjs_samplesheet", default="", help="path to file with predefined Cjs for non-10X supervised mode")

group_technical = parser.add_argument_group('Technical and performance-related')

group_technical.add_argument("--bin_path", default="bin", type=str, help="path to a directory where satc, satc_dump, satc_merge, sig_anch, kmc, kmc_tools, bkc, dsv_manip, gap_shortener, compactors binaries are")
group_technical.add_argument("--tmp_dir", default="", type=str, help="path to a directory where temporary files will be stored")

if EXPERIMENTAL_INPUT_GENERATORS:
    group_technical.add_argument("--n_input_generators", default=0, type=int, help="number of input generators running in parallel, used only when input generators are defined (0 means auto adjustment)")

group_technical.add_argument("--n_threads_stage_1", default=0, type=int, help="number of threads for the first stage, too large value is not recomended because of intensive disk access here, but may be profitable if there is a lot of small size samples in the input (0 means auto adjustment)")
group_technical.add_argument("--n_threads_stage_1_internal", default=0, type=int, help="number of threads per each stage 1 thread (0 means auto adjustment)")
group_technical.add_argument("--n_threads_stage_1_internal_boost", default=1, type=int, help="multiply the value of n_threads_stage_1_internal by this (may increase performance but the total number of running threads may be high)")
group_technical.add_argument("--n_threads_stage_2", default=0, type=int, help="number of threads for the second stage, high value is recommended if possible, single thread will process single bin (0 means auto adjustment)")
group_technical.add_argument("--n_bins", default=128, type=int, help="the data will be split in a number of bins that will be merged later")
group_technical.add_argument("--kmc_use_RAM_only_mode", default=False, action='store_true', help="True here may increase performance but also RAM-usage")
group_technical.add_argument("--kmc_max_mem_GB", default=12, type=int, help="maximal amount of memory (in GB) KMC will try to not extend")
group_technical.add_argument("--dont_clean_up", default=False, action='store_true', help="if set then intermediate files will not be removed")
group_technical.add_argument("--logs_dir", default="logs", type=str, help="director where run logs of each thread will be stored")

group_10X_Visium = parser.add_argument_group('10x/visium processing')
group_10X_Visium.add_argument("--cbc_len", default=16, type=int, help="call barcode length (in case of 10X/Visium data)")
group_10X_Visium.add_argument("--umi_len", default=12, type=int, help="UMI length (in case of 10X/Visium data)")
group_10X_Visium.add_argument("--soft_cbc_umi_len_limit", default=0, type=int, help="allow additional symbols (beyond cbc_len + umi_len in _1.fastq 10X file UMI")
group_10X_Visium.add_argument("--cbc_filtering_thr", default=0, type=int, help="how to filter cbcs, if 0 do the same as umi tools, in the opposite case keep cbcs with freq >= <cbc_filtering_thr>")
group_10X_Visium.add_argument("--cell_type_samplesheet", default="", help="path for mapping barcode to cell type, is used Helmert-based supervised mode is turned on")
group_10X_Visium.add_argument("--export_cbc_logs", default=False, action='store_true', help="use if need cbc log files")
group_10X_Visium.add_argument("--predefined_cbc", default="", help="path to file with predefined CBCs")
group_10X_Visium.add_argument("--export_filtered_input", default=False, action='store_true', help="use if need filtered FASTQ files")
group_10X_Visium.add_argument("--allow_strange_cbc_umi_reads", default=False, action='store_true', help="use to prevent the application from crashing when the CBC+UMI read length is outside the acceptable range (either shorter than CBC+UMI or longer than CBC+UMI+soft_cbc_umi_len_limit)")

group_postprocessing = parser.add_argument_group("Postprocessing")
#parser.add_argument("--postprocessing_item", action="append", default=["postprocessing/blast.json"], type=str, help="path to JSON defining postprocessing, may be defined multiple times, will be executed in the order of provided arguments")
group_postprocessing.add_argument("--postprocessing_item", action="append", default=[], type=str, help="path to JSON defining postprocessing, may be defined multiple times, will be executed in the order of provided arguments")
group_postprocessing.add_argument("--exclude_postprocessing_item", action="append", type=str, help="Path to JSON defining postprocessing to exclude from the default or provided postprocessing items")


#Hidden options
group_10X_Visium.add_argument("--REMOVE_INPUT", default=False, action='store_true', help=argparse.SUPPRESS)

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

bin_path=args.bin_path
if EXPERIMENTAL_INPUT_GENERATORS:
    n_input_generators = args.n_input_generators

n_threads_stage_1 = args.n_threads_stage_1
n_threads_stage_1_internal = args.n_threads_stage_1_internal
n_threads_stage_1_internal_boost = args.n_threads_stage_1_internal_boost
n_threads_stage_2 = args.n_threads_stage_2
anchor_list = args.anchor_list
n_bins = args.n_bins
cell_type_samplesheet=args.cell_type_samplesheet
Cjs_samplesheet=args.Cjs_samplesheet
cbc_len = args.cbc_len
umi_len = args.umi_len
soft_cbc_umi_len_limit = args.soft_cbc_umi_len_limit
allow_strange_cbc_umi_reads = args.allow_strange_cbc_umi_reads
anchor_len = args.anchor_len
target_len = args.target_len
gap_len = args.gap_len
poly_ACGT_len=args.poly_ACGT_len
artifacts = args.artifacts
dont_filter_illumina_adapters = args.dont_filter_illumina_adapters
anchor_unique_targets_threshold = args.anchor_unique_targets_threshold
anchor_count_threshold = args.anchor_count_threshold
anchor_samples_threshold = args.anchor_samples_threshold
anchor_sample_counts_threshold = args.anchor_sample_counts_threshold
cbc_filtering_thr = args.cbc_filtering_thr
export_cbc_logs = args.export_cbc_logs
predefined_cbc = args.predefined_cbc
export_filtered_input = args.export_filtered_input
n_most_freq_targets = args.n_most_freq_targets
n_most_freq_targets_for_stats = args.n_most_freq_targets_for_stats
n_most_freq_targets_for_dump = args.n_most_freq_targets_for_dump
opt_num_inits=args.opt_num_inits
opt_num_iters=args.opt_num_iters
num_rand_cf=args.num_rand_cf
num_splits=args.num_splits
opt_train_fraction=args.opt_train_fraction
kmc_use_RAM_only_mode = args.kmc_use_RAM_only_mode
kmc_max_mem_GB = args.kmc_max_mem_GB
without_alt_max = args.without_alt_max
without_sample_spectral_sum = args.without_sample_spectral_sum
with_effect_size_cts = args.with_effect_size_cts
with_pval_asymp_opt = args.with_pval_asymp_opt
without_seqence_entropy = args.without_seqence_entropy
sample_name_to_id = args.sample_name_to_id
compute_also_old_base_pvals = args.compute_also_old_base_pvals
dump_sample_anchor_target_count_txt = args.dump_sample_anchor_target_count_txt
dump_sample_anchor_target_count_binary = args.dump_sample_anchor_target_count_binary
pvals_correction_col_name = args.pvals_correction_col_name
technology = args.technology
without_compactors = args.without_compactors
compactors_config = args.compactors_config
lookup_table_config = args.lookup_table_config
fdr_threshold = args.fdr_threshold
min_hamming_threshold = args.min_hamming_threshold
keep_top_n_target_entropy = args.keep_top_n_target_entropy
keep_top_n_effect_size_bin = args.keep_top_n_effect_size_bin
keep_significant_anchors_satc = args.keep_significant_anchors_satc
keep_top_target_entropy_anchors_satc = args.keep_top_target_entropy_anchors_satc
keep_top_effect_size_bin_anchors_satc = args.keep_top_effect_size_bin_anchors_satc
satc_merge_dump_format = args.satc_merge_dump_format
supervised_test_samplesheet = args.supervised_test_samplesheet
supervised_test_anchor_sample_fraction_cutoff = args.supervised_test_anchor_sample_fraction_cutoff
supervised_test_num_anchors = args.supervised_test_num_anchors
outname_prefix=args.outname_prefix
dump_Cjs=args.dump_Cjs
max_pval_opt_for_Cjs=args.max_pval_opt_for_Cjs
clean_up=not args.dont_clean_up
logs_dir=args.logs_dir
input_file=args.input_file
tmp_dir=args.tmp_dir
postprocessing_item = args.postprocessing_item
exclude_postprocessing_item = args.exclude_postprocessing_item

REMOVE_INPUT = args.REMOVE_INPUT

def list_remove_duplicates_keep_order(in_list):
    S = set()
    out_list = []
    for item in in_list:
        if item not in S:
            out_list.append(item)
            S.add(item)
    return out_list

postprocessing_item = list_remove_duplicates_keep_order(postprocessing_item) #make postprocessings unique
if exclude_postprocessing_item:
    postprocessing_item = [item for item in postprocessing_item if item not in args.exclude_postprocessing_item]

#determine postprocessings actual paths (may be where splash main script is or in the pwd, if in both pwd overrides)
postprocessing_item_fixed_paths = []
for pp in postprocessing_item:
    script_path = os.path.dirname(__file__)
    if os.path.exists(pp):
        postprocessing_item_fixed_paths.append(pp)
    elif os.path.exists(os.path.join(script_path, pp)):
        postprocessing_item_fixed_paths.append(os.path.join(script_path, pp))
    else:
        print(f"Error: cannot find postprocessing definition file: {pp}")
        sys.exit(1)
postprocessing_item = postprocessing_item_fixed_paths



def get_cur_time():
    return strftime("%Y-%m-%d %H:%M:%S", localtime())

print("Welcome to SPLASH")
print("Version: ", SPLASH_VERSION)
print("Current time:", get_cur_time(), flush=True)

print("----------------  Configuration  ---------------")
for arg, value in vars(args).items():
    print("%s: %s" % (arg, value))
print("------------------------------------------------", flush=True)


if technology not in ("base", "10x", "visium"):
    print("Error: --technology must be one of 'base', '10x', 'visium'")
    sys.exit(1)

def is_10x_or_visium():
    return technology == "10x" or technology == "visium"

if supervised_test_samplesheet != "":
    if dump_sample_anchor_target_count_binary == False:
        print("Warning: --supervised_test_samplesheet forces --dump_sample_anchor_target_count_binary")
        dump_sample_anchor_target_count_binary = True
    if shutil.which("Rscript") is None:
        print("Error: Rscript must be installed to run splash with --supervised_test_samplesheet")
        sys.exit(1)

    if n_most_freq_targets < 2:
        print("Warning: --supervised_test_samplesheet requires at least 2 n_most_freq_targets, setting --n_most_freq_targets 2")
        n_most_freq_targets = 2

if tmp_dir == "":
    tmp_dir="splash-tmp-"+uuid.uuid4().hex

# check for duplicates in sample names
with open(input_file) as f:
    sample_names = set()
    for line in f:
        line = line.strip();
        if line == "": continue
        sample_name = line.split()[0]
        if sample_name in sample_names:
            print(f"Error: duplicated sample name {sample_name} in file {input_file}")
            sys.exit(1)
        sample_names.add(sample_name)

os.makedirs(tmp_dir)

#check if generators are used to make input files
#if so lets use them to create file and create new input file
#without generators for further run

if EXPERIMENTAL_INPUT_GENERATORS:
    input_generators_commands = []
    path_to_input_without_generators = f"{tmp_dir}/input_no_generators.txt"
    paths_to_input_without_generators_for_10x = []
    with open(path_to_input_without_generators, "w") as of:
        with open(input_file) as f:
            line_no = -1
            for line in f:
                line = line.strip()
                if line == "": continue

                line_no += 1

                if is_10x_or_visium():
                    _10x_sample_name, _10x_path = line.split()
                    path_to_input_without_generators_for_10x = f"{tmp_dir}/input_no_generators_10x_{line_no}.txt"
                    paths_to_input_without_generators_for_10x.append(path_to_input_without_generators_for_10x)
                    of.write(f"{_10x_sample_name} {path_to_input_without_generators_for_10x}\n")
                    with open(path_to_input_without_generators_for_10x, "w") as of_10x:
                        with open(_10x_path) as _10x_file:
                            for line2 in _10x_file:
                                line2 = line2.strip()
                                if line2 == "": continue
                                match = re.match(f"(\S+)\s*(.*)?", line2)

                                if not match:
                                    print(f"Error: wrong line in file: {_10x_path}")
                                    sys.exit(1)
                                before_command = match.group(1)
                                if match.group(2):
                                    input_generators_commands.append(match.group(2))

                                of_10x.write(f"{before_command}\n")
                else:
                    match = re.match(r"(\S+\s+\S+)\s*(.*)?", line)
                    if not match:
                        print(f"Error: wrong line in input file: {line}")
                        sys.exit(1)

                    before_command = match.group(1)
                    if match.group(2):
                        input_generators_commands.append(match.group(2))

                    of.write(f"{before_command}\n")

was_error = False

class StopAllBecauseCrititcalError(Exception):
    pass

def check_and_handle_error():
    global was_error
    if was_error:
        if threading.current_thread() is threading.main_thread():
            print("Exiting because of previous error")
            sys.exit(1)
        raise StopAllBecauseCrititcalError()

def wrap_cmd_with_time(cmd):
    if platform == "darwin":
        if shutil.which("gtime") is not None:
            cmd = f"gtime -v {cmd}"
        else:
            cmd = f"/usr/bin/time {cmd}" # no -v for /usr/bin/time on mac os
    else:
        cmd = f"/usr/bin/time -v {cmd}"
    return cmd

def run_cmd(cmd, out, err):
    global was_error
    check_and_handle_error()
    cmd = wrap_cmd_with_time(cmd)
    out.write(get_cur_time() + ": " + cmd + "\n")
    out.flush()
    p = subprocess.Popen(cmd, stderr=err,stdout=out, shell=True)
    p.communicate()
    if p.returncode != 0:
        print(f"Error running command: {cmd}", flush=True)
        if out.name == err.name:
            print(f"For details check {out.name}", flush=True)
        else:
            print(f"For details check {out.name} or {err.name}", flush=True)
        was_error = True
        check_and_handle_error()


class action_at_function_exit:
    def __init__(self, action):
        self.action = action
    def __del__(self):
        self.action()

if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

max_cpus_to_use_in_auto_adjust = min(multiprocessing.cpu_count(), 64)

if EXPERIMENTAL_INPUT_GENERATORS:
    if n_input_generators == 0:
        n_input_generators = max_cpus_to_use_in_auto_adjust

if EXPERIMENTAL_INPUT_GENERATORS:
    if len(input_generators_commands) > 0:
        print("Generating input files by commands")
        print("Current time:", get_cur_time(), flush=True)

        input_generators_queue = queue.Queue()

        def input_generators_worker(tid):
            tid = str(tid)
            while len(tid) < 4:
                tid = "0" + tid

            stdoutfile = open(f"{logs_dir}/input_generator_thread-{tid}.log", "w")
            while True:
                item = input_generators_queue.get()
                if item is None:
                    break
                command = item[0]

                action_at_function_exit(lambda: input_generators_queue.task_done()) # will be called even in a case of exception

                try:
                    run_cmd(command, stdoutfile, stdoutfile)
                except StopAllBecauseCrititcalError: #just consume next task
                    pass

            stdoutfile.close()

        input_generator_threads = []
        for i in range(n_input_generators):
            t  =threading.Thread(target=input_generators_worker, args=(i,))
            t.start()
            input_generator_threads.append(t)

        for command in input_generators_commands:
            input_generators_queue.put((command, ))

        input_generators_queue.join()

        # stop workers
        for i in range(n_input_generators):
            input_generators_queue.put(None)
        for t in input_generator_threads:
            t.join()
        check_and_handle_error()

        print("Generating input files by commands done")

    input_file = path_to_input_without_generators


if n_threads_stage_1 == 0 and n_threads_stage_1_internal == 0:
    n_input_samples = 0
    max_input_sample_size = 0 # maybe I should consider plan vs gz in fasta/fastq variants, but lets keep it simple
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if line == "": continue

            n_input_samples += 1
            path = line.split()[1]
            size = 0
            if is_10x_or_visium():
                with open(path) as f2:
                    for line2 in f2:
                        line2 = line2.strip()
                        if line2 == "": continue
                        paths = line2.split(',')
                        size += os.path.getsize(paths[0]) + os.path.getsize(paths[1])
            else:
                size = os.path.getsize(path)

            if size > max_input_sample_size:
                max_input_sample_size = size

    n_threads_stage_1 = min(max_cpus_to_use_in_auto_adjust, 4, n_input_samples)

    #if there is a lot of samples and all are quite small prefer more n_threads_stage_1
    if n_threads_stage_1 < max_cpus_to_use_in_auto_adjust // 2 and n_input_samples > 64 and max_input_sample_size < 2000000000:
        n_threads_stage_1 = max_cpus_to_use_in_auto_adjust // 2

    n_threads_stage_1_internal = max_cpus_to_use_in_auto_adjust // n_threads_stage_1

    print(f"n_threads_stage_1 auto adjusted to {n_threads_stage_1}", flush=True)
    print(f"n_threads_stage_1_internal auto adjusted to {n_threads_stage_1_internal}", flush=True)
elif n_threads_stage_1_internal == 0:
    n_threads_stage_1_internal = max(1, max_cpus_to_use_in_auto_adjust // n_threads_stage_1)
    print(f"n_threads_stage_1_internal auto adjusted to {n_threads_stage_1_internal}", flush=True)
elif n_threads_stage_1 == 0:
    n_threads_stage_1 = max(1, max_cpus_to_use_in_auto_adjust // n_threads_stage_1_internal)
    print(f"n_threads_stage_1 auto adjusted to {n_threads_stage_1}", flush=True)

if n_threads_stage_1_internal_boost != 1:
    n_threads_stage_1_internal *= n_threads_stage_1_internal_boost
    print(f"n_threads_stage_1_internal after applying boost: {n_threads_stage_1_internal}")

if n_threads_stage_2 == 0:
    n_threads_stage_2 = max_cpus_to_use_in_auto_adjust
    print(f"n_threads_stage_2 auto adjusted to {n_threads_stage_2}", flush=True)

def get_gap_len_from_the_data():
    def get_first_line_from_txt_file(path):
        try:
            with open(path) as f:
                return f.readline().strip()
        except:
            print(f"Error: cannot read first line from {path}")
            sys.exit(1)

    path = get_first_line_from_txt_file(input_file).split()[1]

    if is_10x_or_visium():
        path = get_first_line_from_txt_file(path).split(',')[1]

    try:
        with gzip.open(path,'rt') if path.endswith(".gz") else open(path) as f:
            f.readline() #skip header
            seq_len = len(f.readline().strip())
    except:
        print(f"Error: cannot read first line from {path}")
        sys.exit(1)

    return max(0, (seq_len - (target_len + anchor_len))//2)


if gap_len == "auto":
    gap_len = get_gap_len_from_the_data()
    print("gap_len inferred from input data:", gap_len, flush=True)
else:
    if not gap_len.isdigit():
        print("Error: --gap_len must be 'auto' or integer value")
        sys.exit(1)
    gap_len = int(gap_len)



if opt_train_fraction <= 0.0 or opt_train_fraction >= 1.0:
    print("Error: opt_train_fraction must be in range (0.0;1.0)", flush=True)
    sys.exit(1)

if is_10x_or_visium():
    if cbc_len > 16:
        print("Error: max cbc_len is 16", flush=True)
        sys.exit(1)
    if umi_len > 16:
        print("Error: max umi_len is 16", flush=True)
        sys.exit(1)

def get_prog_or_fail(prog):
    p = os.path.join(bin_path, prog)
    if os.path.exists(p):
        return p

    #maybe it is in the same dir where cur script
    script_path = os.path.dirname(__file__)
    p = os.path.join(script_path, prog)
    if os.path.exists(p):
        return p

    #maybe it is installed
    if shutil.which(prog) is not None:
        return prog

    print(f"Error: {prog} does not exist in {bin_path} nor it is installed.")
    sys.exit(1)

bkc=get_prog_or_fail("bkc")
dsv_manip=get_prog_or_fail("dsv_manip")
satc=get_prog_or_fail("satc")
satc_dump=get_prog_or_fail("satc_dump")
satc_merge=get_prog_or_fail("satc_merge")
sig_anch=get_prog_or_fail("sig_anch")
kmc=get_prog_or_fail("kmc")
kmc_tools=get_prog_or_fail("kmc_tools")
gap_shortener=get_prog_or_fail("gap_shortener")
read_selector=get_prog_or_fail("read_selector")
compactors=get_prog_or_fail("compactors")
lookup_table=get_prog_or_fail("lookup_table")
tsv_to_fasta=get_prog_or_fail("tsv_to_fasta")

if supervised_test_samplesheet != "":
    supervised_test = get_prog_or_fail("supervised_test.R")

#mkokot_TODO: do I enforce compactors if needed?
class SplashPostprocessing:
    def __run_postprocessing_cmd(cmd, out, err):
        #cmd = wrap_cmd_with_time(cmd)
        p = subprocess.Popen(cmd, stderr=err,stdout=out, shell=True)
        p.communicate()
        if p.returncode != 0:
            print(f"Warning: error ocurred running command: {cmd}", flush=True)
            if out.name == err.name:
                print(f"For details check {out.name}", flush=True)
            else:
                print(f"For details check {out.name} or {err.name}", flush=True)

    def __fill_placeholders(command, replacements):
        placeholders = re.findall(r"\{([a-zA-Z0-9_]+)\}", command)
        def replace(match):
            placeholder = match.group(1)
            return replacements.get(placeholder, f"{{{placeholder}}}")  # Keep the original placeholder if not in replacements

        updated_command = re.sub(r"\{([a-zA-Z0-9_]+)\}", replace, command)
        return updated_command


    def __load_postprocessing_item(self, path):
        required_fields = ["name", "command"]
        optional_fields = ["test_run"]

        #field cannot be required and optional at the same time
        assert len(set(optional_fields).intersection(set(required_fields))) == 0
        known_fields = required_fields + optional_fields
        try:
            with open(path) as f:
                res = json.load(f)

                for x in res.keys():
                    if not x in known_fields:
                        print(f"Error: unknown field '{x}', allowed fields:", ", ".join(known_fields))
                        sys.exit(1)
                for x in required_fields:
                    if not x in res:
                        print(f"Error: {path} does not contain required field \"{x}\"")
                        sys.exit(1)

                command = res['command']

                res['path'] = path
                name = res['name']
                invalid_chars = '<>:"/\\|?* '
                if any(c in set(invalid_chars) for c in name):
                    print(f"Error: field name ({name}) is invalid (contains one of '{invalid_chars}')")
                    sys.exit(1)

                if "test_run" in res:
                    try:
                        cmd_res = subprocess.run(res["test_run"], shell=True, check=True, text=True, capture_output=True)
                    except subprocess.CalledProcessError as e:
                        print(f"Error: test_run for {name} failed with code {e.returncode}")
                        print(f"Output:\n", e.stderr)
                        sys.exit()

        except BaseException as e:
            if isinstance(e, SystemExit):
                raise

            print(f"Error: cannot read postprocessing config from {path}, \n{e}")
            sys.exit(1)
        return res

    """postprocessing_items is a list of strings, each representing path to postprocessing definition file, items on this list should be unique"""
    def __init__(self, postprocessing_items : list) -> None:
        self.postprocessing_list = []
        postprocessing_names = {}
        for x in postprocessing_items:
            self.postprocessing_list.append(self.__load_postprocessing_item(x))
            name = self.postprocessing_list[-1]["name"]
            if name in postprocessing_names.keys():
                postprocessing_names[name] += 1
            else:
                postprocessing_names[name] = 1

        for name, cnt in postprocessing_names.items():
            if cnt > 1:
                print(f"Warning: postprocessing with name {name} was defined more {cnt} times, results of all will be in the same output directory. It is recommended to use unique names for postprocessings")

    def __assure_fasta(input_tsv, outpupt_fasta, log_file_path):
        #if fasta was already created do nothing
        if os.path.exists(outpupt_fasta):
            return

        cmd = f"{tsv_to_fasta} \
            {input_tsv} \
            {outpupt_fasta}"

        with open(log_file_path, "w") as stdoutfile:
            run_cmd(cmd, stdoutfile, stdoutfile)

    def RequireTopTargetEntropyExtendors(self):
        for pp in self.postprocessing_list:
            for x in pp["run_on"]:
                if x["type"] == "top_target_entropy_extendors":
                    return True
        return False

    def RequireTopEffectSizeExtendors(self):
        for pp in self.postprocessing_list:
            for x in pp["run_on"]:
                if x["type"] == "top_effect_size_extendors":
                    return True
        return False

    def Run(self):
        if len(self.postprocessing_list) > 0:
            print("Start splash postprocessing")

        paths = {
            "all_extendors_tsv": f"{outname_prefix}.after_correction.all_anchors.tsv",
            "significant_extendors_tsv": f"{outname_prefix}.after_correction.scores.tsv",
            "top_target_entropy_extendors_tsv": f"{outname_prefix}.after_correction.scores.top_target_entropy.tsv",
            "top_effect_size_extendors_tsv": f"{outname_prefix}.after_correction.scores.top_effect_size_bin.tsv",
            "compactors_top_target_entropy_tsv": f"{compactors_dir}/after_correction.scores.top_target_entropy.tsv",
            "compactors_top_effect_size_bin_tsv": f"{compactors_dir}/after_correction.scores.top_effect_size_bin.tsv",

            "all_extendors_fasta": f"{outname_prefix}.after_correction.all_anchors.extendors.fasta",
            "significant_extendors_fasta": f"{outname_prefix}.after_correction.scores.extendors.fasta",
            "top_target_entropy_extendors_fasta": f"{outname_prefix}.after_correction.scores.top_target_entropy.extendors.fasta",
            "top_effect_size_extendors_fasta": f"{outname_prefix}.after_correction.scores.top_effect_size_bin.extendors.fasta",
            "compactors_top_target_entropy_fasta": f"{compactors_dir}.after_correction.scores.top_target_entropy.fasta",
            "compactors_top_effect_size_bin_fasta": f"{compactors_dir}/after_correction.scores.top_effect_size_bin.fasta"
        }

        conver_to_fasta_log_paths = {
            "all_extendors_fasta": f"{logs_dir}/fasta_from.all_anchors.log",
            "significant_extendors_fasta": f"{logs_dir}/fasta_from.scores.log",
            "top_target_entropy_extendors_fasta": f"{logs_dir}/fasta_from.scores.top_target_entropy.log",
            "top_effect_size_extendors_fasta": f"{logs_dir}/fasta_from.scores.top_effect_size_bin.log",
            "compactors_top_target_entropy_fasta": f"{logs_dir}/fasta_from.compactors.scores.top_target_entropy.log",
            "compactors_top_effect_size_bin_fasta": f"{logs_dir}/fasta_from.compactors.scores.top_effect_size_bin.log"
        }

        #INFO: if user uses non absolut path (so most cases) I prepend ../ because it will be run from inside postprocessing directory
        if not os.path.isabs(outname_prefix):
            assert not os.path.isabs(compactors_dir)
            for key in paths:
                paths[key] = os.path.join("..", paths[key])

        if not os.path.isabs(logs_dir):
            for key in conver_to_fasta_log_paths:
                conver_to_fasta_log_paths[key] = os.path.join("..", conver_to_fasta_log_paths[key])



        allowed_placeholders = paths.keys()

        for pp in self.postprocessing_list:
            name = pp["name"]
            path = pp["path"]
            command = pp["command"]

            placeholders_in_command = re.findall(r'{(.*?)}', command)

            used_placeholders = [ph for ph in placeholders_in_command if ph in allowed_placeholders]

            unknown_placeholders = set(placeholders_in_command) - set(used_placeholders)
            for uph in unknown_placeholders:
                print(f"Warning: in command \"{command}\" unknown placeholder ({uph}) was used.")

            pp_dir = f"{outname_prefix}_{name}"
            if not os.path.exists(pp_dir):
                os.makedirs(pp_dir)

            pp_stdout_file = open(f"{pp_dir}/stdout.log", "w")
            pp_stderr_file = open(f"{pp_dir}/stderr.log", "w")
            pp_command_file = open(f"{pp_dir}/command.log", "w")

            #lets run postprocessing in its output directory
            os.chdir(pp_dir)

            for ph in used_placeholders:
                if ph.endswith("_fasta"):
                    prefix = ph[0:-len("_fasta")]

                    input_tsv = [p for p in paths if p.startswith(prefix) and p.endswith("_tsv")]
                    assert len(input_tsv) == 1
                    input_tsv = input_tsv[0]
                    SplashPostprocessing.__assure_fasta(paths[input_tsv], paths[ph], conver_to_fasta_log_paths[ph])

            

            print(f"Start postprocessing: {name}: {path}")
            
            parsed_command = SplashPostprocessing.__fill_placeholders(command, paths)
            pp_command_file.write(parsed_command + '\n');
            pp_command_file.flush()
            SplashPostprocessing.__run_postprocessing_cmd(parsed_command, pp_stdout_file, pp_stderr_file)

            pp_stdout_file.close()
            pp_stderr_file.close()
            pp_command_file.close()
            os.chdir("..")


splash_postprocessing = SplashPostprocessing(postprocessing_item)

for executable_path in [dsv_manip, satc, satc_dump, satc_merge, sig_anch, read_selector, compactors, lookup_table]:
    p = subprocess.Popen(executable_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode("utf-8")
    err = err.decode("utf-8")
    out += err
    version_is_OK = False
    for line in out.split('\n'):
        if line.startswith("splash version:"):
            ver = line.split()[-1]
            if ver != SPLASH_VERSION:
                print(f"Error: {executable_path} version is {ver} while this script is {SPLASH_VERSION}. Please use the same version.")
                sys.exit(1)
            else:
                version_is_OK=True
                break
    if not version_is_OK:
        print(f"Error: cannot find splash version number for {executable_path}.")
        sys.exit(1)


with open(f"{logs_dir}/splash-cmd.log", "w") as f:
    f.writelines(" ".join(sys.argv))

inputs = []

sample_name_to_id_file = open(sample_name_to_id, "w")

def remove_kmc_output(path):
    os.remove(f"{path}.kmc_pre")
    os.remove(f"{path}.kmc_suf")

with open(input_file) as f:
    for line in f:
        sample_id, path = line.strip().split()
        inputs.append([path, sample_id])


def get_extension_gz_aware(path):
    base, ext = os.path.splitext(path)
    if ext == ".gz":
        _, pre_gz = os.path.splitext(base)
        ext = pre_gz + ext
    return ext[1:] # without leading '.'


def get_file_format(path):
    extensions_to_format = {
        "fastq": "fq",
        "fastq.gz": "fq",
        "fq": "fq",
        "fq.gz": "fq",

        "fasta": "fa",
        "fasta.gz": "fa",
        "fa": "fa",
        "fa.gz": "fa",

        "bam": "fbam"
    }
    ext = get_extension_gz_aware(path)
    if ext not in extensions_to_format.keys():
        print(f"Error: Unknown extension of {path}.")
        sys.exit(1)

    return extensions_to_format[ext]

def get_file_format_10x_visium(path):
    with open(path) as f:
        for line in f:
            R1, R2 = line.strip().split(",")
            R1 = R1.strip()
            R2 = R2.strip()

            R1_format = get_file_format(R1)
            R2_format = get_file_format(R2)

            if R1_format != R2_format:
                print(f"Files {R1} and {R2} are in different formats ({R1_format} and {R2_format})")
                sys.exit(1)
            return R1_format

def remove_bkc_input(path):
    with open(path) as f:
        for line in f:
            R1, R2 = line.strip().split(",")
            R1 = R1.strip()
            R2 = R2.strip()

            os.remove(R1)
            os.remove(R2)

def crash_on_wrong_format(file_format, supported, err_msg):
    if not file_format in supported:
        print(err_msg)
        sys.exit(1)

compactors_cfg = {
    "max_length": anchor_len * 6,
    "num_threads": max_cpus_to_use_in_auto_adjust,
    "epsilon": 0.001,
    "beta": 0.5,
    "num_kmers": 2,
}

lookup_table_cfg = {

}

#to have compatibile name with compactors, but this is not parameter, because by default we dont run lookup_tables (at least for now)
#if we decide to run by default (hard to say how, because it require some config) we may change this to parameter
#for now I will just set it to False if config was specified and seems to be correct
without_lookup_table = lookup_table_config == ""

# verify if file format is supported
for path, sample_id in inputs:
    file_format_for_10x_visium = ["fq", "fa"]
    file_format_for_base = ["fq", "fa", "fbam"]
    file_format_for_base_err_msg = ["fq", "fa", "bam"] # because I need fbam for kmc but to the user I want display just bam...

    err_10x_visium = f"Error: for --technology 10x or visium only following input file formats are supported: "+ ", ".join(file_format_for_10x_visium)
    err_base = f"Error: for --technology base only following input file formats are supported: "+ ", ".join(file_format_for_base_err_msg)

    if is_10x_or_visium():
        format = get_file_format_10x_visium(path)
        crash_on_wrong_format(format, file_format_for_10x_visium, err_10x_visium)
    else:
        format = get_file_format(path)
        crash_on_wrong_format(format, file_format_for_base, err_base)

    if format == "fa":
        compactors_cfg["input_format"] = "fasta"
    elif format == "fq":
        compactors_cfg["input_format"] = "fastq"
    else:
        without_compactors = True
        print(f"Warning: file format ({format}) is not supported by compactors -> forcing --without_compactors")

if is_10x_or_visium() and not without_compactors and not export_filtered_input:
    print("Warning: when --without_compactors is not used for 10x/visium data --export_filtered_input is forced")
    export_filtered_input = True

if keep_top_n_target_entropy == 0 and not without_compactors:
    keep_top_n_target_entropy = 10000
    print(f"Warning: when --without_compactors is not used --keep_top_n_target_entropy must be different than zero, auto adsjuting to {keep_top_n_target_entropy}")

if keep_top_n_target_entropy == 0 and splash_postprocessing.RequireTopTargetEntropyExtendors():
    keep_top_n_target_entropy = 10000
    print(f"Warning: one of postprocessings require --keep_top_n_target_entropy different than zero, auto adsjuting to {keep_top_n_target_entropy}")

if keep_top_n_effect_size_bin == 0 and not without_compactors:
    keep_top_n_effect_size_bin = 20000
    print(f"Warning: when --without_compactors is not used --keep_top_n_effect_size_bin must be different than zero, auto adsjuting to {keep_top_n_effect_size_bin}")

if keep_top_n_effect_size_bin == 0 and splash_postprocessing.RequireTopEffectSizeExtendors():
    keep_top_n_effect_size_bin = 20000
    print(f"Warning: one of postprocessings require --keep_top_n_effect_size_bin different than zero, auto adsjuting to {keep_top_n_effect_size_bin}")

if keep_top_n_target_entropy == 0 and not without_lookup_table:
    keep_top_n_target_entropy = 10000
    print(f"Warning: when lookup table is used --keep_top_n_target_entropy must be different than zero, auto adsjuting to {keep_top_n_target_entropy}")

if keep_top_n_effect_size_bin == 0 and not without_lookup_table:
    keep_top_n_effect_size_bin = 20000
    print(f"Warning: when lookup table is used --keep_top_n_effect_size_bin must be different than zero, auto adsjuting to {keep_top_n_effect_size_bin}")

if keep_top_n_target_entropy == 0 and keep_top_target_entropy_anchors_satc:
    keep_top_n_target_entropy = 10000
    print(f"Warning: when --keep_top_target_entropy_anchors_satc is used --keep_top_n_target_entropy must be different than zero, auto adsjuting to {keep_top_n_target_entropy}")

if keep_top_n_effect_size_bin == 0 and keep_top_effect_size_bin_anchors_satc:
    keep_top_n_effect_size_bin = 20000
    print(f"Warning: when --keep_top_effect_size_bin_anchors_satc is used --keep_top_n_effect_size_bin must be different than zero, auto adsjuting to {keep_top_n_target_entropy}")

if not without_compactors and compactors_config != "":
    try:
        with open(compactors_config) as f:
            compactors_cfg.update(json.load(f))
    except:
        print(f"Error: cannot read compactors config from {compactors_config}")
        sys.exit(1)

if lookup_table_config != "":
    try:
        with open(lookup_table_config) as f:
            lookup_table_cfg.update(json.load(f))
            if not "lookup_tables" in lookup_table_cfg:
                print(f"Error: {lookup_table_config} does not contain required field \"lookup_tables\" (array of lookup tables)")
                sys.exit(1)

            for item in lookup_table_cfg["lookup_tables"]:
                if not "name" in item or not "path" in item:
                    print(f"Error: for each lookup table following fields are required: \"name\", \"path\", problem is with {item}")
                    sys.exit(1)

                name = item["name"]
                path = item["path"]

                invalud_chars = '<>:"/\\|?*'
                if any(c in set(invalud_chars) for c in name):
                    print(f"Error: lookup_table name ({name}) is invalid (contains one of '{invalud_chars}')")
                if not os.path.isfile(path):
                    print(f"Error: lookup_table {name} path ({path}) is not a file")
                    sys.exit(1)
    except:
        print(f"Error: cannot read lookup_table config from {lookup_table_config}")
        sys.exit(1)

class Timer:
    def __init__(self):
        self.timings = {}
        self.time_point = time.perf_counter()
        self.total_time = 0

    def start(self):
        self.time_point = time.perf_counter()
        self.timings = {}
        self.total_time = 0

    def catch(self, desc):
        t2 = time.perf_counter()
        time_span = t2 - self.time_point
        self.time_point = t2
        self.total_time += time_span
        self.timings[desc] = time_span
        return time_span

    def get_total(self):
        return self.total_time()

    def get_timings(self):
        return self.timings, self.total_time

def calc_dir_size_bytes(path):
    total = 0
    for entry in os.scandir(path):
        if entry.is_file(follow_symlinks=False):
            total += entry.stat().st_size
        elif entry.is_dir(follow_symlinks=False):
            total += calc_dir_size_bytes(entry.path)
    return total

def calc_dir_size(path):
    return calc_dir_size_bytes(path) / 1024 / 1024 # MiB


#TODO: can be done more performant, in some cases probably these number may be just
#implementation based on: https://stackoverflow.com/questions/845058/how-to-get-the-line-count-of-a-large-file-cheaply-in-python
def get_number_of_lines(path):
    with open(path, 'rb') as f:
        lines = 0
        buf_size = 1024 * 1024
        read_f = f.raw.read

        buf = read_f(buf_size)
        while buf:
            lines += buf.count(b'\n')
            buf = read_f(buf_size)

        return lines

class MemoryLogger:
    def __init__(self) -> None:
        self.logs = {}
    def report(self, name, value):
        self.logs[name] = value

timer = Timer()
memory_logger = MemoryLogger()

# to be called at the end but also ocassionaly, such that if there were a crash at least partial results are present
def store_current_time_disk_ram_logs():
    #prefix T_ means time in sec
    #prefix S_ means size on disk in MiB
    #prefix M_ means memory in MiB
    #prefix N_ means a number
    timings_memory_disk = {}
    timings_memory_disk.update(timer.timings)
    timings_memory_disk.update(memory_logger.logs)

    with open(f"{logs_dir}/stats.json", "w") as f:
        f.write(json.dumps(timings_memory_disk, indent="\t"))



def get_max_maximum_resident_set_size(path):
    with open(path) as f:
        res = 0
        for line in f:
            line = line.strip()
            #this seems to be in KiB
            if "Maximum resident set size (kbytes): " in line:
                mrss = int(line.split(":")[1])
                if mrss > res:
                    res = mrss
        return res / 1024 #in MiB

def get_max_maximum_resident_set_size_many_files(paths):
    res = 0
    for path in paths:
        tmp = get_max_maximum_resident_set_size(path)
        if tmp > res:
            res = tmp
    return res

def get_sum_maximum_resident_set_size_many_files(paths):
    res = 0
    for path in paths:
        res += get_max_maximum_resident_set_size(path)
    return res

###############################################################################
# STAGE 1
# Construct (sample, anchor, target, count) records per each sample
# per each sample there is {n_bins} output files
###############################################################################

print("Starting stage 1")
print("Current time:", get_cur_time(), flush=True)

if export_filtered_input and not is_10x_or_visium():
    print("Warning: exporting filtered input has only effect in case of 10x/visium data")

export_filtered_input_dir = f"{outname_prefix}_filtered_input"
if export_filtered_input and not os.path.exists(export_filtered_input_dir):
    os.makedirs(export_filtered_input_dir)

cbc_dir = ""
if export_cbc_logs and is_10x_or_visium():
    cbc_dir = f"{outname_prefix}_cbc"
    if not os.path.exists(cbc_dir):
        os.makedirs(cbc_dir)

#251 because the max k for kmc is 255, but for efficient use of (k,x)-mer k should be at most 32*n-4 (x=3, 2 bits (1 symbol) for x)
reduce_gap_len = False
if not is_10x_or_visium() and anchor_len + gap_len + target_len > 251:
    reduce_gap_len = True
    reduced_gap_len = 50
    x = 4 # in kmc max_x = 3 but there are 2 more bits (1 symbol) for the value of x= 4, or if newer kmc will be used the max x = 4 and there is not need for additional symbols
    while (anchor_len + reduced_gap_len + target_len + 2 * x) % 32 != 0:
        reduced_gap_len+=1
    print(f"INFO: because the gap_len ({gap_len}) is quite large the input data will be preprocessed such that lower gap len ({reduced_gap_len}) may be used")

def stage_1_task(id, input, out, err):
    _artifacts_param = f"--artifacts {artifacts}" if artifacts != "" else ""

    if is_10x_or_visium():
        fname = input[0]
        sample_name = input[1]

        _apply_filter_illumina_adapters_param = "--apply_filter_illumina_adapters" if not dont_filter_illumina_adapters else ""
        _predefined_cbc_param = f"--predefined_cbc {predefined_cbc}" if predefined_cbc != "" else ""
        _log_path_param = f"--log_name {cbc_dir}/{sample_name}" if export_cbc_logs else ""
        _filtered_input_path_param = f"--filtered_input_path {export_filtered_input_dir}" if export_filtered_input else ""
        _allow_strange_cbc_umi_reads_param = f"--allow_strange_cbc_umi_reads" if allow_strange_cbc_umi_reads else ""

        file_format = get_file_format_10x_visium(fname)
        if file_format == "fa":
            file_format = "fasta"
        elif file_format == "fq":
            file_format = "fastq"
        else:
            print(f"Error: unknown file format {file_format} for bkc")
            sys.exit(1)
        cmd = f"{bkc} \
            --mode pair \
            --input_format {file_format} \
            --n_threads {n_threads_stage_1_internal} \
            --leader_len {anchor_len} \
            --follower_len {target_len} \
            --gap_len {gap_len} \
            --cbc_len {cbc_len} \
            --umi_len {umi_len} \
            --soft_cbc_umi_len_limit {soft_cbc_umi_len_limit} \
            --n_splits {n_bins} \
            --cbc_filtering_thr {cbc_filtering_thr} \
            --output_format splash \
            --leader_sample_counts_threshold {anchor_sample_counts_threshold} \
            --poly_ACGT_len {poly_ACGT_len} \
            --technology {technology} \
            {_filtered_input_path_param} \
            --export_filtered_input_mode second \
            {_artifacts_param} \
            {_apply_filter_illumina_adapters_param} \
            {_predefined_cbc_param} \
            {_allow_strange_cbc_umi_reads_param} \
            {_log_path_param} \
            --input_name {fname} \
            --verbose 2 \
            --output_name {tmp_dir}/{sample_name} \
            --sample_id {id}"
        run_cmd(cmd, out, err)

        if REMOVE_INPUT:
            remove_bkc_input(fname)

    else:    
        fname = input[0]
        sample_name = input[1]
        
        kmc_dir_tmp_name = f"{tmp_dir}/kmc_tmp_{sample_name}"

        if not os.path.exists(kmc_dir_tmp_name):
            os.makedirs(kmc_dir_tmp_name)

        ram_only_param = ""
        if kmc_use_RAM_only_mode:
            ram_only_param="-r"

        file_format = get_file_format(fname)

        gap_len_for_kmc = gap_len
        fname_for_kmc = fname
        reduce_gap_files = []
        if reduce_gap_len:
            #I assume each input filename is different
            base_name = os.path.basename(fname)
            if base_name.endswith(".gz"):
                base_name = base_name[:-len(".gz")]
            base_name = os.path.splitext(base_name)[0]

            fname_reduced_gap_base = f"{kmc_dir_tmp_name}/reduced_{base_name}"

            cmd = f"{gap_shortener} \
                --anchor_len {anchor_len} \
                --gap_len {gap_len} \
                --new_gap_len {reduced_gap_len} \
                --target_len {target_len} \
                --n_output_files {n_threads_stage_1_internal} \
                --n_threads {n_threads_stage_1_internal} \
                --compression_level 3 \
                {fname} \
                {fname_reduced_gap_base}"
            run_cmd(cmd, out, err)

            for i in range(1, n_threads_stage_1_internal + 1):
                reduce_gap_files.append(fname_reduced_gap_base + f"_{i}.fa.gz")

            fname_for_kmc = f"{kmc_dir_tmp_name}/reduced_list.txt"
            with open(fname_for_kmc, "w") as f:
                for p in reduce_gap_files:
                    f.write(f"{p}\n")
            fname_for_kmc = "@" + fname_for_kmc

            file_format = "fa"
            gap_len_for_kmc = reduced_gap_len

        cmd = f"{kmc} -t{n_threads_stage_1_internal} -{file_format} -b -ci1 -cs65535 -k{anchor_len + gap_len_for_kmc + target_len} {ram_only_param} -m{kmc_max_mem_GB} {fname_for_kmc} {tmp_dir}/{sample_name} {kmc_dir_tmp_name}"
        run_cmd(cmd, out, err)

        if clean_up:
            if reduce_gap_len:
                os.remove(f"{kmc_dir_tmp_name}/reduced_list.txt")
                for x in reduce_gap_files:
                    os.remove(x)
            os.rmdir(kmc_dir_tmp_name)

        cmd = f"{kmc_tools} -t{n_threads_stage_1_internal} transform {tmp_dir}/{sample_name} sort {tmp_dir}/{sample_name}.sorted"
        run_cmd(cmd, out, err)

        # if we have small k kmc_tools will omit sorting because it is already sorted, lets just copy
        if not os.path.exists(f"{tmp_dir}/{sample_name}.sorted.kmc_pre"):
            shutil.copy(f"{tmp_dir}/{sample_name}.kmc_pre", f"{tmp_dir}/{sample_name}.sorted.kmc_pre")
            shutil.copy(f"{tmp_dir}/{sample_name}.kmc_suf", f"{tmp_dir}/{sample_name}.sorted.kmc_suf")

        if clean_up:
            remove_kmc_output(f"{tmp_dir}/{sample_name}")
        
        _dont_filter_illumina_adapters_param = "--dont_filter_illumina_adapters" if dont_filter_illumina_adapters else ""
        cmd = f"{satc} \
            --anchor_len {anchor_len} \
            --target_len {target_len} \
            --n_bins {n_bins} \
            --anchor_sample_counts_threshold {anchor_sample_counts_threshold} \
            --min_hamming_threshold {min_hamming_threshold} \
            --poly_ACGT_len {poly_ACGT_len} \
            {_artifacts_param} \
            {_dont_filter_illumina_adapters_param} \
            {tmp_dir}/{sample_name} \
            {tmp_dir}/{sample_name}.sorted {id}"
        run_cmd(cmd, out, err)
        
        if clean_up:
            remove_kmc_output(f"{tmp_dir}/{sample_name}.sorted")


stage_1_queue = queue.Queue()
def stage_1_worker(tid):
    tid = str(tid)
    while len(tid) < 4:
        tid = "0" + tid
    
    stdoutfile = open(f"{logs_dir}/stage_1_thread-{tid}.log", "w")
    while True:
        item = stage_1_queue.get()
        if item is None:
            break
        id = item[0]
        input = item[1]

        action_at_function_exit(lambda: stage_1_queue.task_done()) # will be called even if stage_1_task will raise exception

        try:
            stage_1_task(id, input, stdoutfile, stdoutfile)
        except StopAllBecauseCrititcalError: #just consume next task
            pass

    stdoutfile.close()


stage_1_threads = []
for i in range(n_threads_stage_1):
    t = threading.Thread(target=stage_1_worker, args=(i,))
    t.start()
    stage_1_threads.append(t)

for id, input in enumerate(inputs):
    sample_name = input[1]
    sample_name_to_id_file.write(f"{sample_name} {id}\n")
    stage_1_queue.put((id, input))

stage_1_queue.join()

# stop workers
for i in range(n_threads_stage_1):
    stage_1_queue.put(None)
for t in stage_1_threads:
    t.join()
check_and_handle_error()

stage_1_logs = []
for i in range(n_threads_stage_1):
    tid = str(i)
    while len(tid) < 4:
        tid = "0" + tid
    stage_1_logs.append(f"{logs_dir}/stage_1_thread-{tid}.log")

memory_logger.report("M_stage_1", get_sum_maximum_resident_set_size_many_files(stage_1_logs))
memory_logger.report("S_stage_1", calc_dir_size(tmp_dir))

if export_filtered_input and is_10x_or_visium():
    memory_logger.report(f"S_{export_filtered_input_dir}", calc_dir_size(export_filtered_input_dir))

sample_name_to_id_file.close()
if len(anchor_list):    
    anchor_list_param = f"--anchor_list {anchor_list}"
else:
    anchor_list_param=""

###############################################################################
# STAGE 2
# For each bin merge all samples
# output: file {n_bins} files, where each have
# its (sample, anchor, target, count)
# this stage may be is used to compute pvals and other stats
# and optionaly to dump merged files (in satc or text format)
###############################################################################
timer.catch("T_stage_1")
store_current_time_disk_ram_logs()
print("Stage 1 done")
print("Starting stage 2")
print("Current time:", get_cur_time(), flush=True)

Cjs_dir = f"{outname_prefix}_Cjs"
if dump_Cjs and not os.path.exists(Cjs_dir):
    os.makedirs(Cjs_dir)

if dump_sample_anchor_target_count_txt:
    dump_dir = f"{outname_prefix}_dumps"
    if not os.path.exists(dump_dir):
        os.makedirs(dump_dir)

keep_satc_from_satc_merge = dump_sample_anchor_target_count_binary or keep_significant_anchors_satc or keep_top_target_entropy_anchors_satc or keep_top_effect_size_bin_anchors_satc

if keep_satc_from_satc_merge:
    satc_dir = f"{outname_prefix}_satc"
    if not os.path.exists(satc_dir):
        os.makedirs(satc_dir)

#must be defined for post-processing
compactors_dir = f"{outname_prefix}_compactors"
if not without_compactors:
    if not os.path.exists(compactors_dir):
        os.makedirs(compactors_dir)

if not without_lookup_table:
    lookup_table_dir = f"{outname_prefix}_lookup_table"
    if not os.path.exists(lookup_table_dir):
        os.makedirs(lookup_table_dir)

def stage_2_task(bin_id, out, err):
    satc_merge_inputs = []

    for input in inputs:
        sample_name = input[1]

        if is_10x_or_visium():
            if n_bins > 1:
                name = f"{tmp_dir}/{sample_name}.{bin_id}"
            else:
                name = f"{tmp_dir}/{sample_name}"
        else:
            name = f"{tmp_dir}/{sample_name}.{bin_id}.bin"

        satc_merge_inputs.append(name)
        
    in_file_name = f"{tmp_dir}/files.bin{bin_id}.lst"
    with open(in_file_name, "w") as f:
        for x in satc_merge_inputs:
            f.write(f"{x}\n")
    
    _without_alt_max_param = "--without_alt_max" if without_alt_max else ""
    _without_sample_spectral_sum_param = "--without_sample_spectral_sum" if without_sample_spectral_sum else ""
    _with_effect_size_cts_param = "--with_effect_size_cts" if with_effect_size_cts else ""
    _with_pval_asymp_opt_param = "--with_pval_asymp_opt" if with_pval_asymp_opt else ""
    _without_seqence_entropy_param = "--without_seqence_entropy" if without_seqence_entropy else ""
    _compute_also_old_base_pvals_param = "--compute_also_old_base_pvals" if compute_also_old_base_pvals else ""

    _cjs_out_param = f"--cjs_out {Cjs_dir}/bin{bin_id}.cjs.gz" if dump_Cjs else ""
    _cjs_without_header_param = "--cjs_without_header" if bin_id != 0 else ""

    _cell_type_samplesheet_param = f"--cell_type_samplesheet {cell_type_samplesheet}" if cell_type_samplesheet != "" else ""
    _Cjs_samplesheet_param = f"--Cjs_samplesheet {Cjs_samplesheet}" if Cjs_samplesheet != "" else ""


    _dump_sample_anchor_target_count_txt_param = ""
    if dump_sample_anchor_target_count_txt:
        dump_dir = f"{outname_prefix}_dumps"
        _dump_sample_anchor_target_count_txt_param = f"--dump_sample_anchor_target_count_txt {dump_dir}/bin{bin_id}.satc.dump"

    _dump_sample_anchor_target_count_binary_param = ""
    if keep_satc_from_satc_merge:
        satc_dir = f"{outname_prefix}_satc"
        _dump_sample_anchor_target_count_binary_param = f"--dump_sample_anchor_target_count_binary {satc_dir}/bin{bin_id}.satc"

    cmd=f"{satc_merge} \
    --technology {technology} \
    {_without_alt_max_param} \
    {_without_sample_spectral_sum_param} \
    {_with_effect_size_cts_param} \
    {_with_pval_asymp_opt_param} \
    {_without_seqence_entropy_param} \
    {_compute_also_old_base_pvals_param} \
    {_cjs_out_param} \
    --max_pval_opt_for_Cjs {max_pval_opt_for_Cjs} \
    {_cjs_without_header_param} \
    --anchor_count_threshold {anchor_count_threshold} \
    --anchor_samples_threshold {anchor_samples_threshold} \
    --anchor_unique_targets_threshold {anchor_unique_targets_threshold} \
    --n_most_freq_targets {n_most_freq_targets} \
    --n_most_freq_targets_for_stats {n_most_freq_targets_for_stats} \
    --n_most_freq_targets_for_dump {n_most_freq_targets_for_dump} \
    --opt_num_inits {opt_num_inits} \
    --opt_num_iters {opt_num_iters} \
    --num_rand_cf {num_rand_cf} \
    --num_splits {num_splits} \
    --opt_train_fraction {opt_train_fraction} \
    {_dump_sample_anchor_target_count_txt_param} \
    {_dump_sample_anchor_target_count_binary_param} \
    --sample_names {sample_name_to_id} \
    --format {satc_merge_dump_format} \
    {_cell_type_samplesheet_param} \
    {_Cjs_samplesheet_param} \
    {anchor_list_param} \
    {tmp_dir}/bin{bin_id}.stats.tsv {in_file_name}"
    run_cmd(cmd, out, err)

    if clean_up:
        os.remove(in_file_name)
        for f in satc_merge_inputs:
            os.remove(f)


stage_2_queue = queue.Queue()
def stage_2_worker(tid):
    tid = str(tid)
    while len(tid) < 4:
        tid = "0" + tid
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)

    stdoutfile = open(f"{logs_dir}/stage_2_thread-{tid}.log", "w")
    while True:
        item = stage_2_queue.get()       
        if item is None:
            break
        bin_id = item

        action_at_function_exit(lambda: stage_2_queue.task_done()) # will be called even if stage_2_task will raise exception

        try:
            stage_2_task(bin_id, stdoutfile, stdoutfile)
        except StopAllBecauseCrititcalError: #just consume next task
            pass

    stdoutfile.close()

stage_2_threads = []
for i in range(n_threads_stage_2):
    t = threading.Thread(target=stage_2_worker, args=(i,))
    t.start()
    stage_2_threads.append(t)

for bin_id in range(n_bins):
    stage_2_queue.put(bin_id)

stage_2_queue.join()
    
# stop workers
for i in range(n_threads_stage_2):
    stage_2_queue.put(None)
for t in stage_2_threads:
    t.join()
check_and_handle_error()

stage_2_logs = []
for i in range(n_threads_stage_2):
    tid = str(i)
    while len(tid) < 4:
        tid = "0" + tid
    stage_2_logs.append(f"{logs_dir}/stage_2_thread-{tid}.log")

memory_logger.report("M_stage_2", get_sum_maximum_resident_set_size_many_files(stage_2_logs))
memory_logger.report("S_stage_2_out", calc_dir_size(tmp_dir))

if dump_sample_anchor_target_count_txt:
    memory_logger.report("S_stage_2_satc_dump", calc_dir_size(dump_dir))

if keep_satc_from_satc_merge:
    memory_logger.report("S_stage_2_satc", calc_dir_size(satc_dir))

#merge cjs
if dump_Cjs:
    with open(f"{Cjs_dir}/cjs.tsv.gz", "wb") as all:
        for bin_id in range(n_bins):
            cjs_file_path = f"{Cjs_dir}/bin{bin_id}.cjs.gz"
            with open (cjs_file_path, "rb") as cjs_file:
                shutil.copyfileobj(cjs_file, all)
            if clean_up:
                os.remove(cjs_file_path)
    memory_logger.report(f"S_{Cjs_dir}", calc_dir_size(Cjs_dir))


timer.catch("T_stage_2")
store_current_time_disk_ram_logs()
print("Stage 2 done")
print("Starting stage 3")
print("Current time:", get_cur_time(), flush=True)

files_for_correction = open(f"{tmp_dir}/files_for_correction.list.txt", "w")
for bin_id in range(n_bins):
    files_for_correction.write(f"{tmp_dir}/bin{bin_id}.stats.tsv" + "\n")
files_for_correction.close()


cmd=f"{sig_anch} \
--samplesheet {input_file} \
--infile_bins {tmp_dir}/files_for_correction.list.txt \
--outfile_scores {outname_prefix}.after_correction.scores.tsv \
--outfile_all_anchors_pvals {outname_prefix}.after_correction.all_anchors.tsv \
--fdr_threshold {fdr_threshold} \
--col_name {pvals_correction_col_name}"

log_file = f"{logs_dir}/correction.log"
with open(log_file, "w") as stdoutfile:
    run_cmd(cmd, stdoutfile, stdoutfile)

memory_logger.report(f"N_all_anchors", get_number_of_lines(f"{outname_prefix}.after_correction.all_anchors.tsv"))
memory_logger.report(f"N_scores", get_number_of_lines(f"{outname_prefix}.after_correction.scores.tsv"))

if clean_up:
    for bin_id in range(n_bins):
        os.remove(f"{tmp_dir}/bin{bin_id}.stats.tsv")
    os.remove(f"{tmp_dir}/files_for_correction.list.txt")
check_and_handle_error()

memory_logger.report("M_stage_3", get_max_maximum_resident_set_size(log_file))

timer.catch("T_stage_3")
store_current_time_disk_ram_logs()
print("Stage 3 done")


if keep_top_n_target_entropy != 0:
    print("Start top target entropy select")
    print("Current time:", get_cur_time(), flush=True)

    cmd=f"{dsv_manip} limit \
        --select highest \
        --sort \
        --sort_order desc \
        --value_type double \
        {keep_top_n_target_entropy} \
        target_entropy \
        {outname_prefix}.after_correction.scores.tsv \
        {outname_prefix}.after_correction.scores.top_target_entropy.tsv"

    log_file = f"{logs_dir}/keep_top_n_target_entropy.log"
    with open(log_file, "w") as stdoutfile:
        run_cmd(cmd, stdoutfile, stdoutfile)

    memory_logger.report("M_select_top_target_entropy", get_max_maximum_resident_set_size(log_file))

    timer.catch("T_select_top_target_entropy")
    store_current_time_disk_ram_logs()
    print("Top target entropy select done")

if keep_top_n_effect_size_bin != 0:
    print("Start top effect_size_bin select")
    print("Current time:", get_cur_time(), flush=True)

    cmd=f"{dsv_manip} limit \
        --select highest \
        --sort \
        --sort_order desc \
        --value_type double \
        {keep_top_n_effect_size_bin} \
        effect_size_bin \
        {outname_prefix}.after_correction.scores.tsv \
        {outname_prefix}.after_correction.scores.top_effect_size_bin.tsv"

    log_file = f"{logs_dir}/keep_top_n_effect_size_bin.log"
    with open(log_file, "w") as stdoutfile:
        run_cmd(cmd, stdoutfile, stdoutfile)

    memory_logger.report("M_select_top_effect_size_bin", get_max_maximum_resident_set_size(log_file))
    timer.catch("T_select_top_effect_size_bin")
    store_current_time_disk_ram_logs()
    print("Effect size bin select done")

def keep_satc_for(anchor_tsv_path, out_path, log_file_path):
    satc_list = open(f"{tmp_dir}/satc_list.txt", "w")
    for bin_id in range(n_bins):
        satc_list.write(f"{satc_dir}/bin{bin_id}.satc" + "\n")
    satc_list.close()

    cmd = f"{satc_dump} \
    --n_bins {n_bins} \
    --anchor_list {anchor_tsv_path} \
    --binary \
    {tmp_dir}/satc_list.txt \
    {out_path}"

    with open(f"{log_file_path}", "w") as stdoutfile:
        run_cmd(cmd, stdoutfile, stdoutfile)

    if clean_up:
        os.remove(f"{tmp_dir}/satc_list.txt")

if keep_significant_anchors_satc:
    print("Start keep SATC for significant anchors")
    print("Current time:", get_cur_time(), flush=True)

    log_file = f"{logs_dir}/keep_significant_anchors_satc.log"
    keep_satc_for(f"{outname_prefix}.after_correction.scores.tsv",
                f"{satc_dir}/after_correction.scores.satc",
                log_file)

    memory_logger.report("M_keep_significant_anchors_satc", get_max_maximum_resident_set_size(log_file))
    timer.catch("T_keep_significant_anchors_satc")
    store_current_time_disk_ram_logs()
    print("keep SATC for significant anchors done")

if keep_top_target_entropy_anchors_satc:
    print("Start keep SATC for top target entropy anchors")
    print("Current time:", get_cur_time(), flush=True)

    log_file = f"{logs_dir}/keep_top_target_entropy_anchors_satc.log"
    keep_satc_for(f"{outname_prefix}.after_correction.scores.top_target_entropy.tsv",
                f"{satc_dir}/after_correction.scores.top_target_entropy.satc",
                log_file)

    memory_logger.report("M_keep_top_target_entropy_anchors_satc", get_max_maximum_resident_set_size(log_file))
    timer.catch("T_keep_top_target_entropy_anchors_satc")
    store_current_time_disk_ram_logs()
    print("keep SATC for top target entropy anchors done")


if keep_top_effect_size_bin_anchors_satc:
    print("Start keep SATC for top effect size bin anchors")
    print("Current time:", get_cur_time(), flush=True)

    log_file = f"{logs_dir}/keep_top_effect_size_bin_anchors_satc.log"
    keep_satc_for(f"{outname_prefix}.after_correction.scores.top_effect_size_bin.tsv",
                f"{satc_dir}/after_correction.scores.top_effect_size_bin.satc",
                log_file)

    memory_logger.report("M_keep_top_effect_size_bin_anchors_satc", get_max_maximum_resident_set_size(log_file))
    timer.catch("T_keep_top_effect_size_bin_anchors_satc")
    store_current_time_disk_ram_logs()
    print("keep SATC for top effect size bin anchors done")

#it means we keep them to be able to dump significant or high entropy anchors only, but the user dont want to keep all SATCs
if keep_satc_from_satc_merge and not dump_sample_anchor_target_count_binary and clean_up:
    for bin_id in range(n_bins):
        os.remove(f"{satc_dir}/bin{bin_id}.satc")

if keep_satc_from_satc_merge:
    memory_logger.report(f"S_{satc_dir}", calc_dir_size(satc_dir))

if supervised_test_samplesheet != "":
    _datatype_param = "non10x" if technology == "base" else "10x"
    print("Starting supervised_test")
    print("Current time:", get_cur_time(), flush=True)

    cmd = f"Rscript {supervised_test} \
        {n_bins} \
        {sample_name_to_id} \
        {satc_dump} \
        {outname_prefix}_satc \
        {outname_prefix}.after_correction.scores.tsv \
        . \
        {supervised_test_samplesheet} \
        supervised_test \
        {_datatype_param} \
        {supervised_test_anchor_sample_fraction_cutoff} \
        {supervised_test_num_anchors}"

    log_file = f"{logs_dir}/supervised_test.log"
    stdoutfile = open(log_file, "w")
    run_cmd(cmd, stdoutfile, stdoutfile)
    stdoutfile.close()

    memory_logger.report("M_supervised_test", get_max_maximum_resident_set_size(log_file))
    timer.catch("T_supervised_test")
    store_current_time_disk_ram_logs()
    print("supervised_test done")

if not without_compactors:
    print("Starting compactors")
    print("Current time:", get_cur_time(), flush=True)
    params = ""
    for (key, value) in compactors_cfg.items():
        params += " --" + key + " " + str(value)

    with open(f"{tmp_dir}/compactors_input.txt", "w") as compactors_input_path:
        if is_10x_or_visium():
            with open(input_file) as f:
                for line in f:
                    path = line.strip().split()[1]
                    with open(path) as f2:
                        for line2 in f2:
                            x = line2.strip().split(",")[1]
                            x = os.path.split(x)[-1]
                            ext = get_extension_gz_aware(x)
                            x = x[:-len(ext)] + "dedup." + ext
                            print(os.path.join(export_filtered_input_dir, x), file=compactors_input_path)
        else:
            with open(input_file) as f:
                for line in f:
                    path = line.strip().split()[1]
                    print(path, file=compactors_input_path)


    cmd_for_top_target_entropy = f"{compactors} \
        --keep_temp \
        {params} \
        {tmp_dir}/compactors_input.txt \
        {outname_prefix}.after_correction.scores.top_target_entropy.tsv \
        {compactors_dir}/after_correction.scores.top_target_entropy.tsv"

    logs_files = []
    logs_files.append(f"{logs_dir}/compactors.top_target_entropy.log")
    with open(logs_files[-1], "w") as stdoutfile:
        run_cmd(cmd_for_top_target_entropy, stdoutfile, stdoutfile)

    cmd_for_top_effect_size_bin = f"{compactors} \
        --keep_temp \
        {params} \
        {tmp_dir}/compactors_input.txt \
        {outname_prefix}.after_correction.scores.top_effect_size_bin.tsv \
        {compactors_dir}/after_correction.scores.top_effect_size_bin.tsv"

    logs_files.append(f"{logs_dir}/compactors.top_effect_size_bin.log")
    with open(logs_files[-1], "w") as stdoutfile:
        run_cmd(cmd_for_top_effect_size_bin, stdoutfile, stdoutfile)

    if clean_up:
        os.remove(f"{tmp_dir}/compactors_input.txt")
        memory_logger.report("S_tmp_compactors", calc_dir_size("tmp-compactors"))
        shutil.rmtree("tmp-compactors")

    memory_logger.report("M_compactors", get_max_maximum_resident_set_size_many_files(logs_files))
    memory_logger.report(f"S_{compactors_dir}", calc_dir_size(compactors_dir))
    timer.catch("T_compactors")
    store_current_time_disk_ram_logs()
    print("compactors done")

if not without_lookup_table:
    print("Starting lookup_table")
    print("Current time:", get_cur_time(), flush=True)

    logs_files = []
    for item in lookup_table_cfg["lookup_tables"]:
        name = item["name"]
        path = item["path"]
        with open(f"{tmp_dir}/lookup_table_{name}_input.txt", "w") as f:
            print(f"--input_fmt extendors --output_fmt extendors --report_fmt plain --stats_fmt with_stats \
                    {outname_prefix}.after_correction.scores.top_target_entropy.tsv \
                    {lookup_table_dir}/after_correction.scores.top_target_entropy.{name}.tsv", file=f)

            print(f"--input_fmt extendors --output_fmt extendors --report_fmt plain --stats_fmt with_stats \
                    {outname_prefix}.after_correction.scores.top_effect_size_bin.tsv \
                    {lookup_table_dir}/after_correction.scores.top_effect_size_bin.{name}.tsv", file=f)

            if not without_compactors:
                print(f"--input_fmt compactors --output_fmt compactors --report_fmt plain --stats_fmt with_stats \
                        {compactors_dir}/after_correction.scores.top_target_entropy.tsv \
                        {lookup_table_dir}/after_correction.scores.top_target_entropy.{name}.compactors.tsv", file=f)

                print(f"--input_fmt compactors --output_fmt compactors --report_fmt plain --stats_fmt with_stats \
                        {compactors_dir}/after_correction.scores.top_effect_size_bin.tsv \
                        {lookup_table_dir}/after_correction.scores.top_effect_size_bin.{name}.compactors.tsv", file=f)

        cmd = f"{lookup_table} query_many --n_threads {max_cpus_to_use_in_auto_adjust} {tmp_dir}/lookup_table_{name}_input.txt {path}"
        logs_files.append(f"{logs_dir}/lookup_table.{name}.log")
        with open(logs_files[-1], "w") as stdoutfile:
            run_cmd(cmd, stdoutfile, stdoutfile)

        if clean_up:
            os.remove(f"{tmp_dir}/lookup_table_{name}_input.txt")

    memory_logger.report("M_lookup_table", get_max_maximum_resident_set_size_many_files(logs_files))
    memory_logger.report(f"S_{lookup_table_dir}", calc_dir_size(lookup_table_dir))
    store_current_time_disk_ram_logs()
    timer.catch("T_lookup_table")
    print("lookup_table done")

store_current_time_disk_ram_logs()

if clean_up:
    if EXPERIMENTAL_INPUT_GENERATORS:
        for path_to_input_without_generators_for_10x in paths_to_input_without_generators_for_10x:
            os.remove(path_to_input_without_generators_for_10x)
        os.remove(path_to_input_without_generators)
    os.rmdir(tmp_dir)
check_and_handle_error()

splash_postprocessing.Run()

print("SPLASH finished!")
print("Current time:", get_cur_time(), flush=True)
