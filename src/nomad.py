#!/usr/bin/env python3
import argparse
import os
import sys
from time import gmtime, strftime, localtime
import subprocess
import threading
import queue
import gzip
import shutil
import uuid

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

NOMAD_VERSION="2.0.0"

parser = argparse.ArgumentParser(
                    prog = "nomad",
                    description = "Welcome to NOMAD\nVersion: " + NOMAD_VERSION,
                    #epilog = 'Text at the bottom of help',
                    #formatter_class=SmartFormatter
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                    )

parser.add_argument("input_file", help="path to the file where input samples are defined, the format is: per each line {sample_name}<space>{path}, path is a fastq[.gz] file")

group_base_configuration = parser.add_argument_group('Base configuration')
group_base_configuration.add_argument("--outname_prefix", default="result", type=str, help="prefix of output file names")
group_base_configuration.add_argument("--anchor_len", default=27, type=int, help="anchor length")
group_base_configuration.add_argument("--gap_len", default="0", type=str, help="gap length, if 'auto' it will be inferred from the data, in the opposite case it must be an int")
group_base_configuration.add_argument("--target_len", default=27, type=int, help="target length")
group_base_configuration.add_argument("--anchor_list", default="", type=str, help="list of accepted anchors, this is path to plain text file with one anchor per line without any header")
group_base_configuration.add_argument("--pvals_correction_col_name", default="pval_opt", type=str, help="for which column correction should be applied")

group_filters_and_thresholds = parser.add_argument_group('Filters and thresholds')
group_filters_and_thresholds.add_argument("--poly_ACGT_len", default=8, type=int, help="filter out all anchors containing poly(ACGT) of length at least <poly_ACGT_len> (0 means no filtering)")
group_filters_and_thresholds.add_argument("--anchor_unique_targets_threshold", default=1, type=int, help="filter out all anchors for which the number of unique targets is <= anchor_unique_targets_threshold")
group_filters_and_thresholds.add_argument("--anchor_count_threshold", default=50, type=int, help="filter out all anchors for which the total count <= anchor_count_threshold")
group_filters_and_thresholds.add_argument("--anchor_samples_threshold", default=1, type=int, help="filter out all anchors for which the number of unique samples is <= anchor_samples_threshold")
group_filters_and_thresholds.add_argument("--anchor_sample_counts_threshold", default=5, type=int, help="filter out anchor from sample if its count in this sample is <= anchor_sample_counts_threshold")
group_filters_and_thresholds.add_argument("--n_most_freq_targets_for_stats", default=0, type=int, help="use at most n_most_freq_targets_for_stats for each contignency table, 0 means use all")
group_filters_and_thresholds.add_argument("--fdr_threshold", default=0.05, type=float, help="keep anchors having corrected p-val below this value")

group_additional_out = parser.add_argument_group('Additional output configuration')
group_additional_out.add_argument("--dump_Cjs", default=False, action='store_true', help="output Cjs")
group_additional_out.add_argument("--max_pval_opt_for_Cjs", default=0.10, type=float, help="dump only Cjs for anchors that have pval_opt <= max_pval_opt_for_Cjs")
group_additional_out.add_argument("--n_most_freq_targets", default=2, type=int, help="number of most frequent tragets printed per each anchor in stats mode")
group_additional_out.add_argument("--with_effect_size_cts", default=False, action='store_true', help="if set effect_size_cts will be computed")
group_additional_out.add_argument("--sample_name_to_id", default="sample_name_to_id.mapping.txt", type=str, help="file name with mapping sample name <-> sammpe id")
group_additional_out.add_argument("--dump_sample_anchor_target_count_txt", default=False, action='store_true', help="if set contignency tables will be generated in text format")
group_additional_out.add_argument("--dump_sample_anchor_target_count_binary", default=False, action='store_true', help="if set contignency tables will be generated in binary (satc) format, to convert to text format later satc_dump program may be used, it may take optionally maping from id to sample_name (--sample_names param)")

group_tuning_stats = parser.add_argument_group('Tuning statistics computation')
group_tuning_stats.add_argument("--opt_num_inits", default=10, type=int, help="the number of altMaximize runs")
group_tuning_stats.add_argument("--opt_num_iters", default=50, type=int, help="the number of iteration in altMaximize")
group_tuning_stats.add_argument("--num_rand_cf", default=50, type=int, help="the number of rand cf")
group_tuning_stats.add_argument("--opt_train_fraction", default=0.25, type=float, help="in calc_stats mode use this fraction to create train X from contingency table")

group_technical = parser.add_argument_group('Technical and performance-related')
group_technical.add_argument("--bin_path", default="./", type=str, help="path to a directory where satc, satc_dump, satc_merge, sig_anch, kmc, kmc_tools binaries are (if any not found there nomad will check if installed and use installed)")
group_technical.add_argument("--tmp_dir", default="", type=str, help="path to a directory where temporary files will be stored")
group_technical.add_argument("--n_threads_stage_1", default=4, type=int, help="number of threads for the first stage, too large value is not recomended because of intensive disk access here, but may be profitable if there is a lot of small size samples in the input")
group_technical.add_argument("--n_threads_stage_1_internal", default=8, type=int, help="number of threads per each stage 1 thread")
group_technical.add_argument("--n_threads_stage_2", default=32, type=int, help="number of threads for the second stage, high value is recommended if possible, single thread will process single bin")
group_technical.add_argument("--n_bins", default=128, type=int, help="the data will be split in a number of bins that will be merged later")
group_technical.add_argument("--kmc_use_RAM_only_mode", default=False, action='store_true', help="True here may increase performance but also RAM-usage")
group_technical.add_argument("--kmc_max_mem_GB", default=12, type=int, help="maximal amount of memory (in GB) KMC will try to not extend")
group_technical.add_argument("--dont_clean_up", default=False, action='store_true', help="if set then intermediate files will not be removed")
group_technical.add_argument("--logs_dir", default="logs", type=str, help="director where run logs of each thread will be stored")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

bin_path=args.bin_path
n_threads_stage_1 = args.n_threads_stage_1
n_threads_stage_1_internal = args.n_threads_stage_1_internal
n_threads_stage_2 = args.n_threads_stage_2
anchor_list = args.anchor_list
n_bins = args.n_bins
anchor_len = args.anchor_len
target_len = args.target_len
gap_len = args.gap_len
poly_ACGT_len=args.poly_ACGT_len
anchor_unique_targets_threshold = args.anchor_unique_targets_threshold
anchor_count_threshold = args.anchor_count_threshold
anchor_samples_threshold = args.anchor_samples_threshold
anchor_sample_counts_threshold = args.anchor_sample_counts_threshold
n_most_freq_targets = args.n_most_freq_targets
n_most_freq_targets_for_stats=args.n_most_freq_targets_for_stats
opt_num_inits=args.opt_num_inits
opt_num_iters=args.opt_num_iters
num_rand_cf=args.num_rand_cf
opt_train_fraction=args.opt_train_fraction
kmc_use_RAM_only_mode = args.kmc_use_RAM_only_mode
kmc_max_mem_GB = args.kmc_max_mem_GB
with_effect_size_cts = args.with_effect_size_cts
sample_name_to_id = args.sample_name_to_id
dump_sample_anchor_target_count_txt = args.dump_sample_anchor_target_count_txt
dump_sample_anchor_target_count_binary = args.dump_sample_anchor_target_count_binary
pvals_correction_col_name = args.pvals_correction_col_name
fdr_threshold = args.fdr_threshold
outname_prefix=args.outname_prefix
dump_Cjs=args.dump_Cjs
max_pval_opt_for_Cjs=args.max_pval_opt_for_Cjs
clean_up=not args.dont_clean_up
logs_dir=args.logs_dir
input_file=args.input_file
tmp_dir=args.tmp_dir

if tmp_dir == "":
    tmp_dir="nomad-tmp-"+uuid.uuid4().hex


def get_cur_time():    
    return strftime("%Y-%m-%d %H:%M:%S", localtime())

print("Welcome to NOMAD")
print("Version: ", NOMAD_VERSION)
print("Current time:", get_cur_time(), flush=True)


print("----------------  Configuration  ---------------")
for arg, value in vars(args).items():
    print("%s: %s" % (arg, value))
print("------------------------------------------------", flush=True)


def get_gap_len_from_the_data():
    def get_first_line_from_txt_file(path):
        try:
            with open(path) as f:
                return f.readline().strip()
        except:
            print(f"Error: cannot read first line from {path}")
            sys.exit(1)

    path = get_first_line_from_txt_file(input_file).split()[1]

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

def get_prog_or_fail(prog):
    p = os.path.join(bin_path, prog)
    if os.path.exists(p):
        return p

    #maybe it is installed
    if shutil.which(prog) is not None:
        return prog

    print(f"Error: {prog} does not exist in {bin_path} nor it is installed.")
    sys.exit(1)

satc=get_prog_or_fail("satc")
satc_dump=get_prog_or_fail("satc_dump")
satc_merge=get_prog_or_fail("satc_merge")
sig_anch=get_prog_or_fail("sig_anch")
kmc=get_prog_or_fail("kmc")
kmc_tools=get_prog_or_fail("kmc_tools")

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

def run_cmd(cmd, out, err):
    global was_error
    check_and_handle_error()
    cmd = f"/usr/bin/time -v {cmd}"
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


for executable_path in [satc, satc_dump, satc_merge, sig_anch]:
    p = subprocess.Popen(executable_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode("utf-8")
    err = err.decode("utf-8")
    out += err
    version_is_OK = False
    for line in out.split('\n'):
        if line.startswith("Version:"):
            ver = line.split()[-1]
            if ver != NOMAD_VERSION:
                print(f"Error: {executable_path} version is {ver} while this script is {NOMAD_VERSION}. Please use the same version.")
                sys.exit(1)
            else:
                version_is_OK=True
                break
    if not version_is_OK:
        print(f"Error: cannot find version number for {executable_path}.")
        sys.exit(1)

if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

with open(f"{logs_dir}/nomad-cmd.log", "w") as f:
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

os.makedirs(tmp_dir)

###############################################################################
# STAGE 1
# Construct (sample, anchor, target, count) records per each sample
# per each sample there is {n_bins} output files
###############################################################################

print("Starting stage 1")
print("Current time:", get_cur_time(), flush=True)

def stage_1_task(id, input, out, err):
    fname = input[0]
    sample_name = input[1]
    
    kmc_dir_tmp_name = f"{tmp_dir}/kmc_tmp_{sample_name}"

    if not os.path.exists(kmc_dir_tmp_name):
        os.makedirs(kmc_dir_tmp_name)

    ram_only_param = ""
    if kmc_use_RAM_only_mode:
        ram_only_param="-r"

    cmd = f"{kmc} -t{n_threads_stage_1_internal} -b -ci1 -cs65535 -k{anchor_len + gap_len + target_len} {ram_only_param} -m{kmc_max_mem_GB} {fname} {tmp_dir}/{sample_name} {kmc_dir_tmp_name}"
    run_cmd(cmd, out, err)

    if clean_up:
        os.rmdir(kmc_dir_tmp_name)

    cmd = f"{kmc_tools} -t{n_threads_stage_1_internal} transform {tmp_dir}/{sample_name} sort {tmp_dir}/{sample_name}.sorted"
    run_cmd(cmd, out, err)

    if clean_up:
        remove_kmc_output(f"{tmp_dir}/{sample_name}")
    
    cmd = f"{satc} --anchor_len {anchor_len} --target_len {target_len} --n_bins {n_bins} --anchor_sample_counts_threshold {anchor_sample_counts_threshold} --poly_ACGT_len {poly_ACGT_len} {tmp_dir}/{sample_name} {tmp_dir}/{sample_name}.sorted {id}"
    run_cmd(cmd, out, err)
    
    if clean_up:
        remove_kmc_output(f"{tmp_dir}/{sample_name}.sorted")

class action_at_function_exit:
    def __init__(self, action):
        self.action = action
    def __del__(self):
        self.action()

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
# this stage may be (and is intended to be) used to compute pvals
# the current state of this may be enabled with --run_mode parameter
# in the case when --run_mode is set to just_merge_and_dump
# the format of the dump may be configured with --format parameter
###############################################################################

print("Stage 1 done")
print("Starting stage 2")
print("Current time:", get_cur_time(), flush=True)

def stage_2_task(bin_id, out, err):
    satc_merge_inputs = []

    for input in inputs:
        sample_name = input[1]

        name = f"{tmp_dir}/{sample_name}.{bin_id}.bin"

        satc_merge_inputs.append(name)
        
    in_file_name = f"{tmp_dir}/files.bin{bin_id}.lst"
    with open(in_file_name, "w") as f:
        for x in satc_merge_inputs:
            f.write(f"{x}\n")
    
    _with_effect_size_cts_param = "--with_effect_size_cts" if with_effect_size_cts else ""

    Cjs_dir = f"{outname_prefix}_Cjs"
    if dump_Cjs and not os.path.exists(Cjs_dir):
        os.makedirs(Cjs_dir)

    _cjs_out_param = f"--cjs_out {Cjs_dir}/bin{bin_id}.cjs" if dump_Cjs else ""

    cmd=f"{satc_merge} \
    {_with_effect_size_cts_param} \
    {_cjs_out_param} \
    --max_pval_opt_for_Cjs {max_pval_opt_for_Cjs} \
    --anchor_count_threshold {anchor_count_threshold} \
    --anchor_samples_threshold {anchor_samples_threshold} \
    --anchor_unique_targets_threshold {anchor_unique_targets_threshold} \
    --n_most_freq_targets {n_most_freq_targets} \
    --n_most_freq_targets_for_stats {n_most_freq_targets_for_stats} \
    --opt_num_inits {opt_num_inits} \
    --opt_num_iters {opt_num_iters} \
    --num_rand_cf {num_rand_cf} \
    --opt_train_fraction {opt_train_fraction} \
    --run_mode calc_stats \
    --sample_names {sample_name_to_id} \
    {anchor_list_param} \
    {tmp_dir}/{outname_prefix}.bin{bin_id}.stats.tsv {in_file_name}"
    run_cmd(cmd, out, err)

    if dump_sample_anchor_target_count_txt:
        dump_dir = f"{outname_prefix}_dumps"
        if not os.path.exists(dump_dir):
            os.makedirs(dump_dir)

        cmd=f"{satc_merge} \
        --anchor_count_threshold {anchor_count_threshold} \
        --anchor_samples_threshold {anchor_samples_threshold} \
        --anchor_unique_targets_threshold {anchor_unique_targets_threshold} \
        --run_mode just_merge_and_dump \
        --sample_names {sample_name_to_id} \
        --format satc \
        {anchor_list_param} \
        {dump_dir}/bin{bin_id}.satc.dump {in_file_name}"
        run_cmd(cmd, out, err)

    if dump_sample_anchor_target_count_binary:
        satc_dir = f"{outname_prefix}_satc"
        if not os.path.exists(satc_dir):
            os.makedirs(satc_dir)

        cmd=f"{satc_merge} \
        --anchor_count_threshold {anchor_count_threshold} \
        --anchor_samples_threshold {anchor_samples_threshold} \
        --anchor_unique_targets_threshold {anchor_unique_targets_threshold} \
        --run_mode just_merge \
        {anchor_list_param} \
        {satc_dir}/bin{bin_id}.satc {in_file_name}"
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

print("Stage 2 done")
print("Starting stage 3")
print("Current time:", get_cur_time(), flush=True)

files_for_correction = open(f"{tmp_dir}/files_for_correction.list.txt", "w")
for bin_id in range(n_bins):
    files_for_correction.write(f"{tmp_dir}/{outname_prefix}.bin{bin_id}.stats.tsv" + "\n")
files_for_correction.close()

cmd=f"{sig_anch} \
--samplesheet {input_file} \
--infile_bins {tmp_dir}/files_for_correction.list.txt \
--outfile_scores {outname_prefix}.after_correction.scores.tsv \
--outfile_all_anchors_pvals {outname_prefix}.after_correction.all_anchors.tsv \
--fdr_threshold {fdr_threshold} \
--col_name {pvals_correction_col_name}"

stdoutfile = open(f"{logs_dir}/correction.log", "w")
run_cmd(cmd, stdoutfile, stdoutfile)
stdoutfile.close()

if clean_up:
    for bin_id in range(n_bins):
        os.remove(f"{tmp_dir}/{outname_prefix}.bin{bin_id}.stats.tsv")
    os.remove(f"{tmp_dir}/files_for_correction.list.txt")
    os.rmdir(tmp_dir)
check_and_handle_error()
print("Stage 3 done")
print("NOMAD finished!")
print("Current time:", get_cur_time(), flush=True)
