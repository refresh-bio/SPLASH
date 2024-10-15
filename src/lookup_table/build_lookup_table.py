#!/usr/bin/env python3
import argparse
import sys
from time import gmtime, strftime, localtime
import subprocess
import os
import shutil
from sys import platform
import uuid

def get_main_splash_script():
    names = ["splash", "splash.py"]
    paths = [os.path.dirname(__file__), ".."]
    for p in paths:
        for n in names:
            x = os.path.join(p, n)
            if os.path.exists(x):
                return x
    print("Error: cannot locate main splash script")
    sys.exit(1)

def get_splash_version():
    with open(get_main_splash_script()) as f:
        for line in f:
            line = line.strip()
            if "SPLASH_VERSION" in line:
                return line.split("=")[1].split("\"")[1]
    print("Error: cannot determine splash version")
    sys.exit(1)

SPLASH_VERSION = get_splash_version()

parser = argparse.ArgumentParser(
                    prog = "splash",
                    description = "Welcome to lookup table builder\nVersion: " + SPLASH_VERSION,
                    #epilog = 'Text at the bottom of help',
                    #formatter_class=SmartFormatter
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter
                    )

parser.add_argument("input_file", help="path to the file where input FASTA files are defined, the format is: per each line path to fasta file")

parser.add_argument("--transcriptomes", default="", help="path to additional (optional) file where transcriptome input FASTA files are defined, the format is: per each line path to fasta file")
#*.slt stands for splash lookup table
parser.add_argument("--outname", default="lookup.slt", type=str, help="prefix of output file names")
parser.add_argument("--kmer_len", default=15, type=int, help="k-mer length")
parser.add_argument("--bin_path", default=".", type=str, help="path to a directory where kmc, kmc_tools")
parser.add_argument("--n_threads", default=8, type=int, help="number of threads")
parser.add_argument("--poly_ACGT_len", default=0, type=int, help="all k-mers containing polyACGT of this length will be filtered out (0 means no filtering)")
#parser.add_argument("--variant", default="sbwt", type=str, choices=['sbwt'], help="k-mer k-mer index type (sbwt is current recommendation)")
parser.add_argument("--category_3_threshold", default=1, type=int, help="accept k-mer in category 3 if its present in a given file <=category_3_threshold times")
parser.add_argument("--precomputed_sbwt", default="", type=str, help="path to precomputed sbwt index (if set lookup_table will use it instead of building own sbwt - must be build for exactly the same set of input files!)")
parser.add_argument("--dont_clean_up", default=False, action='store_true', help="if set then intermediate files will not be removed")

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

transcriptomes = args.transcriptomes
outname = args.outname
kmer_len = args.kmer_len
bin_path = args.bin_path
n_threads = args.n_threads
poly_ACGT_len = args.poly_ACGT_len
variant = "sbwt" #variant = args.variant
category_3_threshold = args.category_3_threshold
input_file = args.input_file
precomputed_sbwt = args.precomputed_sbwt
dont_clean_up = args.dont_clean_up

def get_cur_time():    
    return strftime("%Y-%m-%d %H:%M:%S", localtime())

print("Welcome to  lookup table builder")
print("Version: ", SPLASH_VERSION)
print("Current time:", get_cur_time(), flush=True)

print("----------------  Configuration  ---------------")
for arg, value in vars(args).items():
    print("%s: %s" % (arg, value))
print("------------------------------------------------", flush=True)

def wrap_cmd_with_time(cmd):
    if platform == "darwin":
        if shutil.which("gtime") is not None:
            cmd = f"gtime -v {cmd}"
        else:
            cmd = f"/usr/bin/time {cmd}" # no -v for /usr/bin/time on mac os
    else:
        cmd = f"/usr/bin/time -v {cmd}"
    return cmd

def run_cmd(cmd):
    cmd = wrap_cmd_with_time(cmd)
    print(get_cur_time() + ": " + cmd, flush=True)
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()
    if p.returncode != 0:
        print(f"Error running command: {cmd}", flush=True)

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


kmc=get_prog_or_fail("kmc")
kmc_tools=get_prog_or_fail("kmc_tools")
lookup_table=get_prog_or_fail("lookup_table")

for executable_path in [lookup_table]:
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
        print(f"Error: cannot find version number for {executable_path}.")
        sys.exit(1)


tmp_dir="lookup-table-tmp-"+uuid.uuid4().hex

if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

kmc_dbs = []
fasta_paths = []
uniq_dbs = []

with open(input_file) as f:
    for line in f:
        line = line.strip()
        if not line: #skip empty lines
            continue
        fasta_paths.append((line, False)) # False means this is not transcriptome

if transcriptomes != "":
    with open(transcriptomes) as f:
        for line in f:
            line = line.strip()
            if not line: #skip empty lines
                continue
            fasta_paths.append((line, True)) # True means this is transcriptome


for (line, is_transcriptome) in fasta_paths:
    b = os.path.basename(line)
    cmd = f"{kmc} -k{kmer_len} -t{n_threads} -fm -cs65536 -ci1 -t8 -n64 {line} {tmp_dir}/{b} {tmp_dir}"
    kmc_dbs.append(f"{tmp_dir}/{b}")
    run_cmd(cmd)

union_kmc_db_path = os.path.join(tmp_dir, "UNION")

#build a database of all k-mers observed in the input (union all DBs)
with open(f"{tmp_dir}/complex.op", "w") as f:
    f.write("INPUT:\n")
    o = ""
    for i, a in enumerate(kmc_dbs):
        f.write(f"i{i}={a}\n")
        o += f"i{i}+"
    o = o[:-1] #remove last +
    f.write("OUTPUT:\n")
    f.write(f"{union_kmc_db_path}={o}\n")

cmd = f"{kmc_tools} -t{n_threads} complex {tmp_dir}/complex.op"
run_cmd(cmd)

for i, a in enumerate(kmc_dbs):
    uniq_out_op_file = f"{tmp_dir}/uniq_{i}.op"
    uniq_out_path = f"{a}_uniq"
    uniq_dbs.append(uniq_out_path)
    with open(f"{uniq_out_op_file}", "w") as f:
        f.write("INPUT:\n")
        f.write(f"B={a}\n")
        f.write(f"S={union_kmc_db_path}\n")
        f.write("OUTPUT:\n")
        f.write(f"{uniq_out_path}=B-(S~B)\n")
    cmd = f"{kmc_tools} -t{n_threads} complex {uniq_out_op_file}"
    run_cmd(cmd)

assert(len(kmc_dbs) == len(fasta_paths))
assert(len(fasta_paths) == len(uniq_dbs))

with open(f"{tmp_dir}/build_input.txt", "w") as f:
    for ((fasta_path, is_transcriptome), all_db, uniq_db) in zip(fasta_paths, kmc_dbs, uniq_dbs):
        f.write(f"{fasta_path};{all_db};{uniq_db};{is_transcriptome}\n")


_precomputed_sbwt_param = f"--precomputed_sbwt {precomputed_sbwt}" if precomputed_sbwt != "" else ""
cmd = f"{lookup_table} build \
    --n_threads {n_threads} \
    --poly_ACGT_len {poly_ACGT_len} \
    --variant {variant} \
    --category_3_threshold {category_3_threshold} \
    --tmp_dir {tmp_dir} \
    {_precomputed_sbwt_param} \
    {tmp_dir}/build_input.txt {outname}"
run_cmd(cmd)


#../../bin/lookup_table query 1.fa lookup o.txt

if not dont_clean_up:
    shutil.rmtree(tmp_dir)
