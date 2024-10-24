#!/usr/bin/env python3
import subprocess
import fileinput
import sys
import os
import stat

def replace_in_file(file_path, search_text, new_text):
    with fileinput.input(file_path, inplace=True) as file:
        for line in file:
            new_line = line.replace(search_text, new_text)
            print(new_line, end='')

def get_ver(splash_path):
    with open(splash_path) as f:
        for line in f.readlines():
            line = line.strip()
            if "SPLASH_VERSION" in line:
                return line.split("=")[-1].strip().split("\"")[1]
    print("Error: cannot read SPLASH_VERSION")
    sys.exit(1)

def get_os():
    if os.name == 'nt':
        return 'windows'
    elif os.name == 'posix':
        if os.uname()[0] == 'Linux':
            return 'linux'
        elif os.uname()[0] == 'Darwin':
            return 'mac'
        else:
            print("Error: unknown os", os.uname()[0])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def get_hardware():
    if os.name == 'nt':
        return 'x64' # TODO: do a real check and support ARM also...
    elif os.name == 'posix':
        if os.uname()[4] == 'x86_64':
            return 'x64'
        elif os.uname()[4] == 'aarch64' or os.uname()[4] == 'arm64':
            return 'arm64'
        else:
            print("Error: unknown hardware", os.uname()[4])
            sys.exit(1)
    else:
        print("Error: unknown os.name", os.name)
        sys.exit(1)

def run_cmd(cmd):
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

def run_cmd_get_stdout(cmd):
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return p.stdout.decode('utf-8')

if __name__ == "__main__":
    system = get_os()
    hardware = get_hardware()

    run_cmd("git submodule update --init --recursive")
    run_cmd("make clean")
    run_cmd("make -j release")

    run_cmd("mkdir -p bin/example")

    run_cmd("cp example/download.py bin/example/download.py")
    run_cmd("cp example/download_10X.py bin/example/download_10X.py")
    run_cmd("cp example/test_data_cell_barcode_samplesheet.csv bin/example/test_data_cell_barcode_samplesheet.csv")
    run_cmd("cp example/test_non_10X_Cj_samplesheet.csv bin/example/test_non_10X_Cj_samplesheet.csv")

    run_cmd("cp example/input.txt bin/example")

    run_cmd("cp example/input-10X.txt bin/example")
    run_cmd("cp example/S10.txt bin/example")
    run_cmd("cp example/S11.txt bin/example")
    run_cmd("cp example/S12.txt bin/example")
    run_cmd("cp example/S13.txt bin/example")

    with open("bin/example/run-example.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write("./download.py\n")
        f.write("../splash --bin_path .. input.txt\n")

    os.chmod("bin/example/run-example.sh", os.stat("bin/example/run-example.sh").st_mode | stat.S_IEXEC)


    ver = get_ver("bin/splash")

    run_cmd(f"cd bin; tar -c * | pigz > ../splash-{ver}.{system}.{hardware}.tar.gz; cd ..;")
    run_cmd("rm -rf bin")
