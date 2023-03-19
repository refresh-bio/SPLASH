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

def get_ver(nomad_path):
    with open(nomad_path) as f:
        for line in f.readlines():
            line = line.strip()
            if "NOMAD_VERSION" in line:
                return line.split("=")[-1].strip().split("\"")[1]
    print("Error: cannot read NOMAD_VERSION")
    sys.exit(1)

def run_cmd(cmd):    
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

run_cmd("make clean")
run_cmd("make -j32 release")

run_cmd("mkdir -p bin/example")

run_cmd("cp example/download.py bin/example/download.py")

run_cmd("cp example/input.txt bin/example")


with open("bin/example/run-example.sh", "w") as f:
    f.write("#!/bin/bash\n")
    f.write("./download.py\n")
    f.write("../nomad --bin_path .. input.txt\n")

os.chmod("bin/example/run-example.sh", os.stat("bin/example/run-example.sh").st_mode | stat.S_IEXEC)


ver = get_ver("bin/nomad")

run_cmd(f"cd bin; tar -c * | pigz > ../nomad-{ver}.linux.x64.tar.gz; cd ..;")
run_cmd("rm -rf bin")
