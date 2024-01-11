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

def run_cmd(cmd):    
    p = subprocess.Popen(cmd, shell=True)
    p.communicate()

def run_cmd_get_stdout(cmd):
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    return p.stdout.decode('utf-8')

if __name__ == "__main__":
    run_cmd("git submodule init")
    run_cmd("git submodule update")
    run_cmd("make clean")
    run_cmd("make -j32 release")

    run_cmd("mkdir -p bin/example")

    run_cmd("cp example/download.py bin/example/download.py")

    run_cmd("cp example/input.txt bin/example")


    with open("bin/example/run-example.sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write("./download.py\n")
        f.write("../splash --bin_path .. input.txt\n")

    os.chmod("bin/example/run-example.sh", os.stat("bin/example/run-example.sh").st_mode | stat.S_IEXEC)


    ver = get_ver("bin/splash")

    run_cmd(f"cd bin; tar -c * | pigz > ../splash-{ver}.linux.x64.tar.gz; cd ..;")
    run_cmd("rm -rf bin")
