#!/usr/bin/env python3
import subprocess
import fileinput
import sys

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

run_cmd("make clean")
run_cmd("make -j32 release")

run_cmd("mkdir -p bin/example")
run_cmd("cp example/run.py bin/example")
run_cmd("cp example/splash.py bin")

run_cmd("mkdir -p bin/example/data")

run_cmd("cp example/data/download.py bin/example/data/download.py")

run_cmd("cp example/input.txt bin/example")

replace_in_file("bin/example/run.py", "../bin", "..")
replace_in_file("bin/example/run.py", "./splash.py", "../splash.py")

ver = get_ver("bin/splash.py")

run_cmd(f"cd bin; tar -c * | pigz > ../splash.{ver}.linux.tar.gz; cd ..;")
