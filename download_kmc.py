#!/usr/bin/env python3
import os
import sys
import subprocess

bin_dir=sys.argv[1]

os.chdir(bin_dir)

if os.path.exists("kmc") and os.path.exists("kmc_tools"):
    sys.exit(0)

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

system = get_os()
hardware = get_hardware()

tar_name=f"KMC3.2.2.{system}.{hardware}.tar.gz"
URL=f"https://github.com/refresh-bio/KMC/releases/download/v3.2.2/{tar_name}"

run_cmd(f"wget {URL}")
run_cmd(f"tar -xvf {tar_name}")
run_cmd("mv bin/kmc .")
run_cmd("mv bin/kmc_tools .")
run_cmd("rm -rf bin")
run_cmd("rm -rf include")
run_cmd(f"rm {tar_name}")
