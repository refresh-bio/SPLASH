#!/bin/bash

bin_dir=$1
cd $bin_dir
if test -f "kmc" && test -f "kmc_tools"; then
    exit 0
fi

wget https://github.com/refresh-bio/KMC/releases/download/v3.2.2/KMC3.2.2.linux.x64.tar.gz
tar -xvf KMC3.2.2.linux.x64.tar.gz
mv bin/kmc .
mv bin/kmc_tools .
rm -rf bin
rm -rf include
rm KMC3.2.2.linux.x64.tar.gz
