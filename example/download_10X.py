#!/usr/bin/env python3
import requests
import os
import sys

url_base='https://raw.githubusercontent.com/kaitlinchaung/nomad-kmc3/main/test_data/'
files = [
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L001_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L002_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L003_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L004_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L001_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L002_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L003_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L004_R2_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L001_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L002_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L003_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_1_S12_L004_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L001_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L002_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L003_R1_001.fastq.gz',
    'subset_TSP2_Muscle_diaphragm_10X_1_2_S13_L004_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L001_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L002_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L003_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L004_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L001_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L002_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L003_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L004_R2_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L001_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L002_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L003_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_1_S10_L004_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L001_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L002_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L003_R1_001.fastq.gz',
    'subset_TSP2_Muscle_rectusabdominus_10X_1_2_S11_L004_R1_001.fastq.gz'
]

for file in files:
    if not os.path.exists(file):    
        url = url_base + file
        print(f"downloading {url}")
        r = requests.get(url, allow_redirects=True)
        open(file, 'wb').write(r.content)
