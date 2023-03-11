#!/usr/bin/env python3
import requests
import os
import sys

url_base='https://raw.githubusercontent.com/salzman-lab/nomad/main/test_data/'
files = [
    'SRX14565338_SRR18431620_1.fastq.gz',
    'SRX14565342_SRR18431616_1.fastq.gz',
    'SRX14565346_SRR18431612_1.fastq.gz',
    'SRX14565361_SRR18431597_1.fastq.gz'    
]

for file in files:
    if not os.path.exists(file):    
        url = url_base + file
        print(f"downloading {url}")
        r = requests.get(url, allow_redirects=True)
        open(file, 'wb').write(r.content)
