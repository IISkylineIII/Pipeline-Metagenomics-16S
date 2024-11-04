# Pipeline-Metagenomics-16S # 

# IMPORTING DATA FOR ANALYSIS #

import requests

urls = [
    "https://zenodo.org/record/800651/files/F3D0_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D0_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D141_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D141_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D142_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D142_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D143_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D143_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D144_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D144_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D145_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D145_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D146_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D146_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D147_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D147_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D148_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D148_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D149_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D149_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D150_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D150_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D1_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D1_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D2_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D2_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D3_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D3_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D5_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D5_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D6_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D6_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D7_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D7_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D8_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D8_R2.fastq",
    "https://zenodo.org/record/800651/files/F3D9_R1.fastq",
    "https://zenodo.org/record/800651/files/F3D9_R2.fastq",
    "https://zenodo.org/record/800651/files/Mock_R1.fastq",
    "https://zenodo.org/record/800651/files/Mock_R2.fastq",
]

for url in urls:
    file_name = url.split("/")[-1]
    r = requests.get(url)
    with open(file_name, 'wb') as f:
        f.write(r.content)
    print(f"{file_name}


   # LOADING REFERENCE DATA # 
   
import requests


reference_urls = [
    "https://zenodo.org/record/800651/files/HMP_MOCK.v35.fasta",
    "https://zenodo.org/record/800651/files/silva.v4.fasta",
    "https://zenodo.org/record/800651/files/trainset9_032012.pds.fasta",
    "https://zenodo.org/record/800651/files/trainset9_032012.pds.tax",
    "https://zenodo.org/record/800651/files/mouse.dpw.metadata",
]

for url in reference_urls:
    file_name = url.split("/")[-1]
    r = requests.get(url)
    with open(file_name, 'wb') as f:
        f.write(r.content)
    print(f"{file_name} 
   
