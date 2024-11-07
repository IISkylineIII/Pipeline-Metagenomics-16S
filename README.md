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

    
    # Group FASTQ files in pairs # 
    
import os
from collections import defaultdict


paired_files = defaultdict(list)


for file in os.listdir():
    if file.endswith("_R1.fastq") or file.endswith("_R2.fastq"):

        sample_name = file.rsplit("_R", 1)[0]
        paired_files[sample_name].append(file)


for sample, files in paired_files.items():
    if len(files) == 2:
        print(f"Sample: {sample} -> Files: {files[0]}, {files[1]}")
    else:
        print(f"Sample: {sample} -> Incomplete pair: {files}")

   # Pairing files #

   
fastq_files = [
    "F3D0_R1.fastq", "F3D0_R2.fastq",
    "F3D1_R1.fastq", "F3D1_R2.fastq",
    "F3D2_R1.fastq", "F3D2_R2.fastq",
    "F3D3_R1.fastq", "F3D3_R2.fastq",
    "F3D4_R1.fastq", "F3D4_R2.fastq",
    "F3D5_R1.fastq", "F3D5_R2.fastq",
    "F3D6_R1.fastq", "F3D6_R2.fastq",
    "F3D7_R1.fastq", "F3D7_R2.fastq",
    "F3D8_R1.fastq", "F3D8_R2.fastq",
    "F3D9_R1.fastq", "F3D9_R2.fastq",
    "F3D141_R1.fastq", "F3D141_R2.fastq",
    "F3D142_R1.fastq", "F3D142_R2.fastq",
    "F3D143_R1.fastq", "F3D143_R2.fastq",
    "F3D144_R1.fastq", "F3D144_R2.fastq",
    "F3D145_R1.fastq", "F3D145_R2.fastq",
    "F3D146_R1.fastq", "F3D146_R2.fastq",
    "F3D147_R1.fastq", "F3D147_R2.fastq",
    "F3D148_R1.fastq", "F3D148_R2.fastq",
    "F3D149_R1.fastq", "F3D149_R2.fastq",
    "F3D150_R1.fastq", "F3D150_R2.fastq",
    "Mock_R1.fastq", "Mock_R2.fastq"
]


sample_pairs = {}

for file in fastq_files:
    sample_name = file.split('_')[0]
    if sample_name not in sample_pairs:
        sample_pairs[sample_name] = []
    sample_pairs[sample_name].append(file)


renamed_pairs = {f"{name}": files for name, files in sample_pairs.items()}

for sample, files in renamed_pairs.items():
    print(f"Sample: {sample} -> Files: {', '.join(files)}")


with open('paired_samples.txt', 'w') as output_file:
    for sample, files in renamed_pairs.items():
        output_file.write(f"{sample} -> {', '.join(files)}\n")

# Quality assessment of the read #

fastq_files = [
    "F3D0_R1.fastq", "F3D0_R2.fastq",
    "F3D1_R1.fastq", "F3D1_R2.fastq",
    "F3D2_R1.fastq", "F3D2_R2.fastq",
    "F3D3_R1.fastq", "F3D3_R2.fastq",
    "F3D5_R1.fastq", "F3D5_R2.fastq",
    "F3D6_R1.fastq", "F3D6_R2.fastq",
    "F3D7_R1.fastq", "F3D7_R2.fastq",
    "F3D8_R1.fastq", "F3D8_R2.fastq",
    "F3D9_R1.fastq", "F3D9_R2.fastq",
    "F3D141_R1.fastq", "F3D141_R2.fastq",
    "F3D142_R1.fastq", "F3D142_R2.fastq",
    "F3D143_R1.fastq", "F3D143_R2.fastq",
    "F3D144_R1.fastq", "F3D144_R2.fastq",
    "F3D145_R1.fastq", "F3D145_R2.fastq",
    "F3D146_R1.fastq", "F3D146_R2.fastq",
    "F3D147_R1.fastq", "F3D147_R2.fastq",
    "F3D148_R1.fastq", "F3D148_R2.fastq",
    "F3D149_R1.fastq", "F3D149_R2.fastq",
    "F3D150_R1.fastq", "F3D150_R2.fastq",
    "Mock_R1.fastq", "Mock_R2.fastq"
]

from Bio import SeqIO
import os

def evaluate_quality(fastq_files):
    for fastq_file in fastq_files:
        if not os.path.isfile(fastq_file):
            print(f"File not found: {fastq_file}")
            continue
        print(f"Quality statistics for {fastq_file}:")
        total_bases = 0
        total_quality = 0
        read_count = 0

        for record in SeqIO.parse(fastq_file, "fastq"):
            total_bases += len(record.seq)
            total_quality += sum(record.letter_annotations["phred_quality"])
            read_count += 1

        average_quality = total_quality / total_bases if total_bases > 0 else 0
        print(f"Total reads: {read_count}, Average Quality: {average_quality:.2f}")


evaluate_quality(fastq_files)

# Combine forward and reverse reads #

  from Bio import SeqIO
import glob


groups = [
    "F3D0", "F3D1", "F3D2", "F3D3", "F3D5", "F3D6", "F3D7", "F3D8", "F3D9",
    "F3D141", "F3D142", "F3D143", "F3D144", "F3D145", "F3D146", "F3D147",
    "F3D148", "F3D149", "F3D150", "Mock"
]

def calculate_average_quality(fastq_file):
    total_quality = 0
    total_reads = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        total_quality += sum(record.letter_annotations["phred_quality"])
        total_reads += len(record.letter_annotations["phred_quality"])
    return total_quality / total_reads if total_reads > 0 else 0

quality_stats = {}
for group in groups:
    r1_file = f"{group}_R1.fastq"
    r2_file = f"{group}_R2.fastq"


    avg_quality_r1 = calculate_average_quality(r1_file)
    avg_quality_r2 = calculate_average_quality(r2_file)


    quality_stats[group] = {
        "R1_total_reads": sum(1 for _ in SeqIO.parse(r1_file, "fastq")),
        "R1_average_quality": avg_quality_r1,
        "R2_total_reads": sum(1 for _ in SeqIO.parse(r2_file, "fastq")),
        "R2_average_quality": avg_quality_r2
    }


for group, stats in quality_stats.items():
    print(f"Quality statistics for {group}:")
    print(f"  {group}_R1 - Total reads: {stats['R1_total_reads']}, Average Quality: {stats['R1_average_quality']:.2f}")
    print(f"  {group}_R2 - Total reads: {stats['R2_total_reads']}, Average Quality: {stats['R2_average_quality']:.2f}\n")

    # Unifying FASTQ file into all_samples_combined.fasta #

    from Bio import SeqIO


def merge_reads_and_convert_to_fasta(sample_names):
    combined_fasta_records = []

    for sample in sample_names:
        forward_file = f"{sample}_R1.fastq"
        reverse_file = f"{sample}_R2.fastq"

        forward_reads = SeqIO.parse(forward_file, "fastq")
        reverse_reads = SeqIO.parse(reverse_file, "fastq")

       a)
        for f_read, r_read in zip(forward_reads, reverse_reads):
            ado
            combined_seq = f_read.seq + r_read.seq


            combined_record = SeqIO.SeqRecord(
                seq=combined_seq,
                id=f"{f_read.id}|{r_read.id}",  s
                description=""
            )
            combined_fasta_records.append(combined_record)


    with open("all_samples_combined.fasta", "w") as fasta_file:
        SeqIO.write(combined_fasta_records, fasta_file, "fasta")


def generate_quality_statistics(sample_names):
    quality_stats = {}

    for sample in sample_names:
        forward_file = f"{sample}_R1.fastq"
        reverse_file = f"{sample}_R2.fastq"


        forward_reads = list(SeqIO.parse(forward_file, "fastq"))
        reverse_reads = list(SeqIO.parse(reverse_file, "fastq"))

       ia
        total_reads_fwd = len(forward_reads)
        total_reads_rev = len(reverse_reads)


        avg_quality_fwd = sum(sum(f_read.letter_annotations["phred_quality"]) for f_read in forward_reads) / total_reads_fwd
        avg_quality_rev = sum(sum(r_read.letter_annotations["phred_quality"]) for r_read in reverse_reads) / total_reads_rev

       s
        quality_stats[sample] = {
            "R1": (total_reads_fwd, avg_quality_fwd),
            "R2": (total_reads_rev, avg_quality_rev),
        }


    with open("quality_statistics.txt", "w") as stats_file:
        for sample, stats in quality_stats.items():
            stats_file.write(f"Quality statistics for {sample}:\n")
            for read_type, (total_reads, avg_quality) in stats.items():
                stats_file.write(f"  {sample}_{read_type} - Total reads: {total_reads}, Average Quality: {avg_quality:.2f}\n")
            stats_file.write("\n")


sample_names = [
    "F3D0", "F3D1", "F3D2", "F3D3", "F3D5",
    "F3D6", "F3D7", "F3D8", "F3D9", "F3D141",
    "F3D142", "F3D143", "F3D144", "F3D145",
    "F3D146", "F3D147", "F3D148", "F3D149",
    "F3D150", "Mock"
]

merge_reads_and_convert_to_fasta(sample_names)
generate_quality_statistics(sample_names)

print("Processamento concluído. Verifique 'all_samples_combined.fasta' e 'quality_statistics.txt'.")

# Combined files #


      with open('all_samples_combined.fasta', 'r') as fasta_file:
    fasta_content = fasta_file.read()
    print(fasta_content[:500])

# Statistical quality #

with open('quality_statistics.txt', 'r') as stats_file:
    stats_content = stats_file.read()
    print(stats_content)

 #   Creating contigs from forward and reverse reads with overlap # 

 A more sophisticated approach involves performing overlap to generate a consensus sequence, based on the Phred quality score for each position. If there is a discrepancy between the bases, the base with the higher quality is selected


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

def make_contigs(forward_path, reverse_path, output_fasta, output_groups, sample_id="Sample"):
    contigs = []
    group_info = []

    forward_reads = SeqIO.parse(forward_path, "fastq")
    reverse_reads = SeqIO.parse(reverse_path, "fastq")

    for forward, reverse in zip(forward_reads, reverse_reads):
        reverse_seq = reverse.seq.reverse_complement()
        forward_seq = forward.seq

        min_length = min(len(forward_seq), len(reverse_seq))
        consensus_seq = []
        consensus_qual = []

        for i in range(min_length):
            if forward_seq[i] == reverse_seq[i]:
                consensus_seq.append(forward_seq[i])
                consensus_qual.append(max(forward.letter_annotations["phred_quality"][i],
                                          reverse.letter_annotations["phred_quality"][i]))
            else:
                if forward.letter_annotations["phred_quality"][i] > reverse.letter_annotations["phred_quality"][i]:
                    consensus_seq.append(forward_seq[i])
                    consensus_qual.append(forward.letter_annotations["phred_quality"][i])
                else:
                    consensus_seq.append(reverse_seq[i])
                    consensus_qual.append(reverse.letter_annotations["phred_quality"][i])

        consensus_seq += forward_seq[min_length:] + reverse_seq[min_length:]
        consensus_qual += forward.letter_annotations["phred_quality"][min_length:] + reverse.letter_annotations["phred_quality"][min_length:]

        contig_seq = Seq("".join(consensus_seq))
        contig_record = SeqRecord(contig_seq, id=forward.id, description="consensus contig")
        contig_record.letter_annotations["phred_quality"] = consensus_qual
        contigs.append(contig_record)

        group_info.append(f"{forward.id}\t{sample_id}")

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(contigs, output_handle, "fasta")

    with open(output_groups, "w") as group_handle:
        for line in group_info:
            group_handle.write(line + "\n")

    print(f"Total de contigs gerados: {len(contigs)}")


sample_names = [
    "F3D0", "F3D1", "F3D2", "F3D3", "F3D5",
    "F3D6", "F3D7", "F3D8", "F3D9", "F3D141",
    "F3D142", "F3D143", "F3D144", "F3D145",
    "F3D146", "F3D147", "F3D148", "F3D149",
    "F3D150", "Mock"
]


base_path = "./"

for sample in sample_names:
    forward_reads_path = f"{base_path}{sample}_R1.fastq"  # Forward reads
    reverse_reads_path = f"{base_path}{sample}_R2.fastq"  # Reverse reads
    output_contigs_fasta = f"{base_path}{sample}_contigs.fasta"
    output_groups_txt = f"{base_path}{sample}_groups.txt"

    # Call the function with the actual paths
    make_contigs(forward_reads_path, reverse_reads_path, output_contigs_fasta, output_groups_txt, sample_id=sample)

    
 
# Combining the 20 x 2 contigs into a single FASTA file # 

from Bio import SeqIO

def combine_fasta_files(sample_names, output_combined_fasta, base_path="./"):
    combined_records = []

    for sample in sample_names:
        input_fasta_path = f"{base_path}{sample}_contigs.fasta"


        with open(input_fasta_path, "r") as input_handle:
            records = SeqIO.parse(input_handle, "fasta")
            combined_records.extend(records)

    with open(output_combined_fasta, "w") as output_handle:
        SeqIO.write(combined_records, output_handle, "fasta")

    print(f"Total contigs combined into {output_combined_fasta}: {len(combined_records)}")


output_combined_fasta = f"{base_path}combined_contigs.fasta"


combine_fasta_files(sample_names, output_combined_fasta, base_path)

 


