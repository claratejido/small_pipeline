import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Open file
input_file = open(sys.argv[1], "r")


# Declare empty dictionary
occ = {}
# Read FASTA records
index = 0
sequence_part_files = []

# Create data folder if it does not exist
data_folder = Path("./data/")
data_folder.mkdir(exist_ok=True)

while True:
    tag = input_file.readline()
    if not tag: # Read TAG or end-of-file
        break
    # Read sequence
    sequence_raw = input_file.readline()
    sequence = sequence_raw.rstrip('\n')
    if len(sequence) >20:
        trimmed_sequence= sequence[:-20]
    else:
        trimmed_sequence = ""
        print("No sequence written as sequence is shorter than 20 nucleotides")
    # Traverse sequence (3-mers)
    for i in range(len(sequence)-2):
        kmer = sequence[i:i+3]
        if kmer not in occ:
            occ[kmer] = 0
        occ[kmer] += 1
    output_file_name= f"sequences.{index}.fasta"
    output_file_path = "./data/" + output_file_name
    output_file = open(output_file_path, "wt")
    output_file.write("%s%s" % (tag,trimmed_sequence))
    output_file.close()
    sequence_part_files.append(output_file_name)
    index += 1


# Close files
input_file.close()


print(occ)

kmer_freq= occ.values()
kmer_labels=occ.keys()

freq =pd.Series(kmer_freq)
plt.figure()
ax = freq.plot(kind="bar")
ax.set_title("Kmer frequency")
ax.set_ylabel("Counts")
ax.set_xticklabels(kmer_labels)
plt.show()
histogram_filename = "kmer_frequencies.png"
plt.savefig(histogram_filename)
print(f"Kmer histogram saved in file: {histogram_filename}")

#Reference genome:
#https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# Generate BWA index from reference genome
os.system(f"bwa index GCF_000146045.2_R64_genomic.fna")

# Creating a SAM file alignment/fasta file and store them in data folder
sam_files = []
for i,fasta in enumerate(sequence_part_files):
    command = "bwa mem GCF_000146045.2_R64_genomic.fna"
    fasta_part = "./data/" + fasta
    sam_file_name = f"./data/aln-se{i}.sam"
    os.system(f"{command} {fasta_part} > {sam_file_name}")
    sam_files.append(sam_file_name)


# Merging sam files in a unique sam
merged_sam_file =  open("merged.sam", "wt")

for sam in sam_files:
    sam_file = open(sam, "rt")
    for line in sam_file:
        if line[0] != "@":
            merged_sam_file.write(line)

merged_sam_file.close()

merged_sam_file =  open("merged.sam", "rt")
sorted_merged_sam_file =  open("sorted_merged.sam", "wt")
sam_lines = []

for line in merged_sam_file:
    split_line = line.split("\t")
    sam_lines.append(split_line)


merged_sam_file.close()

sam_lines.sort(key = lambda x: (x[2], x[3]))
#sort by third field (chr) and fourth field (position)

aligned_count=0

for list in sam_lines:
    cigar = list[5]
    match = re.search(r'([0-9]+)M',cigar) # look for the number before the "M" in the cigar (bases aligned: match)
    if match:
        aligned_count += int(match.group(1)) #sum all the numbers before the M
    if cigar !=  "*":
        print(cigar)
    string ="\t".join(list)
    sorted_merged_sam_file.write(string)

sorted_merged_sam_file.close()

print (f"{aligned_count} bases aligned")

#sort the matches from the CIGAR (M), we are not interested in the soft clipped sequences (S) because they are the ones that do not align






