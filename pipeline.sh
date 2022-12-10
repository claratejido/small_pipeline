#!/bin/bash

# Create an empty dictionary to count 3-mer occurrences
declare -A occ

# Create data folder if it does not exist
data_folder="./data/"
mkdir -p "$data_folder"

# Read FASTA records from input file and count 3-mer occurrences
index=0
while IFS='' read -r tag || [[ -n "$tag" ]]; do
  # Read sequence
  read -r sequence_raw
  sequence="${sequence_raw//[[:space:]]/}"
  # Trim sequence if it is longer than 20 nucleotides
  if [[ "${#sequence}" -gt 20 ]]; then
    trimmed_sequence="${sequence::-20}"
  else
    trimmed_sequence=""
    echo "No sequence written as sequence is shorter than 20 nucleotides"
  fi
  # Count 3-mer occurrences in sequence
  for ((i=0; i<"${#sequence}"-2; i++)); do
    kmer="${sequence:i:3}"
    if [[ -z "${occ[$kmer]}" ]]; then
      occ[$kmer]=0
    fi
    occ[$kmer]=$((occ[$kmer]+1))
  done
  # Write FASTA record to file
  output_file_name="sequences.${index}.fasta"
  output_file_path="$data_folder$output_file_name"
  echo -e "$tag$trimmed_sequence" > "$output_file_path"
  index=$((index+1))
done < "$1"

# Print 3-mer occurrences
for kmer in "${!occ[@]}"; do
  echo "$kmer: ${occ[$kmer]}"
done

# Generate BWA index from reference genome
bwa index GCF_000146045.2_R64_genomic.fna

# Align FASTA files to reference genome and output SAM files
for fasta_file in "$data_folder"*; do
  bwa mem GCF_000146045.2_R64_genomic.fna "$fasta_file" > "${fasta_file%.*}.sam"
done

# Merge SAM files into a single file
merged_sam_file="$data_folder"merged.sam
for sam_file in "$data_folder"*.sam; do
  # Skip the merged SAM file
  if [[ "$sam_file" == "$merged_sam_file" ]]; then
    continue
  fi
  # Append the content of the SAM file to the merged SAM file
  cat "$sam_file" >> "$merged_sam_file"
done

# Sort merged SAM file by chromosome and position
sort -k 3,3 -k 4,4n "$merged_sam_file" > "${merged_sam_file%.*}_sorted.sam"

# Count the number of aligned bases in the SAM file
aligned_count=$(grep -v '^@' "${merged_sam_file%.*}_sorted.sam" | awk -F 'M' '{sum += $1} END {print sum}')
