#!/bin/bash

# Create an empty dictionary to count 3-mer occurrences
declare -A occ

# Create data folder if it does not exist
data_folder="./data/"
mkdir -p "$data_folder"
# Clear data folder, in case there is old data
rm $data_folder*

filename=$1
extension="${filename##*.}"
echo $extension

if [[ "$extension" == "fastq" ]]; then
  IS_FASTQ=1
else
  IS_FASTQ=0
fi
echo "Is fastq is : $IS_FASTQ"

# Read FASTA records from input file and count 3-mer occurrences
index=0
while IFS='' read -r tag || [[ -n "$tag" ]]; do
  # Read sequence
  read -r sequence_raw
  sequence="${sequence_raw//[[:space:]]/}"

  if [ $IS_FASTQ ]; then
    read -r fastq_separator
    read -r quality
  fi


  # Trim sequence if it is longer than 20 nucleotides
  if [[ "${#sequence}" -gt 20 ]]; then
    trimmed_sequence="${sequence::-20}"
    if [ $IS_FASTQ ]; then
      trimmed_quality="${quality::-20}"
    fi
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
  # Write record to file
  output_file_name="sequences.${index}.$extension"
  output_file_path="$data_folder$output_file_name"
  echo -e "$tag\n$trimmed_sequence" > "$output_file_path"
  if [ $IS_FASTQ ]; then
    echo -e "+\n$trimmed_quality\n" >> "$output_file_path"
  fi
  index=$((index+1))
done < "$1"

# Print 3-mer occurrences
#for kmer in "${!occ[@]}"; do
#  echo "$kmer: ${occ[$kmer]}"
#done

# Generate BWA index from reference genome
bwa index GCF_000146045.2_R64_genomic.fna

# Align files to reference genome and output SAM files
for file in "$data_folder"*.$extension; do
  echo "Writing to file: ${file%.*}.sam"
  bwa mem GCF_000146045.2_R64_genomic.fna "$file" > "${file%.*}.sam"
done

# Merge SAM files into a single file
merged_sam_file="$data_folder"merged.sam
for sam_file in "$data_folder"*.sam; do
  # Skip the merged SAM file
  if [[ "$sam_file" == "$merged_sam_file" ]]; then
    continue
  fi
  # Append the content of the SAM file to the merged SAM file
  cat "$sam_file" | grep -v '^@' >> "$merged_sam_file"
done

# Sort merged SAM file by chromosome and position
sort -k 3,3 -k 4,4n "$merged_sam_file" > "${merged_sam_file%.*}_sorted.sam"

# Count the number of aligned bases in the SAM file
aligned_count=$(grep -v '^@' "${merged_sam_file%.*}_sorted.sam" | cut -f6 | grep -oP "(\d+)M" | tr -d "M" | awk '{sum += $1} END {print sum}')

echo "Total aligned count: $aligned_count"