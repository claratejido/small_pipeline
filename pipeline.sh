#!/bin/bash

DATA_DIR=./data

bwa index GCF_000146045.2_R64_genomic.fna

for file in [$DATA_DIR/*.fasta]
do
    bwa mem GCF_000146045.2_R64_genomic.fna $file > $file%.sam
done