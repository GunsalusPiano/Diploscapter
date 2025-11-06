#!/bin/bash

# first map reads to the genome
# READS = D. coronatus HiFi genomic reads
READS="nxDipCoro1_hifi.fastq"

# GENOME = genome assembly of D. coronatus
GENOME="nxDipCoro1_1.curated_primary.fa"

# map reads to the genome, using minimap2 v2.28-r1209
# Sort SAM then turn into BAM BAM
# Sort and index alignments
ALIGNMENTS="Dco.aln"

minimap2 -ax map-hifi ${GENOME} ${READS} > ${ALIGNMENTS}.sam
samtools sort -@ 16 -o ${ALIGNMENTS}.sorted.bam ${ALIGNMENTS}.sam
samtools index -@ 16 ${ALIGNMENT}.sorted.bam ${ALIGNMENT}.sorted.bam.bai

# Use samtools depth to get the depth
samtools depth -aa -o ${ALIGNMENT}.sorted.samtools_depth.txt -@ 14 ${ALIGNMENT}.sorted.bam

# This depth file is very big as it contains a value for every genomic position
# Take the arithmetic mean of depths every 100,000 bps using Python script samtools_depth_output_simplifier_v4.py
python3 samtools_depth_output_simplifier_v4.py -d ${ALIGNMENT}.sorted.samtools_depth.txt -w 100000

# The output should be ${ALIGNMENT}.sorted.samtools_depth.window100000.txt
# Use R to graph the read coverage