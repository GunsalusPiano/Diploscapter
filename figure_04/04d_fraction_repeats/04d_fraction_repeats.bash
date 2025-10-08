################################################################
# D. coronatus HiC PCA                                         #
################################################################
# run the Arima pipeline
# modified from the original at Arima’s GitHub
# https://github.com/ArimaGenomics/mapping_pipeline
# the output of this pipeline is a BAM file: Dco_hiC_rep1.bam
sbatch JOB_2024-11-18_003_Dco-HiC_Arima_pipe.sh

# The main result of the pipe is a bam file with mappings for both ends of the reads
# in ./deduplicated_files
# Use bam2pairs from pairix to convert into .pairs file
# Then use juicer_tools to convert .pairs into .hic
# nxDipCoro1_1.curated_primary.chrom.sizes contains 3 lines and 2 columns:
# First column being the chromosome name, second the size in bps
conda create -n pairix-env -c bioconda pairix
conda activate pairix-env
bam2pairs -l Dco_hiC_rep1.bam Dco_hiC_rep1
java -Xmx64g -jar juicer_tools.jar pre Dco_hiC_rep1.bsorted.pairs.gz Dco_hiC_rep1.bsorted.hic ./nxDipCoro1_1.curated_primary.chrom.sizes
conda deactivate

# Try principal component analysis with hicexplorer
#
# convert into .cool, resolution 100,000 bps
conda create -n hicexplorer-env -c bioconda hicexplorer
conda activate hicexplorer-env
hicConvertFormat --matrices Dco_hiC_rep1.bsorted.hic --inputFormat hic -r 100000 --outFileName Dco_hiC_rep1.bsorted.cool --outputFormat cool

# convert to h5
hicConvertFormat --matrices Dco_hiC_rep1.bsorted_100000.cool --inputFormat cool -r 100000 --outFileName Dco_hiC_rep1.bsorted_100000.h5 --outputFormat h5

# do a diagnostic before matrix correction
hicCorrectMatrix diagnostic_plot -m Dco_hiC_rep1.bsorted_100000.cool -o Dco_hiC_rep1.bsorted_100000.diag.png

# attempt to correct the matrix, using ICE
hicCorrectMatrix correct --matrix Dco_hiC_rep1.bsorted_100000.h5 --filterThreshold -3 3 --correctionMethod ICE --outFileName Dco_hiC_rep1.bsorted_100000.correctedICE.h5

# do the PCA
# 5 PCs
hicPCA --matrix Dco_hiC_rep1.bsorted_100000.correctedICE.h5 --outputFileName Dco_hiC_rep1.bsorted_100000.correctedICE.pca1.bedgraph Dco_hiC_rep1.bsorted_100000.correctedICE.pca2.bedgraph Dco_hiC_rep1.bsorted_100000.correctedICE.pca3.bedgraph Dco_hiC_rep1.bsorted_100000.correctedICE.pca4.bedgraph Dco_hiC_rep1.bsorted_100000.correctedICE.pca5.bedgraph --whichEigenvectors 1 2 3 4 5 --format bedgraph --ignoreMaskedBins

# The PC that corresponds to long-range inter-/intra-chromosome interactions is PC1
# Take the bedgraph file and convert it into the format similar to all the traits being graphed, ie:
# scaffold_position    PC1
# 50000                0
# 150000               0.01
# ...etc.-
# PCA_bedgraph_converter.py takes the average of all values within 100,000-bp windows and plots the average in the middle of the interval (ie. %100,000 = 50,000)
# nxDipCoro1_1.curated_primary.fa is the D. coronatus genome fasta
# 
# The command below should give the files "Dcor.SUPER_1.PC1.window100000.increment100000.txt" and "Dcor.SUPER_2.PC1.window100000.increment100000.txt"
python3 PCA_bedgraph_converter.py -f nxDipCoro1_1.curated_primary.fa -b Dco_hiC_rep1.bsorted_100000.correctedICE.pca1.bedgraph -x scaffold_position -y PC1 -p Dcor


################################################################
# D. pachys PoreC PCA                                          #
################################################################
# First, PoreC homozygous region matrix transformation
# This is necessary because a homozygous region exists between DpaA and DpaB. The homozygous region receives 2x the number of reads mapping.
# In order to do PCA, we randomly assigned reads mapping to this region to either DpaA or DpaB.
#
# Input file is "NlaIII_run01_Dpa-canu-het00c-YaHS-v202304-hethomosep_unphased.unsorted.pairs"
# which contains reads assigned to DpaA heterzygous region, DpaB heterozygous region, and the homozygous region
#
# This Python script randomly assign reads mapping to the homozygous region to either DpaA (scaffold_1) or DpaB (scaffold_2)
# The output of this script should be "Dpa-canu-het00c-YaHS-v202304.pairs"
python3 Dpa-canu-het00c-YaHS-v202304-hethomosep_PoreC_correction.py

# use pairtools (version 1.0.2) to sort the .pairs file
# pairtools: https://github.com/open2c/pairtools
conda create -n pairtools-env -c conda-forge -c bioconda pairtools
conda activate pairtools-env
pairtools sort Dpa-canu-het00c-YaHS-v202304.pairs

# use juicer_tools.jar pre to generate a hic file from the .sorted.pairs file
# juicer_tools.jar can be found at: https://github.com/aidenlab/JuicerTools
# chrom.sizes is just a tab-delimited chromosome_name, chromosome_size_in_bps
# (chromosomes separated by return character)
java -Xmx2g -jar juicer_tools.jar pre Dpa-canu-het00c-YaHS-v202304.sorted.pairs Dpa-canu-het00c-YaHS-v202304.hic ./chrom.sizes

# try principal component analysis with hicexplorer (version 3.7.2)
# convert into cool, resolution 100,000 bps
conda deactivate
conda activate hicexplorer-env
hicConvertFormat --matrices Dpa-canu-het00c-YaHS-v202304.hic --inputFormat hic -r 100000 --outFileName Dpa-canu-het00c-YaHS-v202304.cool --outputFormat cool

# convert to h5
hicConvertFormat --matrices Dpa-canu-het00c-YaHS-v202304_100000.cool --inputFormat cool -r 100000 --outFileName Dpa-canu-het00c-YaHS-v202304_100000.h5 --outputFormat h5

# do a diagnostic before matrix correction
hicCorrectMatrix diagnostic_plot -m Dpa-canu-het00c-YaHS-v202304_100000.cool -o Dpa-canu-het00c-YaHS-v202304_100000.diag.png

# attempt to correct the matrix, using the ICE algorithm
hicCorrectMatrix correct --matrix Dpa-canu-het00c-YaHS-v202304_100000.h5 --filterThreshold -3.5 3.5 --correctionMethod ICE --outFileName Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.h5

# attempt to do PCA, with 5 PCs
# output as .bedgraph file
hicPCA --matrix Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.h5 --outputFileName Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.pca1.begraph Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.pca2.bedgraph Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.pca3.bedgraph Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.pca4.bedgraph Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.pca5.bedgraph --numberOfEigenvectors 5 --format bedgraph –ignoreMaskedBins
conda deactivate

# Take the bedgraph file and convert it into the format similar to all the traits being graphed, ie:
# scaffold_position    PC3
# 50000                0
# 150000               0.01
# ...etc.-
# PCA_bedgraph_converter.py takes the average of all values within 100,000-bp windows and plots the average in the middle of the interval (ie. %100,000 = 50,000)
# Dpa-canu-het00c-YaHS-v202304.fasta is the genome file
# 
# The command below should give the files "Dpac.scaffold_1.PC3.window100000.increment100000.txt" and "Dpac.scaffold_2.PC3.window100000.increment100000.txt"
python3 PCA_bedgraph_converter.py -p Dpac -f Dpa-canu-het00c-YaHS-v202304.fasta -b Dpa-canu-het00c-YaHS-v202304_100000.correctedICE.pca3.bedgraph -w 100000 -i 100000 -x scaffold_position -y PC3