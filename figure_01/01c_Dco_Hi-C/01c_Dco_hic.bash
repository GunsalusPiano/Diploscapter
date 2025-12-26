# run the Arima pipeline
# modified from the original at Arima’s GitHub
# https://github.com/ArimaGenomics/mapping_pipeline
#
# input files:
#   -the hic library reads, SRA run ERR10228685
#   -reference genome, nxDipCoro1_1.curated_primary.fa
#   -
# requires:
#   -bwa v0.7.17
#   -picard v2.27.5
#   -samtools v1.20
sbatch JOB_2024-11-18_003_Dco-HiC_Arima_pipe.sh

# The main result of the pipe is a bam file with mappings for both ends of the reads
# in ./deduplicated_files
# Use bam2pairs from pairix to convert into .pairs file
# Then use juicer_tools to convert .pairs into .hic
# nxDipCoro1_1.curated_primary.chrom.sizes contains 3 lines and 2 columns:
# First column being the chromosome name, second the size in bps
# requires:
#   -pairix v0.3.8 (https://github.com/4dn-dcic/pairix)
#   -juicer_tools v1.22.01 (https://github.com/aidenlab/JuicerTools)

# bam2pairs
cp ./deduplicated_files/Dco_hiC_rep1.bam .
bam2pairs -l Dco_hiC_rep1.bam Dco_hiC_rep1

# juicer_tools
java -Xmx64g -jar juicer_tools.jar pre Dco_hiC_rep1.bsorted.pairs.gz Dco_hiC_rep1.bsorted.hic ./nxDipCoro1_1.curated_primary.chrom.sizes

# Try principal component analysis and plotting the matrix with hicexplorer
# convert into .cool, resolution 100,000 bps
# requires:
#   -hicexplorer v3.7.5 (https://github.com/deeptools/HiCExplorer)
conda create -n hicexplorer-env -c bioconda hicexplorer
conda activate hicexplorer-env
hicConvertFormat --matrices Dco_hiC_rep1.bsorted.hic --inputFormat hic -r 100000 --outFileName Dco_hiC_rep1.bsorted.cool --outputFormat cool
hicConvertFormat --matrices Dco_hiC_rep1.bsorted.hic --inputFormat hic -r 250000 --outFileName Dco_hiC_rep1.bsorted.cool --outputFormat cool

# plot the hiC matrix at the 250,000 bp resolution
hicPlotMatrix -m Dco_hiC_rep1.bsorted_250000.cool --outFileName Dco_hiC_rep1.bsorted_250000.cool.range_0-2000.colorMap_binary.pdf --vMin 0 --vMax 2000 --colorMap binary

# Do the PCA using the 100,000-bp resolution matrix
# convert to h5
hicConvertFormat --matrices Dco_hiC_rep1.bsorted_100000.cool --inputFormat cool -r 100000 --outFileName Dco_hiC_rep1.bsorted_100000.h5 --outputFormat h5

# do a diagnostic before matrix correction
hicCorrectMatrix diagnostic_plot -m Dco_hiC_rep1.bsorted_100000.cool -o Dco_hiC_rep1.bsorted_100000.diag.png

# attempt to correct the matrix, using ICE
hicCorrectMatrix correct --matrix Dco_hiC_rep1.bsorted_100000.h5 --filterThreshold -3 3 --correctionMethod ICE --outFileName Dco_hiC_rep1.bsorted_100000.correctedICE.h5

# Do a PCA on the hiC data
# grab the first five PCs

hicPCA --matrix Dco_hiC_rep1.bsorted_100000.correctedICE.h5 --outputFileName Dco_hiC_rep1.bsorted_100000.correctedICE.pca1.bedgraph Dco_hiC_rep1.bsorted_100000.correctedICE.pca2.bedgraph Dco
_hiC_rep1.bsorted_100000.correctedICE.pca3.bedgraph Dco_hiC_rep1.bsorted_100000.correctedICE.pca4.bedgraph Dco_hiC_re
p1.bsorted_100000.correctedICE.pca5.bedgraph --whichEigenvectors 1 2 3 4 5 --format bedgraph --ignoreMaskedBins

# the 1st PC is plotted in Figure 4 b, iii and iv