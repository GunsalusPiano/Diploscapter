# For this figure, the nuclear reference genome onto which PoreC reads were mapped was separated into three fragments:
# the homozygous region (1-17839851), DpaA heterozygous region, and DpaB heterozygous region.
# 
# The PoreC reads were mapped to this genome using the PoreC Snakemake pipeline: https://github.com/nanoporetech/Pore-C-Snakemake
# with the appropriate changes to the /config/references.tsv (gives the file path to the reference genome)
# and the /config/basecalls.tsv (gives file path to the PoreC reads)
#
# The reference genome, "Dpa-canu-het00c-YaHS-v202304-hethomosep.fasta", has the homozygous and heterozygous regions separated such that
# there are three FASTA entries in the file:
#     >SCAFFOLD_HOM
#     >SCAFFOLD_1_HET
#     >SCAFFOLD_2_HET
#
# Use the pipeline to also generate matrices for both HiCExplorer (.cool) and Juicer (.hic) suites of programs
snakemake --use-conda --cores 16 juicer
snakemake --use-conda --cores 16 cool

# Plot the HiC matrix
# using the hicexplorer suite of programs (https://github.com/deeptools/HiCExplorer)
conda activate hicexplorer-env

hicPlotMatrix --matrix NlaIII_run01_Dpa-canu-het00c-YaHS-v202304-hethomosep_unphased_250000.cool --outFileName NlaIII_run01_Dpa-canu-het00c-YaHS-v202304-hethomosep_unphased_250000.range_0-35.colorMap_binary.pdf --chromosomeOrder SCAFFOLD_HOM SCAFFOLD_1_HET SCAFFOLD_2_HET --vMin 0 --vMax 35 --colorMap binary
