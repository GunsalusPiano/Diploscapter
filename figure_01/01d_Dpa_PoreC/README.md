# Visualizing the *D. pachys* PoreC data


```01c_Dco_hic.bash``` has the commands and parameters to
  - run ```JOB_2024-11-18_003_Dco-HiC_Arima_pipe.sh``` - mapping pipeline adapted from [Arima](https://github.com/ArimaGenomics/mapping_pipeline). Downloads *D. coronatus* PDL0010 Hi-C reads from SRA (accession [ERR10228685](https://www.ncbi.nlm.nih.gov/sra/ERR10228685)), then maps the reads to the *D. coronatus* genome assembly "nxDipCoro1_1.curated_primary.fa"
  - use ```bam2pairs``` from pairix v0.3.8 (https://github.com/4dn-dcic/pairix) to turn .bam to .pairs
  - use ```juicer_tools pre``` from juicer_tools v1.22.01 (https://github.com/aidenlab/JuicerTools) to turn .pairs into .hic for initial interactive visualization in Juicer.
  - use ```hicConvertFormat``` from hicexplorer v3.7.5 (https://github.com/deeptools/HiCExplorer) to generate publication-ready Hi-C matrix plot
  - use ```hicPCA``` from hicexplorer to do principal component analysis on the Hi-C matrix (for **Fig. 4B**)
Code for producing the PoreC plot in HiCExplorer (https://hicexplorer.readthedocs.io/en/latest/).
Also requires Oxford Nanopore's Pore-C Snakemake (https://github.com/nanoporetech/Pore-C-Snakemake).
