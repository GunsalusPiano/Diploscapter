# Visualizing the *D. pachys* PoreC data


```01d_Dpa_porec.bash``` has the commands and parameters to
  - run PoreC Snakemake - PoreC mapping pipeline adapted from [Oxford Nanopore Technologies](https://github.com/nanoporetech/Pore-C-Snakemake)). For this figure, the reference genome onto which PoreC reads were mapped contained three genome fragments: the homozygous region, the DpaA heterozygous region, and the DpaB heterozygous region. In an earlier draft assembly, scaffold_1 = DpaA, and scaffold_2 = DpaB.
  - ```hicPlotMatrix``` from hicexplorer v3.7.5 (https://github.com/deeptools/HiCExplorer) to generate publication-ready PoreC matrix plot
 
