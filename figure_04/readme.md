# Fig. 4: Signatures of ancestral chromosome arrangements.

The ancestral linkage groups mapped to the the *Diploscapter* chromosomes, and several other molecular traits of the chromosomes.   
- Panel A: 'Nigon painting' done in R based on code by [Gonzalez & al. (2021) *G3*](https://doi.org/10.1093/g3journal/jkaa020)([GitHub link](https://github.com/pgonzale60/vis_ALG))
- Panel B: PCA of Hi-C/PoreC matrices
  - Analysis done with hicPCA from [HiCExplorer](https://github.com/deeptools/HiCExplorer). For _D. pachys_, an additional step was required to reassign PoreC reads mapping to the homozygous region to either DpaA or DpaB.
  - Convert .bed format into x, y coordinates
  - Plot with R ggplot2
- Panel D: custom Python script uses RepeatScout + RepeatMasker for repeat annotations. Plotting with R ggplot2.
- Panels C, E, F: bioinformatic analyses done with Python scripts. Plotting with R ggplot2. 
