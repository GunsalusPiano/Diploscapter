# Principal components in the *Diploscapter* contact data

For *D. coronatus*:
- The ```04b_HiC_PoreC_PCA.bash``` file contains the run commands to get the PCs for contact data using hicexplorer (https://hicexplorer.readthedocs.io/en/latest/).
- PC1 was used for both DcoA and DcoB

For *D. pachys*:
- First, the .pairs file from the PoreC Snakemake pipeline was modified such that reads mapping to the homozygous region (identical between DpaA and DpaB) was reassigned to *either* DpaA or DpaB. This was done with ```Dpa-canu-het00c-YaHS-v202304-hethomosep_PoreC_correction.py```. The process is described in more detail in ```2023-05-04_D_pachys_Pore_C_transformation.pdf```.
- Then, the re-assigned .pairs file was turned into .hic and .cool for visualization and PCA. Commands in ```04b_HiC_PoreC_PCA.bash```

Plotting of this data was done with the R ggplot2 code in ```genomic_features_visualizations.md``` in the main **Fig. 4** folder.
