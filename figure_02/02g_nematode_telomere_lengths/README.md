# telomere lengths in *Diploscapter* compared to *Caenorhabditis*
1.  Use the instructions in ```find_telomere_lengths.bash``` to run ```worm_telomere_lengths_extractor_v3.py``` on long sequencing reads from *D. coronatus* PDL0010, *D. pachys* PF1309, *C. briggsae* QX1410 (Stevens et al. 2022 *PLOS Genet*), *C. elegans* CB4856 (Hawaiian, Kim et al. 2019 *Genom Res*), and *C. elegans* VC2010 (derived from the N2 reference, Yoshimura et al. 2019 *Genom Res*). The Python script
    - scans the sequencing reads for terminal tandem repeats of the telomeric motif.
    - attempts to match these reads to the known subtelomeric sequences
    - estimates the length of the telomeric array. If this coincides with where the subtelomere ends, then report the telomere length estimate.

2.  Use the R code in ```telomere_lengths_violin_plot.Rmd``` to generate a plot of the distribution of telomere lengths. The resulting graph:

<img width="1488" height="637" alt="image" src="https://github.com/user-attachments/assets/37606c5e-f998-4c7d-b91b-cc79e187ac43" />

