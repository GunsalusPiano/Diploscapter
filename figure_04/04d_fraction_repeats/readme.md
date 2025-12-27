# Repeat annotation with RepeatScout, and quantifying the fraction of repetitive DNA in the genome

The .bash file contains the commands to run the custom Python script RepeatScout_repeat_occupancy_table.py to search for and annotate repetitive sequences in the *D. pachys* and the *D. coronatus* genomes. It broadly follows a repeat annotating strategy described by [Yoshida & al. (2023) *Nat Ecol Evol*](https://www.nature.com/articles/s41559-022-01980-z) for *Pristionchus* genomes by using [RepeatScout](https://github.com/Dfam-consortium/RepeatScout) and [RepeatMasker](https://www.repeatmasker.org/).
1. The script calls RepeatScout to identify repetitive DNA (requries RepeatScout to be in $PATH).
2. RepeatMasker is invoked to mark all genomic positions with repeats identified by RepeatScout (Requries RepeatMasker in &PATH).
3. Limit the analysis to repetitive DNA occuring 5 times or more in the genome, and the occupancy of these repeats are counted under two categories: "Simple_repeat" and "Unspecified" (used by RepeatScout).
4. The script outputs the sum of these two values along the length of the chromosome, one chromosome per output file.

Plotting of this data was done with the R ggplot2 code in ```genomic_features_visualizations.md``` in the main **Fig. 4** folder.
