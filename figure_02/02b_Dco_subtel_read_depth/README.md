**Read depths at the *D. coronatus* subtelomere**

Requires the ```samtools depth``` output from **Fig. 1E**.

1. The ```Dco_get_subtelomeric_read_depths.bash``` file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 40,000 bp of DcoA and DcoB.
2. The .Rmd file contains code in R to plot the read depth on top of a subtelomere repeat schematic, using the ```gggenes``` R package (https://github.com/wilkox/gggenes)

NB. At the time of assembly, DcoA was given the name SUPER_2, and DcoB = SUPER_1.
