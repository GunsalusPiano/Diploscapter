**Read depths at the *D. pachys* subtelomere**

Requires the ```samtools depth``` output from **Fig. 1F**.

1. The .bash file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 80,000 bp of DpaA and DpaB.
2. The .Rmd file contains code in R to plot the read depth on top of a subtelomere repeat schematic, using the ```gggenes``` R package (https://github.com/wilkox/gggenes)

NB. "scaffold_1" and "scaffold_2" refer to DpaA and DpaB respectively. The analysis was done on an earlier draft of the genome assembly where an assembly error in DpaB meant the coordinates on the right end were 2364 bp off. This has no bearing on the analysis of subtelomeric repeats. 
