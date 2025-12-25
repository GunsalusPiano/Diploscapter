# Read depths at the *D. pachys* and *D. coronatus* subtelomere

Requires the ```samtools depth``` output from **Fig. 1E** (*D. coronatus*) and **Fig. 1F** (*D. pachys*).


**For *D. pachys***:

The .bash file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 80,000 bp of DpaA and DpaB.

NB. "scaffold_1" and "scaffold_2" were the names of DpaA and DpaB respectively in an earlier draft assembly. This assembly had an error in DpaB, which meant the coordinates on the right end were 2364 bp off. This had no bearing on the analysis of subtelomeric repeats.

For both species:

The .Rmd file contains code in R to plot the read depth on top of a subtelomere repeat schematic, using standard ```ggplot2``` and the ```gggenes``` R package (https://github.com/wilkox/gggenes)
