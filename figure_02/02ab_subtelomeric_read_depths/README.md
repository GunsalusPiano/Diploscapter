# Read depths at the *D. pachys* and *D. coronatus* subtelomeres

Requires the ```samtools depth``` output from **Fig. 1E** (*D. coronatus*) and **Fig. 1F** (*D. pachys*).

## For *D. coronatus*:
The ```Dco_get_subtelomeric_read_depths.bash``` file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 40,000 bp of DcoA and DcoB.

NB. At the time of assembly, DcoA was given the name SUPER_2, and DcoB = SUPER_1.

## For *D. pachys*:

The .bash file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 80,000 bp of DpaA and DpaB.

NB. "scaffold_1" and "scaffold_2" were the names of DpaA and DpaB respectively in an earlier draft assembly. This assembly had an error in DpaB, which meant the coordinates on the right end were 2364 bp off. Other than a coordinate change, this had no bearing on the analysis of subtelomeric repeats.

## For both species:

The R code in ```fig_02_subtel_read_depth_and_repeat_annotations(2024-12).Rmd``` takes the read depths from above, as well as subtelomere unit boundaries annotated by hand, and produces the read depth plots on top of a subtelomere repeat schematics. The R code uses ```ggplot2``` and the ```gggenes``` R packages (https://github.com/wilkox/gggenes).
