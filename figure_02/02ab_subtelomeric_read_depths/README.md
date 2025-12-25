# Read depths at the *D. pachys* and *D. coronatus* subtelomeres

Requires the ```samtools depth``` output from **Fig. 1E** (*D. coronatus*) and **Fig. 1F** (*D. pachys*).

## For *D. coronatus*:
The ```Dco_get_subtelomeric_read_depths.bash``` file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 40,000 bp of DcoA and DcoB, the ranges specified in ```Dco_chromosome_ranges.txt```. Read depth values are reported every 10 bp (-w 10 in ```samtools_depth_output_simplifier_v4.py```).

NB. At the time of assembly, DcoA was given the name "SUPER_2", and DcoB = "SUPER_1".

## For *D. pachys*:

The ```Dpa_get_subtelomeric_read_depths.bash``` file describes using ```samtools_depth_output_simplifier_v4.py``` to extract read depth data for the terminal 80,000 bp of DpaA and DpaB, the ranges specified in ```Dpa_chromosome_ranges.txt```. Read depth values are reported every 10 bp (-w 10 in ```samtools_depth_output_simplifier_v4.py```).

NB. "scaffold_1" and "scaffold_2" were the names of DpaA and DpaB respectively in an earlier draft assembly. This assembly had an error in DpaB, which meant the coordinates on the right end were 2364 bp off. Other than a coordinate change, this had no bearing on the analysis of subtelomeric repeats.

## For both species:

The R code in ```fig_02_subtel_read_depth_and_repeat_annotations-2024-12-.md``` takes the read depths from above, as well as subtelomere unit boundaries annotated by hand (```Dpa_subtelomere_features.txt``` and ```Dco_subtelomere_features.txt```), and produces the read depth plots on top of subtelomere repeat schematics. The R code uses ```ggplot2``` and the ```gggenes``` R packages (https://github.com/wilkox/gggenes).

The plots produced are in ```Dpa_chrom_end_depth_and_subtel_repeats_wide_short.pdf```
<img width="1481" height="516" alt="image" src="https://github.com/user-attachments/assets/3258dd36-66ec-45b0-9d1f-d94d1bc51945" />

and in ```Dco_chrom_end_depth_and_subtel_repeats_skinny_short.pdf```
<img width="1451" height="557" alt="image" src="https://github.com/user-attachments/assets/0d4a0f5e-2946-4d9c-b4fc-8937a3c10d64" />

These were lightly edited for clarity in Adobe Illustrator before publication.
