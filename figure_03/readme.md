# Plotting the chromosome-to-chromosome alignment in *Diploscapter*
The chromosome-to-chromosome alignment was performed using nucmer (from [mummer](https://github.com/mummer4/mummer) suite of programs), then plotted.
- First, separate the four *Diploscapter* chromosomes into their own files: DpaA.fasta, DpaB.fasta, DcoA.fasta, DcoB.fasta.
- List them in the file ```list_of_FASTAs.txt```, one file name per line.
- ```repeat_nucmer_and_delta_simplified_dotplot_prep.py``` was used to automate the nucmer and show-coords commands.
```shell
python3 repeat_nucmer_and_delta_simplified_dotplot_prep.py -f list_of_FASTAs.txt # run command
```
- The Python script
  - runs chromosome-to-chromosome ```nucmer``` comparisons, ignoring self-to-self alignment and using the fasta files listed in ```list_of_FASTAs.txt```
  - uses the ```delta-filter``` program to eliminate spurious alignments, with the -l flag set at 10000.
  - interprets the resulting .delta file, then outputs a more human-readable tabular file (.coords) that lists the start and end coordinates of the alignment query and reference.
- the R code in ```Diploscapter_collinearity.md``` was used to produce the alignment dotplots.


The resulting matrix of dot plots below.
<img width="915" height="913" alt="image" src="https://github.com/user-attachments/assets/27efbc14-a219-49ea-916c-bdf3e9167ca9" />


