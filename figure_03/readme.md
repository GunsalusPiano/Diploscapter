# Plotting the chromosome-to-chromosome alignment in *Diploscapter*
The chromosome-to-chromosome alignment was performed using nucmer (from [mummer](https://github.com/mummer4/mummer) suite of programs), and the resulting alignments were plotted using the ```show-coords``` output.
- First, separate the four *Diploscapter* chromosomes into their own files: DpaA.fasta, DpaB.fasta, DcoA.fasta, DcoB.fasta.
- List them in the file ```list_of_FASTAs.txt```, one file name per line.
- ```repeat_nucmer_and_delta_simplified_dotplot_prep.py``` was used to automate the nucmer and show-coords commands.
  - ```shell
    python3 repeat_nucmer_and_delta_simplified_dotplot_prep.py -f list_of_FASTAs.txt # run command
    ```
  - The script runs chromosome-to-chromosome nucmer comparisons, ignoring self-to-self alignment.
  - The script uses the ```delta-filter``` program to eliminate spurious alignments, with the -l flag set at 10000.
  - T
- 

