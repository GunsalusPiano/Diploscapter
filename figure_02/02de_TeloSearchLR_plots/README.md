# Graphing the occupancy of telomeric repeats at the ends of long sequencing reads

Here, the *Diploscapter* telomeric repeat motif was found using TeloSearchLR ([Chung & al. 2025 *G3*](https://academic.oup.com/g3journal/article/15/6/jkaf062/8102967), [GitHub](https://github.com/gchchung/TeloSearchLR/)), which examines long sequencing reads for repeat motifs and their reverse complement motifs enriched at the ends. TeloSearchLR relies on the tandem repeat detection algorithm TideHunter ([Gao & al 2019 *Bioinformatics*](https://academic.oup.com/bioinformatics/article/35/14/i200/5529224), [GitHub](https://github.com/Xinglab/TideHunter)).

To graph the telomeric repeat occupancies at the ends of reads:
1. Run TeloSearchLR.py on both *D. pachys* and *D. coronatus* long reads. The commands are in ```TeloSearchLR_Dpa_commands.bash```.
2. With the occupancy data (tabular text files), run the R code in ```Fig2DE.md```. 
<img width="1004" height="906" alt="image" src="https://github.com/user-attachments/assets/7b729c03-686d-4bb2-b3c1-0e2bc61cca23" />
3. Move the figure to Adobe Illustrator for light editing for clarity.
