# *Diploscapter* chromosome-to-chromosome comparisons and RAxML analysis

## Chromosome-to-chromosome alignment and percentage identity
We turn to nucmer ([mummer suite](https://github.com/mummer4/mummer)) again for chromosome-scale alignments. But first, prep work.
- For this analysis only, DpaA's inversion at 41921280-42998403 with respect to every other *Diploscapter* chromosome is reverted.
- Separate the chromosomes into individual FASTA files named DpaA.fasta, DpaB.fasta, DcoA.fasta, and DcoB.fasta.

```nucmer_percentage_id_along_length_v5_for_Dpa_vs_Dco.py``` helps run nucmer in a pairwise fashion. To focus only on alignments that are strictly between syntenic regions of the chromosomes, the script applies ```delta-filter -g``` ('global' alignment). It also calculates the percentage identity between two stretches of DNA for every 100,000 bp window of the reference DNA, then writes that percentage to a file. Thus the percentage identity is (number of query bp identical to reference)/(number of reference bp alignable), every 100,000 bp. A problem could arise if the denominator here were very low, perhaps due to the number of Ns (scaffolding gaps), or simply that the query and the reference could not be confidently aligned. In this case, the script limits the reporting of percent identity if the number of alignable bp in the reference, the denominator, is at least 70,000, or 70% of the window. The percentage identities are concatenated in the file ```all_vs_all.comparison.percentage_id_along_length.window_100000.sample_freq_50000.reject_threshold_0.7.txt```, and the R code in ```nucmer_percentage_id_along_length_visualization.md``` was used to plot the data.

## RAxML analysis
Prep work:
- Transcript sequences, by chromosome (DpaA, DpaB, DcoA, DcoB).
- Protein sequences, by chromosome.

For both, also include data from _C. elegans_ (from WormBase WS285) as an outgroup for these analyses.

The Python script ```chromosome_wide_MSA_of_BUSCOs_2024-01-02_update.py``` takes care of the following steps
- Fetch "Complete" orthologs from the BUSCO table (data from **Fig. 4A**) - that is, the ortholog of the BUSCO exists on the chromosome/genome.
- Align the protein sequences of the BUSCO orthologs using [MUSCLE](https://github.com/rcedgar/muscle).
- Turn the amino acid alignments into CDS alignment
- Concatenates all the CDS alignments
- Run [modeltest-ng](https://github.com/ddarriba/modeltest) to test the best fit model for nucleotide evolution. First, second and third nucleotides of codons are considered separately.
- Run [RAxML-ng](https://github.com/amkozlov/raxml-ng) using the models determined above.

Furthermore, the script considers the Dpa homozygous region separately from the heterozygous region. At the end of the runs there should be *.raxml.support files (RAxML trees with bootstrap support, 300 replicates for the homozygous region, 1000 replicates for the heterozygous region). [FigTree](https://tree.bio.ed.ac.uk/software/figtree/) was used to plot the phylogenetic relationships between DpaA, DpaB, DcoA, DcoB.
