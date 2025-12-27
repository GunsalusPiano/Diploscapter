# *Diploscapter* chromosome-to-chromosome comparisons and RAxML analysis

We turn to nucmer ([mummer suite](https://github.com/mummer4/mummer)) again for chromosome-scale alignments. But first, prep work.
- For this analysis only, DpaA's inversion at 41921280-42998403 with respect to every other *Diploscapter* chromosome is reverted.
- Separate the chromosomes into individual FASTA files named DpaA.fasta, DpaB.fasta, DcoA.fasta, and DcoB.fasta.

```nucmer_percentage_id_along_length_v5_for_Dpa_vs_Dco.py``` helps run nucmer in a pairwise fashion. To focus only on alignments that are strictly between syntenic regions of the chromosomes, the script applies ```delta-filter -g``` ('global' alignment). It also calculates the percentage identity between two stretches of DNA for every 100,000 bp window of the reference DNA, then writes that percentage to a file. Thus the percentage identity is (number of query bp identical to reference)/(number of reference bp alignable), every 100,000 bp. A problem could arise if the denominator here is very low, perhaps due to the number of Ns (scaffolding gaps), or simply that the query and the reference cannot be confidently aligned. In this case, the script limits the reporting of percent identity if the number of alignable bp in the reference, the denominator, is at least 70,000, or 70% of the window.     
