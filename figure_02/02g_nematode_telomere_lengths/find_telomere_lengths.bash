# these commands will use worm_telomere_lengths_extractor_v3.py to extract telomere lengths from long-read libraries

# This is how worm_telomere_lengths_extractor_v3.py extracts C. elegans VC2010/PD1074 telomere lengths from the Yoshimura PacBio dataset
# (SRA run accession SRR7594465, https://pubmed.ncbi.nlm.nih.gov/31123080/)
#
# (1) worm_telomere_lengths_extractor_v3.py first downloads SRR7594465 as a FASTA file (-a flag).
# (2) worm_telomere_lengths_extractor_v3.py runs TideHunter (https://github.com/Xinglab/TideHunter) to identify tandem repeats in the reads and their positions
#     - produces "SRR7594465.TideHunter_6-mers.txt" output to record the 6-mer repeats and their positions.
#     - the telomeric motif is set by the -t flag.
# (3) worm_telomere_lengths_extractor_v3.py assigns a 'teloscore' to each read, which is a dot product of occupancy and penalty_reward arrays, given "TTAGGC" telomeric motif. For example:
#     read sequence:      a  c  t  T  T  A  G  G  C  T  T  A  G  G  C  t  g  t  g  c  a  t  c  t  a  a  a  c  c  a  a  a  t  t  t  a  T  T  A  G  G  C  T  T  A  G  G  C
#     occupancy:        [ 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
#     penalty_reward:   [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
#
#     teloscore = occupancy ⋅ penalty_reward = 0
#
#     the extent of the -1 and the 1 in the penalty_reward array is indicated by the -n flag. This value is 15 for the example above (15 elements of -1 and +1)
#     The more positive this dot product is, the more likely it is a telomeric read. We call this dot product 'teloscore'.\
#     The teloscore cutoff to consider a read as telomeric can be set using the -c flag.
# (4) retain all reads with a teloscore greater than the -c flag.
# (5) use nucmer to find the extent of the subtelomere. The subtelomeric sequences are stored in the FASTA file indicated by the -s flag.
# (6) scan the reads using two adjacent windows (window width -w flag): where the occupancies of the telomeric motif (TTAGGC) are the most different is the extent of the telomere.
# (7) figure out if the extent of the telomere from step (5) matches step (6). Because of sequence errors etc they are often slightly different. This difference/leeway is specified
#     using the -d flag. Only if the difference is < this value is the telomere length considered.
# (8) the telomere lengths are reported for reads using values from step 6, in the *.telomereLengths.txt file.
# (9) extra info about the library can be provided in a text file using the -i (info) flag.

python3 worm_telomere_lengths_extractor_v3.py -a SRR7594465 -t TTAGGC -n 1000 -c 20 -s c_elegans.PRJNA13758.WS285.genomic.all_subtel2000.fasta -w 20 -d 50 -i library_info.SRR7594465.txt


# similarly for C. elegans CB4856 (SRA runs SRR8599835 to SRR8599843, Kim & al. https://pubmed.ncbi.nlm.nih.gov/31123081/). Concatenate these libraries first > SRR8599835-43.fasta and
# save as a local library. We can use this local library by invoking the -l flag.
# Analyze the C. briggsae QX1410 reads similarly (SRA run SRR17074503, Stevens & al. https://pubmed.ncbi.nlm.nih.gov/35348662/)
python3 worm_telomere_lengths_extractor_v3.py -l SRR8599835-43.fasta -t TTAGGC -n 1000 -c 20 -s Celegans_CB4856_GCA_004526295.1_ASM452629v1_genomic.identifiable_subtel2000.fasta -w 20 -d 50 -i library_info.SRR8599835-43.txt
python3 worm_telomere_lengths_extractor_v3.py -a SRR17074503 -t TTAGGC -n 1000 -c 20 -s Cbriggsae_QX1410_GCA_021491975.1_ASM2149197v1_genomic.identifiable_subtels.fasta -w 20 -d 50 -i library_info.SRR17074503.txt


# for Diploscapter, the worm_telomere_lengths_extractor_v3.py has an -r flag that enables the search for repetitive subtelomeres - as is the case in D. pachys and D. coronatus.
# for D. pachys, we used the Porechopped D. pachys PF1309 ONT R9.4 genomic reads.
# the -T flag specifies that the telomeric motifs are TAAGGG and TAAGGC, rather than the TAAGGGTAAGGC that TideHunter used to find tandem repeats.
python3 worm_telomere_lengths_extractor_v3.py -l Dpachys_porechopped.fasta -t TAAGGGTAAGGC -n 1000 -c 20 -s Dpa_subtelomeres.fasta -w 20 -d 50 -T TAAGGG,TAAGGC -i library_info.Dpachys.txt -r

# for D. coronatus, we used the D. coronatus PDL0010 PacBio HiFi reads.
python3 worm_telomere_lengths_extractor_v3.py -l nxDipCoro1_hifi.fasta -t TAAGGGTAAGGC -n 1000 -c 20 -s Dco_subtelomeres.fasta -w 20 -d 50 -T TAAGGG,TAAGGC -i library_info.Dcoronatus.txt -r

# concatenate all the telomere lengths results
cat *.telomereLengths.headerless.txt > all_nematode_telomeres.headerless.txt

# proceed to plot the distribution of telomere lengths in "all_nematode_telomeres.headerless.txt" using R
