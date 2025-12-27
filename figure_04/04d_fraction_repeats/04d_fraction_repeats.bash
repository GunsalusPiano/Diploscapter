##################################################################
# RepeatScout analysis of D. pachys nuclear genome               #
##################################################################
# 
# We used the script RepeatScout_repeat_occupancy_table.py to automate the search for repetitive sequences in the Dpa genome.
# Our use of RepeatScout broadly follows a repeat annotating method described by Kohta Yoshida & al. (2023) in Nat Ecol Evol
# for Pristionchus genomes. First, the script calls RepeatScout to identify repetitive DNA. RepeatMasker is invoked to mark
# all genomic positions with repeats identified by RepeatScout. The script then limits the analysis to repetitive DNA occuring
# 5 times or more in the genome, and the occupancy of these repeats are tabulated in the format below, where "Simple_repeat"
# and "Unspecified" are categories used by RepeatScout, and "sum" is just the sum of these two values.
# 
# coordinate    Simple_repeat    Unspecified      sum
# 50000	              0.00255        0.63467  0.63722
# 150000              0.00854        0.28424  0.29278
# 250000              0.00804        0.30466   0.3127
# ...etc.
#
# Occupancy (repetitive_DNA/all_DNA) is calculated for every 100,000-bp interval, and plotted at %100000 = 50,000 bp.
#
# arguments:
#     -f: scaffold_1.fasta = DpaA, scaffold_2.fasta = DpaB. Repeat libraries were generated separately per chromosome.
#     -p: output file prefix.
#     -l: l-mer table constructed by RepeatScout (=14).
#     -t: the threshold of repeat occurrence to consider (=5).
#     -w: window size (in bp) to consider for repeat occupancy (=100,000)
#     -i: increment (in bp) of the window (=100,000)
#
# requires:
#     RepeatScout in PATH
#     RepeatMasker in PATH
#
# Each D. pachys chromosome is analyzed separately, so the -f parameter has the individual D. pachys scaffolds listed separated by a comma
python3 RepeatScout_repeat_occupancy_table.py -f scaffold_1.fasta,scaffold_2.fasta -p Dpac -l 14 -t 5 -w 100000 -i 100000

##################################################################
# RepeatScout analysis of D. coronatus nuclear genome            #
##################################################################
#
# Same as above, but with D. coronatus scaffolds SUPER_1.fasta (DcoB) and SUPER_2.fasta (DcoA).

python3 RepeatScout_repeat_occupancy_table.py -f SUPER_1.fasta,SUPER_2.fasta -p Dcor -l 14 -t 5 -w 100000 -i 100000

