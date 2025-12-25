# samtools_depth_output_simplifier_v4.py takes a samtools depth output (specified with -d) and simplifies it
# by taking the average every 'w' nucleotides (window, -w flag). So a samtools depth file with 1,000,000 entries
# corresponding to 1,000,000 bp) can be simplied to 1,000,000/w entries.
# Furthermore, a file specifying the ranges of coordinates is supplied (-c flag) with the following tab-delim format
# 
# scaffold/contig	start_position	end_position
#
# the position values are 1-based, inclusive on both ends.

# from aligning Dpa ONT R9.4 genomic reads to the assembly (data from Fig. 1F, = Dpa.aln.sorted.samtools_depth.txt))
# use samtools_depth_output_simplifier_v4.py to simplify read depth values, and only focus on the terminal 80 kb
# specified in the Dco_chromosome_ranges.txt
python3 samtools_depth_output_simplifier_v4.py -d Dpa.aln.sorted.samtools_depth.txt -w 10 -c Dpa_chromosome_ranges.txt

# proceed to plot the read depths using R ggplot