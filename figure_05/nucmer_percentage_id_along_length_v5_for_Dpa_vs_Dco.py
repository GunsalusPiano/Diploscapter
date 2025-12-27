#############################################################################################################
# nucmer_percentage_id_along_length.py
# determines the percentage identity of two sequences aligned by nucmer
# 
# inputs:
#    - the individual FASTA files of the alignment
#    - the alignment delta (likely filtered by -q and -r ?)
#    - the alignment show-snps output (which identifies all the mismatches and indels in the alignment)
#    - the alignment show-coords output
#
# params:
#    - window size: the range from which we can calculate % identity
#
# outputs:
#    - a text file with % identity along the length of the 'reference' 
#
#
# Identity calculations:
#    reference: ATACGG-ACC
#               || ||  |||
#    query:     AT-CGATACC
#    This would be considered identity = 0.7 with a window_size = 9. Define the denominator as the total number (10) of nucleotide spaces to match
#    even though both reference and query have 9 nucleotides. Thus 7 identical out of 10 spaces, or 7/10 = 0.7
#
#
# Search strategy:
#    With repeated regions, nucmer has a few options as to which alignment to use. Some observations:
#       - With delta-filter options -q or -r, the filter picks the best length*identity. This may seem like a good option, but sometimes the best
#         length*identity match is elsewhere on the genome, rather than the collinear/homologous locus on the genome.
#       - We can avoid repeats altogether by using the show-snps with option -C. Problem with this is that it's hard to tell what range of nucleotides
#         were aligned. The show-snps program does have a -S function to allow the show-coords output to be piped into show-snps, but it is POORLY DOCUMENTED
#         and the program says "ERROR: Could not parse input from stdin" when I attempted to pipe in a show-coords output.
#       - An alternative: use delta-filter with option -g (genome mode). Then generate the snps from this filtered delta.
#            - Scan through show-snps output
#            - Save the alignments as a nested list:
#                 - [ [# of positions, # of matches]_1, [# of positions, # of matches]_2, [# of positions, # of matches]_3 ... etc ]   
#            - Calculate the % identity etc.
##################################################################################################################

from Bio import SeqIO
import math
import os

# generate_identity_array(args)
# 
def generate_identity_array(show_coords_output, show_snps_output, nucmer_reference_fasta):
	#len(SeqIO.read(nucmer_reference_fasta, "fasta")) + 1
	
	# identity array code:
	#   "N" = Not aligned
	#   "M" = aligned and Matched
	#   "S" = aligned but miSmatched (SNP)
	#   "I" = aligned but Insertion on ref, deletion on query
	#   "D" = aligned but Deletion on ref, insertion on query
	
	# build identity array: first fill with "N". WE will build an array 1 unit larger than the sequence length to take advantage of 1-based index
	
	identity_array = ["N" for item in range(len(SeqIO.read(nucmer_reference_fasta, "fasta")) + 1)]
	
	# first modify the 0th element to become a dictionary - this will store all the positions with deletions and the number of nucleotide positions for alignment
	
	identity_array[0] = {}
	
	
	
	# then, modify the array to change "N" at all positions that are aligned to "M".
	
	with open(show_coords_output, "r") as show_coords_fh:
		show_coords_fh.readline()
		show_coords_fh.readline()
		show_coords_fh.readline()
		show_coords_fh.readline()
		
		show_coords_line = show_coords_fh.readline()
		while show_coords_line:
			alignment_start = int(show_coords_line.split()[0])
			alignment_end = int(show_coords_line.split()[1])
			
			for position in range(alignment_start, alignment_end + 1):
				identity_array[position] = "M"
				
			show_coords_line = show_coords_fh.readline()
	
	# now, modify the array to change "M" at all mismatched positions to "S", "I", or "D"
	
	with open(show_snps_output, "r") as show_snps_fh:
		show_snps_fh.readline()
		show_snps_fh.readline()
		show_snps_fh.readline()
		show_snps_fh.readline()
		
		show_snps_line = show_snps_fh.readline()
		
		while show_snps_line:
			snp_position = int(show_snps_line.split()[0])
			
			# mismatch (letter to letter, not indel)
			if show_snps_line.split()[1] != "." and show_snps_line.split()[2] != ".":
				identity_array[snp_position] = "S"
				show_snps_line = show_snps_fh.readline()
				
			# if insertion on the reference, deletion on the query. show_snps_line.split()[1] should be a letter but show_snps_line.split()[2] should be a "."
			elif show_snps_line.split()[1] != "." and show_snps_line.split()[2] == ".":
				identity_array[snp_position] = "I"
				show_snps_line = show_snps_fh.readline()
				
			# if deletion on the reference, insertion on the query. show_snps_line.split()[1] should be "."
			elif show_snps_line.split()[1] == ".":
				identity_array[snp_position] = "D"
				
				number_of_positions = 2
				
				next_snp_line = show_snps_fh.readline()
				next_snp_position = int(next_snp_line.split()[0])
				
				while next_snp_position == snp_position:
					number_of_positions = number_of_positions + 1
					next_snp_line = show_snps_fh.readline()
					next_snp_position = int(next_snp_line.split()[0])
				
				identity_array[0][snp_position] = number_of_positions
				show_snps_line = next_snp_line
			#print(snp_position)
			
		#print("done modifying id array")
	return identity_array

# identity_array_to_percent_id
# takes the identity array generated by the function above
# reports the % identity at every sample_frequency using window_size
# reject threshold: the fraction of the window that is aligned. If < reject threshold, then do not calculate the %id
def identity_array_to_percent_id(identity_array, window_size, sample_frequency, reject_threshold, output_file):
	# step from 1x sample_frequency, 2 x sample_frequency, 3 x sample_frequency etc. until end of identity array
	# take value from sample_frequency +/- window_size/2
	number_of_positions = 0
	number_of_identities = 0
	
	length_of_sequence = len(identity_array)
	
	with open(output_file, "w") as output_fh:
		sample_position = sample_frequency
		
		sample_lower_bound = 0
		sample_upper_bound = 0
		
		if window_size % 2 == 0:
			sample_lower_bound = int(window_size/2)
			sample_upper_bound = int(window_size/2)
		else:
			sample_lower_bound = math.floor(window_size/2)
			sample_upper_bound = math.floor(window_size/2) + 1
			
		position_dictionary = {"N": 0, "M": 1, "S": 1, "I": 1, "D": 0}
		identity_dictionary = {"N": 0, "M": 1, "S": 0, "I": 0, "D": 1}
				
		while sample_position < length_of_sequence:
			for position in range(max(1, sample_position - sample_lower_bound), min(sample_position + sample_upper_bound, length_of_sequence)):
				if identity_array[position] == "D":
					position_dictionary["D"] = identity_array[0][position]
				else:
					pass
				number_of_positions = number_of_positions + position_dictionary[identity_array[position]]
				number_of_identities = number_of_identities + identity_dictionary[identity_array[position]]
			
			line_to_print = ""
			
			if number_of_positions < reject_threshold * window_size:
				line_to_print = str(sample_position) + "\t" + str(number_of_positions) + "\t" + str(number_of_identities) + "\n"
			else:
				line_to_print = str(sample_position) + "\t" + str(number_of_positions) + "\t" + str(number_of_identities) + "\t" + str(100*number_of_identities/number_of_positions) + "\n"
			
			output_fh.write(line_to_print)
			
			sample_position = sample_position + sample_frequency
			number_of_positions = 0
			number_of_identities = 0
	
def write_identity_array_to_files(identity_array, deletion_positions_filename, identity_array_filename):
	with open(identity_array_filename, "w") as identity_array_fh, open(deletion_positions_filename, "w") as deletion_positions_fh:
		for deletion_position in identity_array[0]:
			line_to_print = str(deletion_position) + "\t" + str(identity_array[0][deletion_position]) + "\n"
			deletion_positions_fh.write(line_to_print)
		
		for index, item in enumerate(identity_array):
			if index == 0:
				pass
			else:
				line_to_print = item + "\n"
				identity_array_fh.write(line_to_print)

def write_files_to_identity_array(deletion_positions_filename, identity_array_filename):
	identity_array = []
	identity_array.append({})
	with open(identity_array_filename, "r") as identity_array_fh, open(deletion_positions_filename, "r") as deletion_positions_fh:
		line = deletion_positions_fh.readline().rstrip()
		while line:
			key = int(line.split()[0])
			value = int(line.split()[1])
			identity_array[0][key] = value
			line = deletion_positions_fh.readline().rstrip()
		
		line = identity_array_fh.readline().rstrip()
		
		while line:
			identity_array.append(line)
			line = identity_array_fh.readline().rstrip()
		
		return identity_array
		
# nucmer_percentage_id_along_length_v3_output_prepper_for_R
# the main prepper of text files for R
#   reference_chromosome is the reference chromosome name
#   list_of_chromosomes_compared is the list of chromosomes that have been nucmer-ed against reference chromosome
#   comparison_info is all the text in the .txt output filename
# THIS VERSION IS EXPLICITLY FOR Dpa_scaffold_2 IN Dpa-canu-het00c-YaHS where Dpa_scaffold_2 is actually missing the left portion
# DO NOT USE OTHERWISE
#def nucmer_percentage_id_along_length_v3_output_prepper_for_R(reference_chromosome, list_of_chromosomes_compared, comparison_info):
#	
#	output_filename = reference_chromosome + "_vs_all." + comparison_info
#	
#	# this was necessary because Dpa_scaffold_2 was missing the left portion (homozygous)
#	if reference_chromosome == "Dpa_scaffold_2":
#		output_filename = reference_chromosome + "_adjusted_vs_all." + comparison_info
#	else:
#		pass
#	
#	with open(output_filename, 'w') as output_fh:
#		line_to_print = "reference\tquery\tcoordinate\tevaluated\tmatched\tpercentage_id\n"
#		output_fh.write(line_to_print)
#		
#		for query_chromosome in list_of_chromosomes_compared:
#			input_filename = reference_chromosome + "_vs_" + query_chromosome + "." + comparison_info
#			with open(input_filename, 'r') as input_fh:
#				line = input_fh.readline()
#				
#				while line:
#					if reference_chromosome != "Dpa_scaffold_2":
#						line_to_print = reference_chromosome + "\t" + query_chromosome + "\t" + line
#						output_fh.write(line_to_print)
#					else:
#						line_items = line.split()
#						
#						coordinate = int(line_items[0])
#						coordinate = coordinate + 17800000
#						
#						line_to_print = ""
#						
#						if len(line_items) == 3:
#							line_to_print = reference_chromosome + "\t" + query_chromosome + "\t" + str(coordinate) + "\t" + line_items[1] + "\t" + line_items[2] + "\n"
#						elif len(line_items) == 4:
#							line_to_print = reference_chromosome + "\t" + query_chromosome + "\t" + str(coordinate) + "\t" + line_items[1] + "\t" + line_items[2] + "\t" + line_items[3] + "\n"
#						
#						output_fh.write(line_to_print)
#					
#					line = input_fh.readline()


# nucmer_percentage_id_along_length_v3_output_prepper_for_R
# the main prepper of text files for R
#	Takes the nucmer_percentage_id_along_length output file and addes a header:
#	reference    query    coordinate    evaluated    matched    percentage_id_along_length
#
#	and two columns:
#   	reference_chromosome: is the reference chromosome name
#		query_chromosome: query chromosome name
#
#  	list_of_chromosomes_compared is the list of query chromosomes that have been nucmer-ed against reference chromosome
#   comparison_info is all the text in the .txt output filename
		
def nucmer_percentage_id_along_length_v5_output_prepper_for_R(ref_scaffolds_list, query_scaffolds_list, comparison_info):
	
	# define the output file name
    output_filename = "all_vs_all." + comparison_info
        
    # write to the output file
    with open(output_filename, 'w') as output_fh:
		
		# appropriate header for R to recognize
        line_to_print = "reference\tquery\tcoordinate\tevaluated\tmatched\tpercentage_id\n"
        output_fh.write(line_to_print)
		
		# go through all the original nucmer_percentage_id_along_length output files
		# and write those additional 2 columns (reference, query) for R
        
        for ref_scaffold in ref_scaffolds_list:
            for query_scaffold in query_scaffolds_list:
                if ref_scaffold == query_scaffold:
                    pass
                else:
                    input_filename = ref_scaffold + "_vs_" + query_scaffold + "." + comparison_info
                    
                    with open(input_filename, "r") as input_fh:
                        line = input_fh.readline()
                        
                        while line:
                            line_to_print = ref_scaffold + "\t" + query_scaffold + "\t" + line
                            output_fh.write(line_to_print)
                            
                            line = input_fh.readline()


def main():
	list_of_scaffolds = ["DpaA", "DpaB", "DcoA", "DcoB"]
	
	# NOTE
	# for Dpa-canu-het00c-YaHS-v202304 and for this comparison only,
	#   Dpa1 actually contains an inversion spanning from ~42 Mbps to ~43 Mbps
	#   for this analysis ONLY we will be using a Dpa1 with this inversion REVERTED (re-inverted)
	#   so that nucmer's delta-filter -g will consider it
	
	
	window_size = 100000
	sample_frequency = 50000
	reject_threshold = 0.7
	
	# run all the necessary alignment work using nucmer and mummer utilities
	# wrapper script to run a bunch of nucmers because I'm lazy
	# but also because this will standardize all the file names making parsing a lot easier
	for reference_scaffold in list_of_scaffolds:
		for query_scaffold in list_of_scaffolds:
			# pass if ref = query (not informative for this purpose)
			if reference_scaffold == query_scaffold:
				pass
			
			# otherwise, if ref is not the query
			else:
				# run nucmer
				pair = reference_scaffold + "_vs_" + query_scaffold
				nucmer_reference_fasta = reference_scaffold + ".fasta"
				nucmer_query_fasta = query_scaffold + ".fasta"
				
				command_to_use = "nucmer -t 10 " + nucmer_reference_fasta + " " + nucmer_query_fasta + " -p " + pair
				print("Now running: ", command_to_use)
				os.system(command_to_use)
				
				# run delta-filter
				command_to_use = "delta-filter -g " + pair + ".delta > " + pair + ".delta-filter_option-g.delta"
				print("Now running: ", command_to_use)
				os.system(command_to_use)
				
				# run show-coords
				command_to_use = "show-coords -rT " + pair + ".delta-filter_option-g.delta > " + pair + ".delta-filter_option-g.show-coords_option-r-T.coords"
				print("Now running: ", command_to_use)
				os.system(command_to_use)
				
				# run show-snps
				command_to_use = "show-snps -CrT " + pair + ".delta-filter_option-g.delta > " + pair + ".delta-filter_option-g.show-snps_option-C-r-T.snps"
				print("Now running: ", command_to_use)
				os.system(command_to_use)
				
				
	# now parse the alignments, coordinates, snps etc
	for reference_scaffold in list_of_scaffolds:
		for query_scaffold in list_of_scaffolds:
		
			# pass if ref == query (self-to-self nucmer, not run)
			if reference_scaffold == query_scaffold:
				pass
			
			# otherwise, if ref is not the query
			else:
				# define all the filenames from the nucmer (and utilities) runs first
				show_coords_output = reference_scaffold + "_vs_" + query_scaffold + ".delta-filter_option-g.show-coords_option-r-T.coords"
				show_snps_output = reference_scaffold + "_vs_" + query_scaffold + ".delta-filter_option-g.show-snps_option-C-r-T.snps"
				identity_array_filename = reference_scaffold + "_vs_" + query_scaffold + ".comparison.identity_array.txt"
				deletion_positions_filename = reference_scaffold + "_vs_" + query_scaffold + ".comparison.deletion_positions.txt"
				nucmer_reference_fasta = reference_scaffold + ".fasta"
				
				identity_array = generate_identity_array(show_coords_output, show_snps_output, nucmer_reference_fasta)
				
				# this next bit is optional but will make debugging a lot easier
				write_identity_array_to_files(identity_array, deletion_positions_filename, identity_array_filename)
				
				
				# now figure out percentage id (one ref vs one query) along the length of the scaffold
				# percent_id_output_file
				window_size = 100000
				sample_frequency = 50000
				reject_threshold = 0.7
				
				percent_id_output_file = reference_scaffold + "_vs_" + query_scaffold + ".comparison.percentage_id_along_length.window_" + str(window_size) + ".sample_freq_" + str(sample_frequency) + ".reject_threshold_" + str(reject_threshold) + ".txt"
				
				identity_array_to_percent_id(identity_array, window_size, sample_frequency, reject_threshold, percent_id_output_file)
				
	
	# now concatenate and format the percent_id_output_file s for R
	comparison_info = "comparison.percentage_id_along_length.window_" + str(window_size) + ".sample_freq_" + str(sample_frequency) + ".reject_threshold_" + str(reject_threshold) + ".txt"
	nucmer_percentage_id_along_length_v5_output_prepper_for_R(list_of_scaffolds, list_of_scaffolds, comparison_info)


	
	
if __name__ == '__main__':
    main()
    

