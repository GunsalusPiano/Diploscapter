# worm_telomere_lengths_extractor_v2.py
# extracts telomere lengths for publication
# makes the code a bit eaiser to adapt to other libraries


from Bio import SeqIO
import os, sys, getopt

version_number = "0.0.2"
TideHunter_path = "/mnt/c/Ubuntu_apps/TideHunter-v1.4.2_x64-linux/"
# fetch_reads_from_SRA(list_of_SRA_accessions)
# Make sure sra-tools v3.0.5 is installed
# download as FASTA (quality info not required)
def fetch_reads_from_SRA(list_of_SRA_accessions):
	for SRA_accession in list_of_SRA_accessions:
		string_to_write = "fasterq-dump --fasta " + SRA_accession
		print(string_to_write)
		os.system(string_to_write)

def eliminate_path(filepath_and_name):
    filename_bits = filepath_and_name.split("/")
    return filename_bits[-1]

def base_name(filename):
    filename_bits = filename.split(".")
    filename_bits.pop(-1)
    return ".".join(filename_bits)

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch

def consolidate_reads(list_of_read_files, consolidated_filename):
	string_to_write = "cat "
	for read_file in list_of_read_files:
		string_to_write = string_to_write + read_file + " "
	
	string_to_write = string_to_write + "> " + consolidated_filename
	os.system(string_to_write)

def split_read_file(large_read_file, split_factor):
	
	# count the number of reads in the file
	total_number_of_reads = 0
	
	for read in SeqIO.parse(large_read_file, "fasta"):
		total_number_of_reads = total_number_of_reads + 1
	
	print("total number of reads:", total_number_of_reads)
	
	number_of_reads_per_smaller_file = total_number_of_reads//(split_factor-1)
	
	list_of_reads = []
	count = 1
	output_filename = large_read_file.replace(".fasta", "") + "_1.fasta"
	for read in SeqIO.parse(large_read_file, "fasta"):
		list_of_reads.append(read)
		count = count + 1
		print("read:", count)
		if count % number_of_reads_per_smaller_file == 0:
			SeqIO.write(list_of_reads, output_filename, "fasta")
			list_of_reads = []
			output_filename = output_filename.replace(".fasta", "") + "_" + str(count/number_of_reads_per_smaller_file + 1) + ".fasta"
		else:
			pass
				
	SeqIO.write(list_of_reads, output_filename, "fasta")			

def cyclical_permute(sequence):
    return sequence[1:] + sequence[0:1]
	
def cyclical_permute_list(sequence):
    sequence = sequence.upper()
    return_matrix = [sequence]
    next_cyclical_permutation = cyclical_permute(sequence)
    while next_cyclical_permutation != sequence:
        return_matrix.append(next_cyclical_permutation)
        next_cyclical_permutation = cyclical_permute(return_matrix[len(return_matrix)-1])

    return return_matrix

def reverse_complement(arg):
    revcomp_dictionary = {"A": "T", \
                         "a": "T", \
                         "C": "G", \
                         "c": "G", \
                         "G": "C", \
                         "g": "C", \
                         "T": "A", \
                         "t": "A"}
    sequence = str(arg)
    reverse_complement = ""
    for letter in sequence[::-1]:
        reverse_complement = reverse_complement + revcomp_dictionary[letter]
    
    return reverse_complement

# penalty_reward_array = [-1, -1, -1, -1, -1, -1, .... 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
def weighted_score_of_repeated_sequence_in_read(repeat_occupancy_array, penalty, penalty_zone_size, reward, reward_zone_size):
	length_of_repeat_occupancy_array = len(repeat_occupancy_array)
	
	if length_of_repeat_occupancy_array >= penalty_zone_size + reward_zone_size:
		penalty_reward_array = [penalty]*penalty_zone_size + [0]*(length_of_repeat_occupancy_array-penalty_zone_size-reward_zone_size) + [reward]*reward_zone_size
		
		dot_product = 0
		
		for array_index in range(length_of_repeat_occupancy_array):
			dot_product = repeat_occupancy_array[array_index]*penalty_reward_array[array_index] + dot_product
		
		return dot_product
	else:
		return -1000000

def TideHunter_telomere_score(list_of_telomere_sequences, list_of_tel_revcomp_sequences, TideHunter_output, penalty, penalty_zone_size, reward, reward_zone_size, TideHunter_telomere_score_output):
	with open(TideHunter_output, "r") as TideHunter_output_fh, open(TideHunter_telomere_score_output, "w") as putative_telomere_TideHunter_reads_fh:
		line = TideHunter_output_fh.readline()
		old_read_name = ""
		repeat_occupancy_array = []
		coverage = 0
		read_length = 0
		repeat_start = 0
		repeat_end = 0
		telomere_read_score = 0
		revcomp = ""
		
		putative_telomere_TideHunter_reads_fh.write("read_name\tlikely_orientation\tread_length\tcoverage\ttelomere_read_score\n")
		
		while line:
			line_items = line.split()
			
			repeat_sequence = line_items[10]
			read_length = int(line_items[2])
			
			if read_length >= penalty_zone_size + reward_zone_size and (repeat_sequence in list_of_telomere_sequences or repeat_sequence in list_of_tel_revcomp_sequences):
				
				new_read_name = line_items[0]
				
				if new_read_name != old_read_name:
				
					if old_read_name != "":
						coverage = sum([abs(occupancy) for occupancy in repeat_occupancy_array])
						telomere_read_score = weighted_score_of_repeated_sequence_in_read(repeat_occupancy_array, penalty, penalty_zone_size, reward, reward_zone_size)
						if sum(repeat_occupancy_array) < 0 and telomere_read_score > 0:
							revcomp = "revcomp"
						else:
							revcomp = "as-is"
						string_to_write = old_read_name + "\t" + revcomp + "\t" +str(read_length) + "\t" + str(coverage) + "\t" + str(telomere_read_score) + "\n"
						putative_telomere_TideHunter_reads_fh.write(string_to_write)
						#putative_telomere_TideHunter_reads_fh.write(str(repeat_occupancy_array))
					
					else:
						pass
					
					read_length = int(line_items[2])
					repeat_start = int(line_items[3])
					repeat_end = int(line_items[4])
					
					repeat_occupancy_array = [0]*read_length
					
					# record the occupancy of the repeats
					for i in range(repeat_start - 1, repeat_end):
						if repeat_sequence in list_of_telomere_sequences:
							repeat_occupancy_array[i] = repeat_occupancy_array[i] + 1
						elif repeat_sequence in list_of_tel_revcomp_sequences:
							repeat_occupancy_array[i] = repeat_occupancy_array[i] - 1
						else:
							pass
					
					# set old_read_name to new_read_name
					old_read_name = new_read_name
					
				else:
					repeat_start = int(line_items[3])
					repeat_end = int(line_items[4])
					for i in range(repeat_start - 1, repeat_end):
						if repeat_sequence in list_of_telomere_sequences:
							repeat_occupancy_array[i] = repeat_occupancy_array[i] + 1
						elif repeat_sequence in list_of_tel_revcomp_sequences:
							repeat_occupancy_array[i] = repeat_occupancy_array[i] - 1
						else:
							pass
			line = TideHunter_output_fh.readline()
		
		coverage = sum([abs(occupancy) for occupancy in repeat_occupancy_array])
		telomere_read_score = weighted_score_of_repeated_sequence_in_read(repeat_occupancy_array, penalty, penalty_zone_size, reward, reward_zone_size)
		if sum(repeat_occupancy_array) < 0 and telomere_read_score > 0:
			revcomp = "revcomp"
		else:
			revcomp = "as-is"
		string_to_write = old_read_name + "\t" + revcomp + "\t" + str(read_length) + "\t" + str(coverage) + "\t" + str(telomere_read_score) + "\n"
		putative_telomere_TideHunter_reads_fh.write(string_to_write)
		#putative_telomere_TideHunter_reads_fh.write(str(repeat_occupancy_array))

def TideHunter_telomere_score_greater_than(TideHunter_telomere_score_cutoff, TideHunter_telomere_score_output, TideHunter_telomere_score_greater_than_output):
	with open(TideHunter_telomere_score_output, "r") as TideHunter_telomere_score_output_fh, open(TideHunter_telomere_score_greater_than_output, "w") as TideHunter_telomere_score_greater_than_output_fh:
		# transfer the header line
		line = TideHunter_telomere_score_output_fh.readline()
		TideHunter_telomere_score_greater_than_output_fh.write(line)
		
		# now evaluate: only transfer lines over to the new file if score is > TideHunter_telomere_score_cutoff
		line = TideHunter_telomere_score_output_fh.readline()
		
		while line:
			line_items = line.split()
			telomere_read_score = int(line_items[4].rstrip())
			
			if telomere_read_score > TideHunter_telomere_score_cutoff:
				TideHunter_telomere_score_greater_than_output_fh.write(line)
			else:
				pass
			
			line = TideHunter_telomere_score_output_fh.readline()


# get reads listed in the first column of a 'file' starting from the second line from 'reads_fasta_filename'
'''def get_reads(file, reads_fasta_filename, new_reads_fasta_filename):
    seqrecord_dict = SeqIO.index(reads_fasta_filename, "fasta")
    with open(file, "r") as fh, open(new_reads_fasta_filename, "w") as new_reads_fasta_fh:
        fh.readline()
        line = fh.readline()
        
        while line:
            read_id = line.split()[0]
            SeqIO.write(seqrecord_dict[read_id], new_reads_fasta_fh, "fasta")
            line = fh.readline()
'''
def get_reads(file, reads_fasta_filename, new_reads_fasta_filename):
	with open(file, "r") as fh, open(new_reads_fasta_filename, "w") as new_reads_fasta_fh:
		fh.readline()
		line = fh.readline()
		for read in SeqIO.parse(reads_fasta_filename, "fasta"):
			if line:
				if read.id == line.split()[0]:
					SeqIO.write(read, new_reads_fasta_fh, "fasta")
					line = fh.readline()
				else:
					pass
			else:
				break

# be sure to do nucmer ref = subtelomere repeat, query = reads
def sort_reads_by_Dpachys_subtelomeres(nucmer_showcoords_output, sort_reads_by_Dpachys_subtelomeres_output):
	list_of_reads_to_keep = []
	
	with open(nucmer_showcoords_output, "r") as nucmer_showcoords_output_fh, open(sort_reads_by_Dpachys_subtelomeres_output, "w") as sort_reads_by_Dpachys_subtelomeres_output_fh:
		# skip the first 4 lines of nucmer show-coords
		for i in range(4):
			sort_reads_by_Dpachys_subtelomeres_output_fh.write(nucmer_showcoords_output_fh.readline())
		
		list_of_lines = []
		
		line = nucmer_showcoords_output_fh.readline()
		
		keep_read_condition = False
		
		old_read = line.split()[10].rstrip()
		
		
				
		while line:
			current_read = line.split()[10].rstrip()
			
			if current_read == old_read:
				# first, save the read_start position and the line  
				read_start = int(line.split()[2])
				list_of_lines.append([read_start, line])
				
			else:
				if keep_read_condition == True:
					# sort the list_of_lines
					# write to sort_reads_by_Dpachys_subtelomeres_output_fh
				
					list_of_lines.sort(key = lambda x: x[0])
					
					for item in list_of_lines:
						sort_reads_by_Dpachys_subtelomeres_output_fh.write(item[1])
				
				else:
					pass
				
				# now keep the new line and its info
				read_start = int(line.split()[2])
				list_of_lines = [[read_start, line]]
				
				# let keep_read_condition be re-evaluated
				keep_read_condition = False
				
				# re-set the old_read name
				old_read = current_read
				
			# evaluate if we should keep the read (ie. does it match roughly position 1 to the end of the subtelomeric repeat)
            # this block of code ensures AT LEAST one almost full repeat (ie. matches from beginning to end, with 5 bp wiggle room at each end)
			subtel_repeat_start = int(line.split()[0])
			subtel_repeat_end = int(line.split()[1])
			subtel_repeat_length = int(line.split()[7])
				
			if subtel_repeat_start < 5 and subtel_repeat_end > subtel_repeat_length - 5 and subtel_repeat_end <= subtel_repeat_length:
				keep_read_condition = keep_read_condition or True
			else:
				keep_read_condition = keep_read_condition or False
			
			line = nucmer_showcoords_output_fh.readline()
		
		# print any remaining reads in the list
		if keep_read_condition == True:
			# sort the list_of_lines
			# write to sort_reads_by_Dpachys_subtelomeres_output_fh
				
			list_of_lines.sort(key = lambda x: x[0])
					
			for item in list_of_lines:
				sort_reads_by_Dpachys_subtelomeres_output_fh.write(item[1])
		else:
			pass
			
			
			

# There is a problem: sometimes the subtel sequence hits the reads multiple times.
# Try to get the very last match (orientation -> telomere) or first match (orientation telomere <-) on the read if multiple matches
# This is for situations where the subtelomere has a repetitive structure (ie. Diploscapter)
# So if the structure of the subtelomere were
# [--subtel_1-->[--subtel_2-->[--subtel_3-->[--subtel_4-->[--subtel_5-->[--su_6>
# the this code keeps subtel_6, even if it's a partial match
def eliminate_multiple_nucmer_subtel_matches(nucmer_showcoords_output, eliminate_multiple_nucmer_subtel_matches_output):
	# now figure out the orientation of the read
	# then pick the TideHunter entry with the last match (--->telomere) or first match (telomere<----)
	with open(nucmer_showcoords_output, "r") as nucmer_showcoords_output_fh, open(eliminate_multiple_nucmer_subtel_matches_output, "w") as eliminate_multiple_nucmer_subtel_matches_output_fh:
		for i in range(4):
			line = nucmer_showcoords_output_fh.readline()
			eliminate_multiple_nucmer_subtel_matches_output_fh.write(line)
		
		line = nucmer_showcoords_output_fh.readline()
		old_read_name = ""
		old_min_query_start = 0
		old_max_query_end = 0
		old_line = ""
		
		while line:
			line_items = line.split()
			
			if line_items[10] != old_read_name:
				if old_read_name == "":
					#print("in block 1")
					pass
				else:
					eliminate_multiple_nucmer_subtel_matches_output_fh.write(old_line)
					#print(old_line)
				old_read_name = line_items[10]
				old_min_query_start = int(line_items[2])
				old_max_query_end = int(line_items[3])
				old_line = line
			else:
				# evaluate if its orientation
				# if left-to-right, see if the max query end needs to be updated
				if int(line_items[3]) > int(line_items[2]):
					#print("in block 2")
					if int(line_items[3]) > old_max_query_end:
						#print("in block 3")
						old_max_query_end = int(line_items[3]) 
						old_line = line
					else:
						#print("in block 4")
						pass
				
				# else if right-to-left, see if the MINIMUM query end needs to be updated
				else:
					if int(line_items[3]) < old_max_query_end:
						#print("in block 5")
						old_max_query_end = int(line_items[3])
						old_line = line
					else:
						#print("in block 6")
						pass
			line = nucmer_showcoords_output_fh.readline()
		
		eliminate_multiple_nucmer_subtel_matches_output_fh.write(old_line)
				
		


# finds the last nucmer match that can be confidently identified as 'subtelomere'
# This is for situations where the subtelomere is reasonably non-repetitive (eg. the last ~n kb of every assembled end is unique)
# Hence we're looking for full-length matches to the subtelomere, rather than spurious sequence matches
# So if the subtelomere had this structure:
# -----[-frag->-----------[--subtel-->[--tel-->
# the code will ignore the subtel-like frag
def find_nucmer_subtel_match(nucmer_showcoords_output, nucmer_leeway, find_nucmer_subtel_match_output):
	with open(nucmer_showcoords_output, "r") as nucmer_showcoords_output_fh, open(find_nucmer_subtel_match_output, "w") as find_nucmer_subtel_match_output_fh:
		# transfer the first 4 lines of nucmer show-coords output to the new file
		for i in range(4):
			line = nucmer_showcoords_output_fh.readline()
			find_nucmer_subtel_match_output_fh.write(line)
		
		line = nucmer_showcoords_output_fh.readline()
		
		while line:
			line_items = line.split()
			ref_end = int(line_items[1])
			subtel_length = int(line_items[7])
            
			# consider line only if the reference end (the end of the subtelomere) is matched.
			# this should be position 2000
			# however because the reads are noisy, we give it a 'nucmer_leeway'
			if ref_end in range(subtel_length - nucmer_leeway, subtel_length+1):
				ref_start = int(line_items[0])
				
				# if the entire subtelomere of subtel_length bps are matched (with leeway)
				# transfer the line to the new file
				if ref_start in range(nucmer_leeway + 1):
					find_nucmer_subtel_match_output_fh.write(line)
				# Otherwise
				else:
					# if the read starts in the middle of the 2000 bps of subtelomere
					query_start = int(line_items[2])
					query_end = int(line_items[3])
					
					query_length = int(line_items[8])
					
					# figure out the orientation (query_end - query_start)
					# if telomere on the right (query_end > query_start)
					# allow only if query_start is near the beginning of the read (matching to the subtelomere)
					
					if query_end > query_start and query_start in range(nucmer_leeway+1):
						find_nucmer_subtel_match_output_fh.write(line)
					
					# else if the telomere is on the left
					# allow only if the query end is almost at the right end of the read
					elif query_end < query_start and query_end in range(query_length - nucmer_leeway, query_length + 1):
						find_nucmer_subtel_match_output_fh.write(line)
					
					# otherwise, discard the line
					else:
						pass
			# if only part of the subtel is matched (and not the end), discard the line
			else:
				pass
			
			line = nucmer_showcoords_output_fh.readline()
            
# returns all positions where a substring is found on a string`
def find_all_substring_positions(arg_string, arg_substring):
    positions = []
    found = True
    current_position = 0
    while found:
        if arg_string.upper().find(arg_substring.upper(), current_position) >= 0:
            positions.append(arg_string.find(arg_substring.upper(), current_position))
            current_position += len(arg_substring)
        else:
            found = False
    return positions 

# find the largest dropoff of the telomere_text occupancy within sequence_text
# arg_text_telomere is a comma-delimited string of telomere variants (no space)
def find_maximum_dropoff(sequence_text, arg_text_telomere, scanning_window):
    # occupancy_array
    # 0 for every position not occupied by telomere motif
    # 1 for every position occupied by telomere motif
    occupancy_array = [0]*len(sequence_text)
    
    # dropoff_array
    # for every position:
    #   calculate the difference in the occupancy in ranges of [0:position] and [position:end]
    # First and Last elements are 0 because the occupancy cannot be divided by 0
    dropoff_array = [0]*len(sequence_text)
    
    cyclical_permutations = []
    for variant in arg_text_telomere.split(","):
        cyclical_permutations = cyclical_permutations + cyclical_permute_list(variant)
    
    for telomere_permutation in cyclical_permutations:
        found_positions = find_all_substring_positions(sequence_text, telomere_permutation)
        
        for position in found_positions:
            for i in range(len(telomere_permutation)):
                occupancy_array[position+i] = 1
    #print(occupancy_array)
    
    for position in range(scanning_window, len(dropoff_array)-scanning_window):
        occupancy_before_position = sum(occupancy_array[position - scanning_window:position])/scanning_window
        occupancy_after_position = sum(occupancy_array[position:position + scanning_window])/scanning_window
        
        dropoff_array[position] = abs(occupancy_before_position - occupancy_after_position)
    
    #print(dropoff_array)
    
    max_dropoff_position = dropoff_array.index(max(dropoff_array))
    max_dropoff = dropoff_array[max_dropoff_position]
    
    return max_dropoff_position + 1, max_dropoff
    


def locate_telomere_occupancy_drop(find_nucmer_subtel_match_output, new_reads_fasta_filename, scanning_window, arg_text_telomere, locate_telomere_occupancy_drop_output):
    #telomeric_repeat = telomeric_repeat.upper()
    #telomeric_repeat_revcomp = telomeric_repeat_revcomp.upper()
	
    # CONSTRUCT	A dictionary of reads first
    read_dictionary = {}
    for read in SeqIO.parse(new_reads_fasta_filename, "fasta"):
        read_dictionary[read.id] = read.seq.upper()
        
    with open(find_nucmer_subtel_match_output, "r") as find_nucmer_subtel_match_output_fh, open(locate_telomere_occupancy_drop_output, "w") as locate_telomere_occupancy_drop_output_fh:
        # skip the first 4 lines of find_nucmer_subtel_match_output
        for i in range(4):
            find_nucmer_subtel_match_output_fh.readline()
			
		
        # write the headers for locate_telomere_occupancy_drop_output
        locate_telomere_occupancy_drop_output_fh.write("read_name\tread_length\tsubtel_start\tsubtel_end\tmax_tel_dropoff\tmax_tel_dropoff_position\tmax_rev_dropoff\tmax_rev_dropoff_position\tdifference\n")
		
		# get read_name
        line = find_nucmer_subtel_match_output_fh.readline()
		
        
        while line:
            line_items = line.split()
            read_name = line_items[10].rstrip()
            read_length = int(line_items[8])
            read_sequence = read_dictionary[read_name]
			
            subtel_match_start = int(line_items[2])
            subtel_match_end = int(line_items[3])
            
            difference = 0
            #print(read_name)
            max_tel_repeat_dropoff_position, max_tel_repeat_dropoff = find_maximum_dropoff(read_sequence, arg_text_telomere, scanning_window)
            
            #print(max_tel_repeat_dropoff_position, max_tel_repeat_dropoff)
            #input()
            
            telomere_revcomp_variants = []
            for telomeric_motif in arg_text_telomere.split(","):
                telomere_revcomp_variants.append(reverse_complement(telomeric_motif))
            
            
            max_rev_repeat_dropoff_position, max_rev_repeat_dropoff = find_maximum_dropoff(read_sequence, ",".join(telomere_revcomp_variants), scanning_window)
            #print(max_rev_repeat_dropoff_position, max_rev_repeat_dropoff)
            #input()
            if subtel_match_end > subtel_match_start:
                difference = abs(subtel_match_end - max_tel_repeat_dropoff_position)
            else:
                difference = abs(subtel_match_end - max_rev_repeat_dropoff_position)
            
            line_to_write = read_name + "\t" + str(read_length) + "\t" + str(subtel_match_start) + "\t" + str(subtel_match_end) + "\t" + str(max_tel_repeat_dropoff) + \
                            "\t" + str(max_tel_repeat_dropoff_position) + "\t" + str(max_rev_repeat_dropoff) + "\t" + str(max_rev_repeat_dropoff_position) + "\t" + str(difference) + "\n"
            locate_telomere_occupancy_drop_output_fh.write(line_to_write)
            
            line = find_nucmer_subtel_match_output_fh.readline()	
				


def telomere_lengths_use_occupancy_dropout_coordinate(locate_telomere_occupancy_drop_output, position_difference_leeway, library_accession_number, species, strain, description, telomere_lengths_output, file_header):
    with open(locate_telomere_occupancy_drop_output, "r") as locate_telomere_occupancy_drop_output_fh, open(telomere_lengths_output, "w") as telomere_lengths_output_fh:
        # write header line of telomere_lengths_output
        if file_header == True:
            telomere_lengths_output_fh.write("accession\tspecies\tstrain\tdescription\tread_name\ttelomere_length\n")
        # skip the first line in locate_telomere_occupancy_drop_output
        locate_telomere_occupancy_drop_output_fh.readline()
        line = locate_telomere_occupancy_drop_output_fh.readline()
        
        while line:
            line_items = line.split()
            read_name = line_items[0]
            read_length = int(line_items[1])
            subtel_start = int(line_items[2])
            subtel_end = int(line_items[3])
            difference_in_position = int(line_items[8])
            max_tel_dropoff_position = int(line_items[5])
            max_rev_dropoff_position = int(line_items[7])
            telomere_length = 0
            
            if difference_in_position < position_difference_leeway:
                if subtel_end > subtel_start:
                    telomere_length = read_length - max_tel_dropoff_position
                else:
                    telomere_length = max_rev_dropoff_position - 1
                
                line_to_write = library_accession_number + "\t" + species + "\t" + strain + "\t" + description + "\t" + read_name + "\t" + str(telomere_length) + "\n"
                telomere_lengths_output_fh.write(line_to_write)
            else:
                pass
            line = locate_telomere_occupancy_drop_output_fh.readline()
		
def read_library_info(arg_library_info_file):
    library_info = {}
    with open(arg_library_info_file, "r") as info_fh:
        line = info_fh.readline()
        while line:
            attribute = line.rstrip().split("\t")[0]
            value = line.rstrip().split("\t")[1]
            library_info[attribute] = value
            
            line = info_fh.readline()
    return library_info

def main(argv):
    # to run this script, we need
    # a) a library of long reads (-l LIBRARY), or an accession number to download reads from (-a ACCESSSION)
    # b) a telomeric motif that can be identified by TideHunter (v1.5.5) (-t TELOMERIC_MOTIF)
    # c) the terminal n nucleotides to consider to look for telomeric motifs (aka reward_zone_size) (-n NUM_OF_NUCLEOTIDES)
    # d) a telomere score cutoff: a dot product ranking of how likely a read is a telomeric read (-c CUTOFF_TELOSCORE)
    # e) a fasta file containing all the known/well-resolved subtelomeres (-s SUBTELOMERE_FASTA)
    # f) a telomeric motif that can be used in a text-based seearch (rather than using TideHunter) (-T TELOMERIC_TEXT optional, if not used = -t)
    #
    # optionally, if the subtelomere is highly repetitive, indicate if we need to consider the last (most distal) subtelomeric repeat (-r)
    #
    # software required
    # sra_toolkit
    # TideHunter
    # nucmer (from mummer)
    
    arg_library = ""
    arg_accession = ""
    arg_tidehunter_telomere = ""
    arg_n = 0
    arg_cutoff_teloscore = 0
    arg_subtelomere = ""
    arg_text_telomere = ""
    scanning_window = 20
    position_difference_leeway = 1000
    arg_library_info_file = ""
    arg_repetitive = False
    
    opts, args = getopt.getopt(argv[1:], "hvl:a:t:n:c:s:w:d:T:i:r", ["help", "version", "library=", "accession=", "TideHunter_telomeric_motif=", "num_of_terminal_nucleotides=",
                                                                   "cutoff_teloscore=", "subtelomere_FASTA=", "window_to_scan_for_telomeres=", "position_Difference_leeway=",
                                                                   "telomeric_motif_text=", "library_info_file=", "repetitive_subtel"])
    
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(0)
        elif opt in ("-v", "--version"):
            print(version_number)
            sys.exit(0)
        elif opt in ("-l", "--library"):
            arg_library = arg
        elif opt in ("-a", "--accession"):
            arg_accession = arg
            arg_library = arg_Library + ".fasta"
        elif opt in ("-t", "--TideHunter_telomeric_motif"):
            arg_tidehunter_telomere = arg
        elif opt in ("-n", "--num_of_terminal_nucleotides"):
            arg_n = int(arg)
        elif opt in ("-c", "--cutoff_teloscore"):
            arg_cutoff_teloscore = int(arg)
        elif opt in ("-s", "--subtelomere_FASTA"):
            arg_subtelomere = arg
        elif opt in ("-w", "--window_to_scan_for_telomeres"):
            scanning_window = int(arg)
        elif opt in ("-d", "--position_Difference_leeway"):
            position_difference_leeway = int(arg)
        elif opt in ("-T", "--telomeric_motif_text"):
            arg_text_telomere = arg
        elif opt in ("-i", "--library_info_file"):
            arg_library_info_file = arg
        elif opt in ("-r", "--repetitive_subtel"):
            arg_repetitive = True
        
    
    if arg_accession != "":
        if os.path.exists(arg_accession + ".fasta"):
            os.system(f"fasterq-dump --fasta {arg_accession}")
    elif arg_library != "":
        arg_accession = base_name(eliminate_path(arg_library))
    else:
        print("no library indicated")
        exit(0)
    
    telomeric_motif_length = len(arg_tidehunter_telomere)
    
    TideHunter_output = f"{arg_accession}.TideHunter_{telomeric_motif_length}-mers.txt"
    info_header = f"[{arg_accession}]: "
    
    # first, run TideHunter v1.4.2
    if os.path.exists(TideHunter_output):
        pass
    else:
        string_to_write = f"TideHunter -p {telomeric_motif_length} -P {telomeric_motif_length} -t 15 -f 2 {arg_library} > {TideHunter_output}"
        print(f"{info_header}Now running {string_to_write}")
        os.system(f"{TideHunter_path}{string_to_write}")
	
	# next, assign each TideHunter analyzed read with a telomere score
    list_of_telomere_sequences = cyclical_permute_list(arg_tidehunter_telomere)
    list_of_tel_revcomp_sequences = cyclical_permute_list(reverse_complement(arg_tidehunter_telomere))
    
    penalty = -1
    penalty_zone_size = arg_n
    reward = 1
    reward_zone_size = arg_n
    TideHunter_telomere_score_output = arg_accession + ".telomereScore.txt"
    prefix = arg_accession + ".telomereScore" + str(arg_cutoff_teloscore)
    if os.path.exists(TideHunter_telomere_score_output):
        pass
    else:
        print(f"{info_header}Now running TideHunter_telomere_score().")
        TideHunter_telomere_score(list_of_telomere_sequences, list_of_tel_revcomp_sequences, TideHunter_output, penalty, penalty_zone_size, reward, reward_zone_size, TideHunter_telomere_score_output)
    
    # next, only work with those reads with scores of >20
    TideHunter_telomere_score_greater_than_output = prefix + ".txt"
    
    if os.path.exists(TideHunter_telomere_score_greater_than_output):
        pass
    else:
        print(f"{info_header}Now running TideHunter_telomere_score_greater_than().")
        TideHunter_telomere_score_greater_than(arg_cutoff_teloscore, TideHunter_telomere_score_output, TideHunter_telomere_score_greater_than_output)
    
    # grab the reads with scores >99
    get_reads_output = prefix + ".fasta"
    
    if os.path.exists(get_reads_output):
        pass
    else:
        print(info_header + "Now running get_reads().")
        get_reads(TideHunter_telomere_score_greater_than_output, arg_library, get_reads_output)
	
    # do some nucmers
    for subtelomere in SeqIO.parse(arg_subtelomere, "fasta"):
        '''
        nucmer_ref_names = ["nucmerSubtel2807", "nucmerSubtel3008", "nucmerSubtel3200"]
        nucmer_refs_dict = {"nucmerSubtel2807": "Dpa1L2L_subtel01_2807bp_endadj77.fasta", \
                            "nucmerSubtel3008": "Dpa2R_subtel02_3008bp_endadj77.fasta", \
                            "nucmerSubtel3200": "Dpa1R_subtel01_3200bp_endadj77.fasta"}
		'''
        delta_file = prefix + "." + subtelomere.id + ".delta"
        coords_file = prefix + "." + subtelomere.id + ".coords"
        #nucmer_reference = "tmp.nucmer.ref.fasta"
        SeqIO.write([subtelomere], f"{subtelomere.id}.fasta", "fasta")
        if os.path.exists(delta_file):
            pass
        else:
            string_to_write = f"nucmer -p {prefix}.{subtelomere.id} -l 10 -c 20 {subtelomere.id}.fasta {get_reads_output}"
            print(f"{info_header}Now running {string_to_write}")
            os.system(string_to_write)
		
        if os.path.exists(coords_file):
            pass
        else:
            string_to_write = f"show-coords -l -T {delta_file} > {coords_file}"
            print(f"{info_header}Now running {string_to_write}")
            os.system(string_to_write)
            
        os.system(f"rm {subtelomere.id}.fasta")
        
        # sort some .coords file
        sort_reads_by_Dpachys_subtelomeres_output = prefix + "." + subtelomere.id + ".sorted.coords"
        if os.path.exists(sort_reads_by_Dpachys_subtelomeres_output):
            pass
        else:
            print(f"{info_header}Now running sort_reads_by_Dpachys_subtelomeres().")
            sort_reads_by_Dpachys_subtelomeres(coords_file, sort_reads_by_Dpachys_subtelomeres_output)
        
        # find the most telomeric nucmer match
        eliminate_multiple_nucmer_subtel_matches_output = prefix + "." + subtelomere.id + ".sorted.multiple_hits_eliminated.coords"
        if os.path.exists(eliminate_multiple_nucmer_subtel_matches_output):
            pass
        else:
            if arg_repetitive:
                print(f"{info_header}Now running eliminate_multiple_nucmer_subtel_matches()")
                eliminate_multiple_nucmer_subtel_matches(sort_reads_by_Dpachys_subtelomeres_output, eliminate_multiple_nucmer_subtel_matches_output)
            else:
                print(f"{info_header}Now running find_nucmer_subtel_match()")
                nucmer_leeway = 5
                find_nucmer_subtel_match(sort_reads_by_Dpachys_subtelomeres_output, nucmer_leeway, eliminate_multiple_nucmer_subtel_matches_output)
            
        # find where the occupancy of the telomeric repeats drop
		# here we do not look for the 12-mer TAAGGGTAAGGC
		# we look for instead TAAGGG and TAAGGC
        locate_telomere_occupancy_drop_output = f"{prefix}.{subtelomere.id}.sorted.multiple_hits_eliminated.scanningWindow{scanning_window}.telOccupancyDrop.txt"
        
        if os.path.exists(locate_telomere_occupancy_drop_output):
            pass
        else:
            if arg_text_telomere == "":
                arg_text_telomere = arg_tidehunter_telomere
            print(f"{info_header}Now running locate_telomere_occupancy_drop()")
            #input()
            locate_telomere_occupancy_drop(eliminate_multiple_nucmer_subtel_matches_output, get_reads_output, scanning_window, arg_text_telomere, locate_telomere_occupancy_drop_output)
        
        # the subtelomere ends can vary wildly in Diploscapter (in where they connect to the bona fide telomere)
        # hence, we will set the position_difference_leeway pretty high (ie. the telomere may start up to 1000 bps away from the last matched subtelomere)
        
        library_info = {}
        species = ""
        strain = ""
        description = ""
        
        if arg_library_info_file != "":
            library_info = read_library_info(arg_library_info_file)
            
        if arg_accession == "":
            if "accession" in library_info:
                arg_accession = library_info["accession"]
        
        if "species" in library_info:
            species = library_info["species"]
		
        if "strain" in library_info:
            strain = library_info["strain"]
        
        if "description" in library_info:
            description = library_info["description"]
        
        telomere_lengths_output = f"{prefix}.{subtelomere.id}.sorted.multiple_hits_eliminated.scanningWindow{scanning_window}.positionLeeway{position_difference_leeway}.telomereLengths."
        if os.path.exists(telomere_lengths_output):
            pass
        else:
            print(info_header + "Now running telomere_lengths_use_occupancy_dropout_coordinate()")
            telomere_lengths_use_occupancy_dropout_coordinate(locate_telomere_occupancy_drop_output, position_difference_leeway, arg_accession, species, strain, description, telomere_lengths_output + "txt", True)
            telomere_lengths_use_occupancy_dropout_coordinate(locate_telomere_occupancy_drop_output, position_difference_leeway, arg_accession, species, strain, description, telomere_lengths_output + "headerless.txt", False)
    
    all_telomeres_summary = f"{arg_accession}.positionLeeway{position_difference_leeway}.allTelomereLengths.headerless.txt"
    if os.path.exists(all_telomeres_summary):
        pass
    else:
        os.system(f"cat *.positionLeeway{position_difference_leeway}.telomereLengths.headerless.txt > {all_telomeres_summary}")
    
if __name__ == "__main__":
    main(sys.argv)