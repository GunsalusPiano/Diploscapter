################################################################
# percent_GC_RepeatMasker_masked_repeats_2023-10-19.py
# Calculates the repeats content (as identifed by earlGrey)
# along the lengths of sequences in a fasta file in which
# earlGrey did the repeat analysis.
#
# REQUIRES:
#    FASTA sequence file of multiple chromosomes
#    GFF output file of the same FASTA from RepeatMasker
#    GTF output file of BRAKER2/Augustus
# Allows the user to adjust the window size and increment
# 
################################################################


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from statistics import mean
import os
import getopt
import sys
#import math

def index_formatter(file_index):
    if file_index > 0 and file_index < 10:
        return "00" + str(file_index)
    elif file_index >=10 and file_index < 100:
        return "0" + str(file_index)
    else:
        return str(file_index)

# filename_extension_remover(filename)
# removes the extension (last string after .) of a file name
def filename_extension_remover(filename):
	filename_split = filename.split(".")
	
	filename_without_extension = ""
	
	for i in range(len(filename_split)-1):
		filename_without_extension = filename_without_extension + filename_split[i] + "."
	
	return filename_without_extension[0:len(filename_without_extension)-1]

# quick_mask(fasta, masking_gff)
# uses the masking gff at the end of the earlGrey pipeline to mask the original genome
def quick_mask(fasta_file, masking_gff):
    new_fasta_file = filename_extension_remover(fasta_file) + ".masked.fasta"
    
    new_sequence_text = ""
    
    with open(new_fasta_file, "w") as new_fasta_fh:
        for old_sequence in SeqIO.parse(fasta_file, "fasta"):
            # make masked_position_array
            # every position is 1 where not masked
            # 0 where masked
            masked_position_array = [1]*(len(old_sequence)+1)
            
            with open(masking_gff, "r") as masking_gff_fh:
                line = masking_gff_fh.readline()
                
                while line:
                    sequence_id = line.split("\t")[0]
                    
                    if sequence_id == old_sequence.id:
                        start = int(line.split("\t")[3])
                        end = int(line.split("\t")[4])
                        
                        for position in range(start, end+1):
                            masked_position_array[position] = 0
                    else:
                        pass
                    
                    line = masking_gff_fh.readline()
            
            masked_position_array.pop(0)
            new_sequence_text = ""
            
            for i in range(len(old_sequence)):
                if masked_position_array[i] == 0:
                    new_sequence_text = new_sequence_text + "N"
                else:
                    new_sequence_text = new_sequence_text + old_sequence[i]
            SeqIO.write(SeqRecord(Seq(new_sequence_text), id=old_sequence.id, description=old_sequence.description, name=old_sequence.name), new_fasta_fh, "fasta")

# make sure the bash script, rmOutToGFF3custom.bash, by Daren Card (and modified by me) is in the folder    
def turn_RepeatMasker_out_to_gff(RepeatMasker_out_filename):
    line_to_write = "./rmOutToGFF3custom.bash -o " + RepeatMasker_out_filename + " > " + filename_extension_remover(RepeatMasker_out_filename) + ".gff"
    print(line_to_write)
    os.system(line_to_write)

def replace_earlGrey_ctg_names(earlGrey_RepeatMasker_gff_filename, fasta_filename):
    ctg_dict = {}
    ctg_index = 1
    for seq in SeqIO.parse(fasta_filename, "fasta"):
        ctg_name = "ctg_" + str(ctg_index)
        ctg_dict[ctg_name] = seq.id
        ctg_index = ctg_index + 1
    
    for ctg_name in ctg_dict:
        line_to_write = "sed -i 's/"+ ctg_name + "/" + ctg_dict[ctg_name] + "/g' " + earlGrey_RepeatMasker_gff_filename
        os.system(line_to_write)
'''    
def RepeatMasker_feature_tally(GFF_filename):
    features_list_filename = filename_extension_remover(GFF_filename) + ".featuresList.txt"
    
    repeats_dict = {}
    occupancy_dict = {}
    
    with open(GFF_filename, "r") as GFF_fh, open(features_list_filename, "w") as features_list_fh:
        line = GFF_fh.readline()
        
        while line:
            if line[0] == "#":
                pass
            else:
                attribute = line.split("\t")[8]
                # attribute has the pattern:
                # Target=Simple_repeat/(CCTTAC)n 1 87
                # we can get rid of "Target=", then split the string by space.
                # then take the first element
                repeat_type = attribute[7:].split(" ")[0].rstrip()
                start = int(line.split("\t")[3])
                end = int(line.split("\t")[4])
                
                if repeat_type in repeats_dict:
                    repeats_dict[repeat_type] = repeats_dict[repeat_type] + 1
                    occupancy_dict[repeat_type] = occupancy_dict[repeat_type] + end - start + 1
                else:
                    repeats_dict[repeat_type] = 1
                    occupancy_dict[repeat_type] = end - start + 1
            line = GFF_fh.readline()
            
        features_list_fh.write("repeat_type\toccurrence\toccupancy\n")
        for repeat_type in repeats_dict:
            line = repeat_type + "\t" + str(repeats_dict[repeat_type]) + "\t" + str(occupancy_dict[repeat_type]) + "\n"
            features_list_fh.write(line)
'''

def extract_gff_attribute(gff_line, attribute_type):
    attributes = gff_line.split("\t")[8].rstrip()
    
    attribute_dict = {}
    for attribute in attributes.split(";"):
        attribute_dict[attribute.split("=")[0]] = attribute.split("=")[1]
    
    if attribute_type in attribute_dict:
        return attribute_dict[attribute_type]
    else:
        return "attribute_NA"

def combined_gff_feature_tally(GFF_filename, feature_tally_filename):
    #features_list_filename = filename_extension_remover(GFF_filename) + ".featuresList.txt"
    
    occurrence_dict = {}
    occupancy_dict = {}
    
    with open(GFF_filename, "r") as GFF_fh, open(feature_tally_filename, "w") as features_list_fh:
        line = GFF_fh.readline()
        
        while line:
            if line[0] == "#":
                pass
            else:
                repeat_type = line.split("\t")[2]+ "#" + extract_gff_attribute(line, "ID")
                start = int(line.split("\t")[3])
                end = int(line.split("\t")[4])
                
                if repeat_type in occurrence_dict:
                    occurrence_dict[repeat_type] = occurrence_dict[repeat_type] + 1
                    occupancy_dict[repeat_type] = occupancy_dict[repeat_type] + end - start + 1
                else:
                    occurrence_dict[repeat_type] = 1
                    occupancy_dict[repeat_type] = end - start + 1
            line = GFF_fh.readline()
            
        features_list_fh.write("repeat_type\toccurrence\toccupancy\n")
        for repeat_type in occurrence_dict:
            line = repeat_type + "\t" + str(occurrence_dict[repeat_type]) + "\t" + str(occupancy_dict[repeat_type]) + "\n"
            features_list_fh.write(line)


'''
# count_RepeatMasker_masked_content
#   gives a fraction of RepeatMasker-annotated repeats in a certain window (window_size)
#   this also can now pick out specific repeat types to count (listed in type_of_repeat_list_file)
#     type_of_repeat_list_file will have a header line beginning with ">" that describes the CLASS of repeats (repeat classification)
#     followed by types of repeats in that CLASS, separated by return characters.
def count_RepeatMasker_masked_content(seq, RepeatMasker_gtf_output, type_of_repeat_list_file, window_size, increment, xaxis_label, file_index):
    # extract the types of repeats in the type_of_repeat_list_file
    # repeat_dict has the dictionary structure:
    # key:value = repeat_type:repeat_classification
    repeat_dict = {}
    repeat_classification = ""
    
    # repeat_occupancy_dict has this structure:
    #   key:value = repeat_type:repeat_occupancy_array
    #   repeat_occupancy_array should have length_of_seq + 1 elements (to make it index_1)
    repeat_occupancy_dict = {}
    length_of_seq = len(seq)
    
    with open(type_of_repeat_list_file, "r") as type_of_repeat_list_fh:
        line = type_of_repeat_list_fh.readline()
        
        while line:
            if line[0] == ">":
                repeat_classification = line.rstrip().replace(">", "")
                repeat_occupancy_dict[repeat_classification] = [0]*(length_of_seq + 1)
            elif line == "\n":
                pass
            else:
                repeat_dict[line.rstrip()] = repeat_classification
                
            line = type_of_repeat_list_fh.readline()
       
    
    #print(repeat_dict)
    # repeat_occupancy_dict[repeat_classification]
    # records every position in DNA sequence that has a repeat annotated
    # 1 for a repeat annotated
    # 0 for no repeat annotated
    # make this a index_1 array (rather than index_0) because all the sequence annotations are index_1
    # therefore the 0th element doesn't have a corresponding nucleotide associated
    
    with open(RepeatMasker_gtf_output, "r") as RepeatMasker_gtf_output_fh:
        
        # read RepeatMasker output line by line
        line = RepeatMasker_gtf_output_fh.readline()
        
        while line:
            # ignore lines starting with #
            
            if line[0] == "#":
                pass
            else:
                line_items = line.split()
                scaffold = line_items[0]
                repeat_type = line_items[8].split(" ")[0].rstrip().replace("Target=", "")
                
                # only consider lines annotating the current seq (scaffold)
                # ignore all the rest
                if scaffold == seq.id and repeat_dict[repeat_type] in repeat_occupancy_dict:
                    start = int(line_items[3])
                    end = int(line_items[4])
                
                    # we want to turn the 0s in positions start to end to turn into 1s
                    # start to end INCLUSIVE (python range() is exclusive on the right end)
                    for x in range(start, end+1):
                        repeat_occupancy_dict[repeat_dict[repeat_type]][x] = 1
                else:
                    pass
            line = RepeatMasker_gtf_output_fh.readline()

    # start recording down % repeat content in the window, then increment through the length of the scaffold
    # because strings (how we're dealing with sequence data here) are index_0, let's turn the repeat_occupancy_array back to index_0
    
    for repeat_classification in repeat_occupancy_dict:
        repeat_occupancy_dict[repeat_classification].pop(0)


    start = 0
    end = 0
    
    if window_size > length_of_seq:
        end = length_of_seq
    else:
        end = window_size       
    
    temp_output_filename = ""
    temp_tidy_output_filename = ""
    
    if file_index == "":
        temp_output_filename = "tmp." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".txt"
        temp_tidy_output_filename = "tmp." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".tidy.txt"
    else:
        temp_output_filename = index_formatter(file_index) + "." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".txt"
        temp_tidy_output_filename = index_formatter(file_index) + "." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".tidy.txt"
    
    with open(temp_output_filename, "w") as temp_output_fh, open(temp_tidy_output_filename, "w") as temp_tidy_output_fh:
        temp_output_fh.write(xaxis_label)
            
        for repeat_classification in repeat_occupancy_dict:
            temp_output_fh.write("\t"+repeat_classification)
        
        temp_output_fh.write("\n")
        
        temp_tidy_output_fh.write(xaxis_label + "\trepeat_fraction\trepeat_class\n")
        
        while start < length_of_seq:
            # First, write down the chromosomal coordinate
            temp_output_fh.write(str((start + end)//2))
            
            
            # Write down the repeated fraction for each repeat_classification
            for repeat_classification in repeat_occupancy_dict:
                temp_tidy_output_fh.write(str((start + end)//2))
                
                # count the number of nucleotides that are annotated as repeats
                number_of_repeated_nucleotides = sum(repeat_occupancy_dict[repeat_classification][start:end])
            
                # count the number of nucleotides that are "N"
                number_of_N = seq.seq.upper().count("N", start, end)
                
                #line_to_print = "start="+str(start) + " end=" + str(end) + " N=" + str(number_of_N) + "\n"
                #print(line_to_print)
                
            
                # write these two values to the temp output file
                if end - start - number_of_N == 0:
                    temp_output_fh.write("\t\t")
                    temp_tidy_output_fh.write("\n")
                else:
                    percentage_repeated = number_of_repeated_nucleotides/(end - start - number_of_N)
                    temp_output_fh.write("\t" + str(percentage_repeated))
                    temp_tidy_output_fh.write("\t" + str(percentage_repeated) + "\t" + repeat_classification + "\n")
            
            temp_output_fh.write("\n")
            
            start = start + increment
        
            if end + increment > length_of_seq:
                end = length_of_seq
            else:
                end = end + increment
    
    #return percent_repeat_list
'''

def combine_gff(fasta_file, gff1, gff2, new_gff, subtel_repeat_names_sep_by_commas):
    scaffold_dict = {}
    for seqRecord in SeqIO.parse(fasta_file, "fasta"):
        if seqRecord.id in scaffold_dict:
            pass
        else:
            scaffold_dict[seqRecord.id] = []
    
    with open(gff1, "r") as gff1_fh, open(gff2, "r") as gff2_fh:
        line = gff1_fh.readline()
        while line:
            if line[0] == "#":
                pass
            else:
                scaffold = line.split("\t")[0]
                start = int(line.split("\t")[3])
                
                scaffold_dict[scaffold].append([start, line])
            line = gff1_fh.readline()
            
        line = gff2_fh.readline()
        while line:
            if line[0] == "#":
                pass
            else:
                scaffold = line.split("\t")[0]
                start = int(line.split("\t")[3])
                
                scaffold_dict[scaffold].append([start, line])
            line = gff2_fh.readline()
    
    subtel_repeat_names_list = subtel_repeat_names_sep_by_commas.split(",")
    
    with open(new_gff, "w") as new_gff_fh:
        for scaffold in scaffold_dict:
            if len(scaffold_dict[scaffold]) == 0:
                pass
                
            else:
                scaffold_dict[scaffold].sort(key = lambda x: x[0])
                
                #print(subtel_repeat_names_list)
                for line in scaffold_dict[scaffold]:
                    #print(extract_gff_attribute(line[1], "ID"))
                    #input()
                    if subtel_repeat_names_sep_by_commas != "" and extract_gff_attribute(line[1], "ID") in subtel_repeat_names_list:
                        line_items = line[1].split("\t")
                        line_items[2] = "Subtelomere"
                        "\t".join(line_items)
                        new_gff_fh.write("\t".join(line_items))
                    else:
                        new_gff_fh.write(line[1])
                        
                    

# uses an old GFF file (old_gff) and featuresList file (features_list_file),
# retains only lines that are "simple_repeats" (any number of times) or contain RepeatScout repeats ("Unspecified" on the GFF)
# that occur a number of times (occurrence_threshold) or more.
def discard_rare_repeats(old_gff, features_list_file, occurrence_threshold, new_gff):
    # make a dictionary to look up how many times a repeat occurs
    
    repeat_occurrence_dict = {}
    
    with open(features_list_file, "r") as features_list_fh:
        # skip the header line
        line = features_list_fh.readline()
    
        line = features_list_fh.readline()
        
        while line:
            repeat_name = line.split("\t")[0].split("#")[1]
            occurrence = int(line.split("\t")[1])
            repeat_occurrence_dict[repeat_name] = occurrence
            
            line = features_list_fh.readline()
            
    with open(old_gff, "r") as old_gff_fh, open(new_gff, "w") as new_gff_fh:
        line = old_gff_fh.readline()
        
        while line:
            if line[0] == "#":
                pass
            else:            
                repeat_class = line.split("\t")[2]
                repeat_name = extract_gff_attribute(line, "ID")
                
                
                if repeat_class == "Unspecified" and repeat_occurrence_dict[repeat_name] >= int(occurrence_threshold):
                    new_gff_fh.write(line)
                elif repeat_class == "Unspecified" and repeat_occurrence_dict[repeat_name] < int(occurrence_threshold):
                    pass
                else:
                    new_gff_fh.write(line)
            
            line = old_gff_fh.readline()
            
            
def auto_repeat_classification(features_list_file, repeat_classification_file):
    features_dict = {}
    
    with open(features_list_file, "r") as features_list_fh:
        # read the header line
        features_list_fh.readline()
        
        line = features_list_fh.readline()
        
        while line:
            # repeat_class is the text before the #
            # repeat_name is the text after the #
            repeat_class = line.split("\t")[0].split("#")[0]
            repeat_name = line.split("\t")[0].split("#")[1]
            
            # occurrence_threshold - only record down repeat classes which are "Unspecified" and occur more than 5x
            occurrence = int(line.split("\t")[1])
            
            if repeat_class in features_dict:
                features_dict[repeat_class].append(repeat_name)
            else:
                features_dict[repeat_class] = [repeat_name]
            
            line = features_list_fh.readline()
                
    with open(repeat_classification_file, "w") as repeat_classification_fh:
        for features_dict_key in features_dict:
            repeat_classification_fh.write(">"+features_dict_key+"\n")
            
            for repeat_name in features_dict[features_dict_key]:
                repeat_classification_fh.write(repeat_name+"\n")
                
def count_RepeatMasker_masked_content_v20231026(seq, RepeatMasker_gtf_output, type_of_repeat_list_file, window_size, increment, xaxis_label, file_prefix):
    # extract the types of repeats in the type_of_repeat_list_file
    # repeat_dict has the dictionary structure:
    # key:value = repeat_type:repeat_classification
    repeat_dict = {}
    repeat_classification = ""
    
    # repeat_occupancy_dict has this structure:
    #   key:value = repeat_type:repeat_occupancy_array
    #   repeat_occupancy_array should have length_of_seq + 1 elements (to make it index_1)
    repeat_occupancy_dict = {}
    length_of_seq = len(seq)
    
    with open(type_of_repeat_list_file, "r") as type_of_repeat_list_fh:
        line = type_of_repeat_list_fh.readline()
        
        while line:
            if line[0] == ">":
                repeat_classification = line.rstrip().replace(">", "")
                repeat_occupancy_dict[repeat_classification] = [0]*(length_of_seq + 1)
            elif line == "\n":
                pass
            else:
                repeat_dict[line.rstrip()] = repeat_classification
                
            line = type_of_repeat_list_fh.readline()
       
    
    #print(repeat_dict)
    # repeat_occupancy_dict[repeat_classification]
    # records every position in DNA sequence that has a repeat annotated
    # 1 for a repeat annotated
    # 0 for no repeat annotated
    # make this a index_1 array (rather than index_0) because all the sequence annotations are index_1
    # therefore the 0th element doesn't have a corresponding nucleotide associated
    
    with open(RepeatMasker_gtf_output, "r") as RepeatMasker_gtf_output_fh:
        
        # read RepeatMasker output line by line
        line = RepeatMasker_gtf_output_fh.readline()
        
        while line:
            # ignore lines starting with #
            
            if line[0] == "#":
                pass
            else:
                line_items = line.split()
                scaffold = line_items[0]
                repeat_type = extract_gff_attribute(line, "ID")
                
                # only consider lines annotating the current seq (scaffold)
                # ignore all the rest
                if scaffold == seq.id and repeat_dict[repeat_type] in repeat_occupancy_dict:
                    start = int(line_items[3])
                    end = int(line_items[4])
                
                    # we want to turn the 0s in positions start to end to turn into 1s
                    # start to end INCLUSIVE (python range() is exclusive on the right end)
                    for x in range(start, end+1):
                        repeat_occupancy_dict[repeat_dict[repeat_type]][x] = 1
                else:
                    pass
            line = RepeatMasker_gtf_output_fh.readline()

    # start recording down % repeat content in the window, then increment through the length of the scaffold
    # because strings (how we're dealing with sequence data here) are index_0, let's turn the repeat_occupancy_array back to index_0
    
    for repeat_classification in repeat_occupancy_dict:
        repeat_occupancy_dict[repeat_classification].pop(0)


    start = 0
    end = 0
    
    if window_size > length_of_seq:
        end = length_of_seq
    else:
        end = window_size       
    
    temp_output_filename = ""
    temp_tidy_output_filename = ""
    
    if file_prefix == "":
        temp_output_filename = "tmp." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".txt"
        temp_tidy_output_filename = "tmp." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".tidy.txt"
    else:
        temp_output_filename = file_prefix + "." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".txt"
        temp_tidy_output_filename = file_prefix + "." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".tidy.txt"
    
    with open(temp_output_filename, "w") as temp_output_fh, open(temp_tidy_output_filename, "w") as temp_tidy_output_fh:
        temp_output_fh.write(xaxis_label)
            
        for repeat_classification in repeat_occupancy_dict:
            temp_output_fh.write("\t"+repeat_classification)
        
        temp_output_fh.write("\tsum\n")
        
        temp_tidy_output_fh.write(xaxis_label + "\trepeat_fraction\trepeat_class\n")
        
        while start < length_of_seq:
            # First, write down the chromosomal coordinate
            temp_output_fh.write(str((start + end)//2))
            
            # set the sum of repeat occupancy percentages to 0
            sum_of_repeat_occupancy = 0
            
            # set the number of Ns in an interval to 0
            number_of_N = 0
            
            # Write down the repeated fraction for each repeat_classification
            for repeat_classification in repeat_occupancy_dict:
                # count the number of nucleotides that are annotated as repeats
                number_of_repeated_nucleotides = sum(repeat_occupancy_dict[repeat_classification][start:end])
                sum_of_repeat_occupancy = sum_of_repeat_occupancy + number_of_repeated_nucleotides
            
                # count the number of nucleotides that are "N"
                number_of_N = seq.seq.upper().count("N", start, end)
                
                #line_to_print = "start="+str(start) + " end=" + str(end) + " N=" + str(number_of_N) + "\n"
                #print(line_to_print)
                
            
                # write these two values to the temp output file
                if end - start - number_of_N == 0:
                    temp_output_fh.write("\t")
                    #temp_tidy_output_fh.write("\n")
                else:
                    percentage_repeated = number_of_repeated_nucleotides/(end - start - number_of_N)
                    #sum_of_repeat_occupancy = sum_of_repeat_occupancy + percentage_repeated
                    temp_output_fh.write("\t" + str(percentage_repeated))
                    temp_tidy_output_fh.write(str((start + end)//2) + "\t" + str(percentage_repeated) + "\t" + repeat_classification + "\n")
            if end - start - number_of_N == 0:
                temp_output_fh.write("\n")
            else:
                percentage_repeated = sum_of_repeat_occupancy/(end - start - number_of_N)
                temp_output_fh.write("\t"+ str(percentage_repeated) + "\n")
                temp_tidy_output_fh.write(str((start + end)//2) + "\t" + str(percentage_repeated) + "\t" + "sum\n")
            
            start = start + increment
        
            if end + increment > length_of_seq:
                end = length_of_seq
            else:
                end = end + increment                

def count_all_RepeatMasker_masked_content(seq, RepeatMasker_gtf_output, window_size, increment, xaxis_label, file_index):
    # repeat_occupancy_array
    # records every position in DNA sequence that has a repeat annotated
    # 1 for a repeat annotated
    # 0 for no repeat annotated
    # make this a index_1 array (rather than index_0) because all the sequence annotations are index_1
    # therefore the 0th element doesn't have a corresponding nucleotide associated
    
    repeat_occupancy_array = [0]*(len(seq) + 1)
    with open(RepeatMasker_gtf_output, "r") as RepeatMasker_gtf_output_fh:
        
        # read RepeatMasker output line by line
        line = RepeatMasker_gtf_output_fh.readline()
        
        while line:
            # ignore lines starting with #
            
            if line[0] == "#":
                pass
            else:
                line_items = line.split()
                scaffold = line_items[0]
                
                # only consider lines annotating the current seq (scaffold)
                # ignore all the rest
                if scaffold == seq.id:
                    start = int(line_items[3])
                    end = int(line_items[4])
                
                    # we want to turn the 0s in positions start to end to turn into 1s
                    # start to end INCLUSIVE (python range() is exclusive on the right end)
                    for x in range(start, end+1):
                        repeat_occupancy_array[x] = 1
                else:
                    pass
            line = RepeatMasker_gtf_output_fh.readline()

    # start recording down % repeat content in the window, then increment through the length of the scaffold
    # because strings (how we're dealing with sequence data here) are index_0, let's turn the repeat_occupancy_array back to index_0
    
    repeat_occupancy_array.pop(0)

    start = 0
    end = 0
    
    if window_size > len(seq):
        end = len(seq)
    else:
        end = window_size       
    
    temp_output_filename = ""
    
    if file_index == "":
        temp_output_filename = "tmp." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".txt"
    else:
        temp_output_filename = index_formatter(file_index) + "." + seq.id + ".repeats.window" + str(window_size) + ".increment" + str(increment) + ".txt"
    
    with open(temp_output_filename, "w") as temp_output_fh:
        temp_output_fh.write(xaxis_label + "\tall_repeat_content\n")
        
        while start < len(seq):
            # First, write down the chromosomal coordinate
            temp_output_fh.write(str((start + end)//2))
            
            number_of_repeated_nucleotides = sum(repeat_occupancy_array[start:end])
            
            # count the number of nucleotides that are "N"
            number_of_N = seq.seq.upper().count("N", start, end)
                
            # write these two values to the temp output file
            if end - start - number_of_N == 0:
                temp_output_fh.write("\t")
            else:
                percentage_repeated = number_of_repeated_nucleotides/(end - start - number_of_N)
                temp_output_fh.write("\t" + str(percentage_repeated))
            
            temp_output_fh.write("\n")
            
            start = start + increment
        
            if end + increment > len(seq):
                end = len(seq)
            else:
                end = end + increment  



def main(argv):
    arg_fasta = ""
    arg_gff = ""
    arg_lmer = "14"
    arg_threshold = "5"
    arg_window = 100000
    arg_increment = 100000
    
    arg_help = "\n\nRepeatScout_repeat_occupancy_table.py\n" + \
                "Uses the RepeatScout pipeline to identify de novo repeats in a genomic fasta.\n" + \
                "Ensure RepeatScout and the associated 'filter-stage-1.prl' are in the PATH.\n\nUsage: " + \
                "{0} [options] -f <genome fasta> -p <genome prefix> -l <l-mer l> -t <repeat occurrence threshold> -w <window size> -i <increment>".format(argv[0]) + \
                "\n(window size and increment default to 0, RepeatScout l-mer defaults to 14-mer, repeat occurrence threshold defaults to 5.)\n\n" + \
                "If -f argument is a list of FASTA file names separated by commas, their repeat libraries will be generated separately,\n" + \
                "but the combined library will be used to scan (RepeatMask) the concatenated sequence.\n\n"
    
    opts, args = getopt.getopt(argv[1:], "hf:p:l:t:w:i:", ["fasta=", "prefix=", "lmer=", "threshold=", "window=", "increment="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta = arg
        elif opt in ("-p", "--prefix"):
            arg_prefix = arg
        elif opt in ("-l", "--lmer"):
            arg_lmer = arg
        elif opt in ("-t", "--threshold"):
            arg_threshold = arg
        elif opt in ("-w", "--windowSize"):
            arg_window = int(arg)
        elif opt in ("-i", "--increment"):
            arg_increment = int(arg)
    line_to_write = "COMMAND: " + " ".join(argv)
    print(line_to_write)
    
    list_of_fasta = arg_fasta.split(",")
    
    print("### RepeatScout Repeat Annotation: Repeat Occupancy Pipeline ###")
    
    print("### Step 1: generating RepeatScout l-mer table...")
    if len(list_of_fasta) == 1:
        lmer_table_filename = arg_prefix + ".01.l-mer_table.out"
        if os.path.exists(lmer_table_filename):
            print("###    INFO:    l-mer table already generated. Skipping.")
        else:
            line_to_write = "build_lmer_table -l " + arg_lmer + " -sequence " + arg_fasta + " -freq " + lmer_table_filename + " > " + filename_extension_remover(lmer_table_filename) + ".log 2>&1"
            print("###    COMMAND: " + line_to_write)
            os.system(line_to_write)
    else:
        for individual_fasta in list_of_fasta:
            lmer_table_filename = arg_prefix + ".01." + filename_extension_remover(individual_fasta) + ".l-mer_table.out"
            if os.path.exists(lmer_table_filename):
                print("###    INFO:    l-mer table already generated. Skipping.")
            else:
                line_to_write = "build_lmer_table -l " + arg_lmer + " -sequence " + individual_fasta + " -freq " + lmer_table_filename + " > " + filename_extension_remover(lmer_table_filename) + ".log 2>&1"
                print("###    COMMAND: " + line_to_write)
                os.system(line_to_write)
    
    print("### Step 2: generating RepeatScout output...")
    if len(list_of_fasta) == 1:
        RepeatScout_output_filename = arg_prefix + ".02.RepeatScout.out"
        if os.path.exists(RepeatScout_output_filename):
            print("###    INFO:    RepeatScout output already generated. Skipping.")
        else:
            line_to_write = "RepeatScout -sequence " + arg_fasta + " -output " + RepeatScout_output_filename + " -freq " + lmer_table_filename + " -l " + arg_lmer + " > " + filename_extension_remover(RepeatScout_output_filename) + ".log 2>&1"
            print("###    COMMAND: " + line_to_write)
            os.system(line_to_write)
    else:
        for individual_fasta in list_of_fasta:
            lmer_table_filename = arg_prefix + ".01." + filename_extension_remover(individual_fasta) + ".l-mer_table.out"
            RepeatScout_output_filename = arg_prefix + ".02." +  filename_extension_remover(individual_fasta) + ".RepeatScout.out"
            if os.path.exists(RepeatScout_output_filename):
                print("###    INFO:    RepeatScout output already generated. Skipping.")
            else:
                line_to_write = "RepeatScout -sequence " + individual_fasta + " -output " + RepeatScout_output_filename + " -freq " + lmer_table_filename + " -l " + arg_lmer + " > " + filename_extension_remover(RepeatScout_output_filename) + ".log 2>&1"
                print("###    COMMAND: " + line_to_write)
                os.system(line_to_write)
        
    print("### Step 3: filtering RepeatScout output with filter-stage-1.prl ...")
    if len(list_of_fasta) == 1:
        RS_filter_stage_1_filename = arg_prefix + ".03.RepeatScout.filter-stage-1.out"
        if os.path.exists(RS_filter_stage_1_filename):
            print("###    INFO:    filter-stage-1.prl output already generated. Skipping.")
        else:
            line_to_write = "filter-stage-1.prl " + RepeatScout_output_filename + " > " + RS_filter_stage_1_filename + " 2> " + filename_extension_remover(RS_filter_stage_1_filename) + ".err"
            print("###    COMMAND: " + line_to_write)
            os.system(line_to_write)
    else:
        for individual_fasta in list_of_fasta:
            RepeatScout_output_filename = arg_prefix + ".02." +  filename_extension_remover(individual_fasta) + ".RepeatScout.out"
            RS_filter_stage_1_filename = arg_prefix + ".03." + filename_extension_remover(individual_fasta) + ".RepeatScout.filter-stage-1.out"
            if os.path.exists(RS_filter_stage_1_filename):
                print("###    INFO:    filter-stage-1.prl output already generated. Skipping.")
            else:
                line_to_write = "filter-stage-1.prl " + RepeatScout_output_filename + " > " + RS_filter_stage_1_filename + " 2> " + filename_extension_remover(RS_filter_stage_1_filename) + ".err"
                print("###    COMMAND: " + line_to_write)
                os.system(line_to_write)
        
    print("### Step 4: RepeatMasker masking using the RepeatScout filter-stage-1 library...")
    if len(list_of_fasta) > 1:
         # concatenate the genomic DNA fasta, and reassign arg_fasta to the concatenated genomic DNA fasta file
        print("###    INFO:    Concatenating the genome.")
        line_to_write = "cat "
        for fasta in list_of_fasta:
            line_to_write = line_to_write + fasta + " "
        arg_fasta = arg_prefix + ".genomic.fasta"
        line_to_write = line_to_write + "> " + arg_fasta
        os.system(line_to_write)
        line_to_write = "###    INFO:    Concatenated genome saved to '" + arg_fasta + "'."
        print(line_to_write)
        # if the combined library doesn't exist, figure out what's wrong
        if not os.path.exists(arg_prefix + ".03.combined_RepeatScout_libraries.out"):
            print("###    INFO:    Concatenating libraries. Please wait.")
            list_of_libraries = []
            # check to see if the filter-stage-1 libraries from individual fastas are present
            # Replace all the "R= " with unique identifier so that RepeatMasker won't get confused
            for individual_fasta in list_of_fasta:
                RS_filter_stage_1_filename = arg_prefix + ".03." + filename_extension_remover(individual_fasta) + ".RepeatScout.filter-stage-1.out"
                RS_filter_stage_1_filename_Rname_change = filename_extension_remover(RS_filter_stage_1_filename) + ".nameChange.out"
                list_of_libraries.append(RS_filter_stage_1_filename_Rname_change)
                line_to_write = "sed 's/>R=/>" + filename_extension_remover(individual_fasta) + "_R=/g' " + RS_filter_stage_1_filename + " > " + RS_filter_stage_1_filename_Rname_change
                os.system(line_to_write)
            # concatenate the libraries, and reassign RS_filter_stage_1_filename to the new concatenated library
            line_to_write = "cat "
            for library in list_of_libraries:
                line_to_write = line_to_write + library + " "
            RS_filter_stage_1_filename = arg_prefix + ".03.combined_RepeatScout_libraries.out"
            line_to_write = line_to_write + "> " + RS_filter_stage_1_filename
            os.system(line_to_write)
           
        else:
            line_to_write = "###    INFO:    Concatenated library, '" + arg_prefix + ".03.combined_RepeatScout_libraries.out', exists. RepeatMasking."
            RS_filter_stage_1_filename = arg_prefix + ".03.combined_RepeatScout_libraries.out"
            print(line_to_write)
    else:
        pass
    RM_dir = arg_prefix + "_RepeatMasker"
    RM_outfile = RM_dir + "/" + arg_fasta + ".out"
    RM_gff = arg_prefix + ".04.RepeatMasker.gff"
    if os.path.exists(RM_outfile) and os.path.exists(RM_gff):
        print("###    INFO:    RepeatMasker appears to have been completed and the gff generated. Skipping.")
    elif os.path.exists(RM_outfile):
        print("###    INFO:    RepeatMasker appears to have been completed previously. Generating GFF file...")
        line_to_write = "bash rmOutToGFF3custom_20231021.bash -o " + RM_outfile + " > " + RM_gff
        print("###    COMMAND: " + line_to_write)
        os.system(line_to_write)
        line_to_write = "sed -i 's/R=/R_/g' " + RM_gff
        os.system(line_to_write)
    else:
        if len(list_of_fasta) > 1:
            print("###    INFO:    Running RepeatMasker with the filtered and concatenated library from Step 3...")
        else:
            print("###    INFO:    Running RepeatMasker with the filtered library from Step 3...")
        RM_log = arg_prefix + ".04.RepeatMasker.log"
        line_to_write = "RepeatMasker -e ncbi -pa 12 -lib " + RS_filter_stage_1_filename + " -dir " + RM_dir + " " + arg_fasta + " > " + RM_log + " 2>&1"
        print("###    COMMAND: " + line_to_write)
        os.system(line_to_write)
        print("###    INFO:    RepeatMasker completed. Generating GFF file...")
        line_to_write = "bash rmOutToGFF3custom_20231021.bash -o " + RM_outfile + " > " + RM_gff
        print("###    COMMAND: " + line_to_write)
        os.system(line_to_write)
        line_to_write = "sed -i 's/R=/R_/g' " + RM_gff
        os.system(line_to_write)
        
    
    print("### Step 5: determining the repeat occurrence and total repeat occupancy for each repeat type...")
    RM_features_list = arg_prefix + ".05.RepeatMasker.featuresList.txt"
    if os.path.exists(RM_features_list):
        print("###    INFO:    File detailing all the features has already been generated. Skipping.")
    else:
        combined_gff_feature_tally(RM_gff, RM_features_list)
    
    line_to_write = "### Step 6: limiting annotation to repeats that occur more than the threshold = " + arg_threshold + "..."
    print(line_to_write)
    new_gff = arg_prefix + ".06.RepeatMasker.repeatsOccurringMoreThan" + arg_threshold + "Times.gff"
    if os.path.exists(new_gff):
        print("###    INFO:    The new GFF with repeats has already been generated. Skipping.")
    else:
        discard_rare_repeats(RM_gff, RM_features_list, arg_threshold, new_gff)
    
    
    print("### Step 7: Generating a table classifying the repeats...")
    repeat_classification_file = arg_prefix + ".07.repeat_classfication.txt"
    if os.path.exists(repeat_classification_file):
        print("###    INFO:    Repeat classification file exists. Skipping.")
    else:
        auto_repeat_classification(RM_features_list, repeat_classification_file)
    
    print("### Step 8: Generating repeat occupancy tables (one file per chromosome/scaffold/contig)...")
    file_index = arg_prefix + ".08"
    for seq in SeqIO.parse(arg_fasta, "fasta"):
        count_RepeatMasker_masked_content_v20231026(seq, new_gff, repeat_classification_file, arg_window, arg_increment, "coordinate", file_index)
    
    '''print("### Step 9: Organizing files...")
    line_to_write = "mkdir " + arg_prefix + "_analysis"
    os.system(line_to_write)
    line_to_write = "mv " + arg_prefix + ".* " + arg_prefix + "_analysis"
    os.system(line_to_write)
    # move the RepeatMasker folder as well
    line_to_write = "mv " + arg_prefix + "_RepeatMasker " + arg_prefix + "_analysis"
    os.system(line_to_write)
    '''
    print("### Pipeline finished.\n\n")
    
if __name__ == "__main__":
    main(sys.argv)