#BUSCO_protein_full_table_to_BUSCO_genomic_full_table.py
#this script takes a BUSCO protine full_table.tsv
#and turns it into a BUSCO genomic full_table equivalent
#by importing information from a provided GTF (genomic annotation) file

import sys, getopt

def extract_StringTie_attribute(line, attribute_type):
    line_items = line.split("\t")
    # line_items[8] are the attributes
    # first eliminate all the spaces
    attributes_list = line_items[8].replace(" ", "").rstrip()
    #print(attributes_list)
    #input(attributes_list + " top of extract function")
    
    # initiate an attributes_dict dictionary
    attributes_dict = {}
    
    # go through each attribute, which is separated by ";"
    for item in attributes_list.split(";"):
        #input(item + " in the for loop")
        if len(item.split('"')) > 1:
            #input(item.split('"')[0] + "\t" + item.split('"')[1])
            attributes_dict[item.split('"')[0]] = item.split('"')[1]
        else:
            pass
    
    if attribute_type in attributes_dict:
        return attributes_dict[attribute_type]
    else:
        return "null"
        
    

def turn_protein_BUSCO_to_genomic_BUSCO(protein_busco_table, gtf, output_file_prefix):
    #First, generate a dictionary from data in the GTF
    #where key is the transcript name
    #the value is a list containing [scaffold_name, start, end, strand] values
    
    genomic_busco_filename = output_file_prefix + ".genomic_busco_equiv_of_" + protein_busco_table
    gtf_dictionary = {}
    
    with open(gtf, 'r') as gtf_filehandle:
        line = gtf_filehandle.readline()
        
        while line:
            line_items = line.split()
            
            #split the line: pertinent information below
            #line_items[0] = scaffold name
            #line_items[2] = type of feature
            #line_items[3] = start position
            #line_items[4] = end position
            #line_items[6] = strandedness (+ or -)
            #line_items[8] = transcript attributes (gxxxxxxx.t1)
            
            if line_items[2] == "transcript":
                scaffold_name = line_items[0]
                start_position = line_items[3]
                end_position = line_items[4]
                strandedness = line_items[6]
                transcript_name = extract_StringTie_attribute(line, "transcript_id")
                #print(scaffold_name, start_position, end_position, strandedness, transcript_name)
                
                gtf_dictionary[transcript_name] = [scaffold_name, start_position, end_position, strandedness]
            else:
                pass
            
            line = gtf_filehandle.readline()
    
    #now, read through the full_table.tsv generated from the protein BUSCO
    #and write lines to a new 'genomic BUSCO' table equivalent
    with open(protein_busco_table, 'r') as protein_busco_filehandle, open(genomic_busco_filename, 'w') as genomic_busco_filehandle:
        line = protein_busco_filehandle.readline()
        line_to_write = ""
        while line:
            if line[0] == "#":
                genomic_busco_filehandle.write(line)
            else:
                line_items = line.split("\t")
                #line_items[0] is the busco id
                #line_items[1] is the match status (complete, duplicated, fragmented, or missing)
                #line_items[2] is the transcript name
                #line_items[3] is the match score
                #line_items[4] is the length of match
                #line_items[5] and [6] may not be present
                
                if line_items[1].rstrip() == "Missing":
                    line_to_write = line
                else:
                    line_to_write = line_items[0] + "\t" + line_items[1] + "\t" + gtf_dictionary[line_items[2]][0] + "\t" + gtf_dictionary[line_items[2]][1] + "\t" + gtf_dictionary[line_items[2]][2] + "\t" + gtf_dictionary[line_items[2]][3] + "\t" + line_items[3] + "\t" + line_items[4]
                
                    if len(line_items) > 5:
                        line_to_write = line_to_write + "\t" + line_items[5] + "\t" + line_items[6]
                    else:
                        pass
            
            genomic_busco_filehandle.write(line_to_write)
            line = protein_busco_filehandle.readline()
            

def main(argv):
    arg_busco = ""
    arg_gff = ""
    arg_prefix = ""
    
    arg_help = "\n\n{0}\n".format(argv[0]) + \
                "Generates a genomic version of a BUSCO table using a protein BUSCO table and a GFF.\n" + \
                "\n\nUsage: " + \
                "{0} [options] -b <protein busco table tsv> -g <gff containing gene annotations> -p <output file prefix>".format(argv[0]) + \
                "\n\nall arguments are required.\n\n"
    
    opts, args = getopt.getopt(argv[1:], "hb:g:p:", ["busco=", "gff=", "prefix=", "xaxis="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-b", "--busco"):
            arg_busco = arg
        elif opt in ("-g", "--gff"):
            arg_gff = arg
        elif opt in ("-p", "--prefix"):
            arg_prefix = arg
    line_to_write = "COMMAND: " + " ".join(argv)
    print(line_to_write)
	
	#protein_busco_table = "full_table.tsv"
	#protein_busco_table = "Dpa1_full_table.tsv"
	#protein_busco_table = "Dpa2_full_table.tsv"
	#protein_busco_table = "Dco1_full_table.tsv"
	
	
	#gtf = "augustus.hints.gtf"
	#gtf = "temp.v202304_scaffold_1.gtf.t1_transcripts_only.gtf"
	#gtf = "temp.v202304_scaffold_2.gtf.t1_transcripts_only.gtf"
	#gtf = "Dco.augustus.hints.gtf.t1_transcripts_only.gtf"
	
    turn_protein_BUSCO_to_genomic_BUSCO(arg_busco, arg_gff, arg_prefix)
	


if __name__ == "__main__":
	main(sys.argv)