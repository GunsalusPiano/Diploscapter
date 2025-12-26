# repeat_nucmer_and_delta_simplified_dotplot_prep_for_fig01.py
# marks every position of a match
# outputs in the format of show-coords -T

from Bio import SeqIO
import sys
import getopt
import os

# filename_extension_remover(filename)
# removes the extension (last string after .) of a file name
def filename_extension_remover(filename):
	filename_split = filename.split(".")
	
	filename_without_extension = ""
	
	for i in range(len(filename_split)-1):
		filename_without_extension = filename_without_extension + filename_split[i] + "."
	
	return filename_without_extension[0:len(filename_without_extension)-1]

def main(argv):
    arg_fasta_file = ""
    arg_help = "{0} [options] -f <list of fasta files in text file>".format(argv[0]) + \
                "\n\nGenerates nucmer delta files from every pairwise comparison of FASTAs in the -f file." + \
                "\nThe -f file should contain one fasta file name per line.\n\n"
    opts, args = getopt.getopt(argv[1:], "hf:", ["help", "fasta="])

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta_file = arg
    
    list_of_fastas = []
    
    with open(arg_fasta_file, "r") as arg_fasta_fh:
        line = arg_fasta_fh.readline()
        while line:
            list_of_fastas.append(line.rstrip())
            line = arg_fasta_fh.readline()
    
    # run nucmer for every pairwise comparison
    for reference in list_of_fastas:
        for query in list_of_fastas:
            if reference == query:
                pass
            else:
                # set delta file prefix
                delta_prefix = filename_extension_remover(reference) + "_vs_" + filename_extension_remover(query)
                
                # run nucmer
                line_to_write = "nucmer -p " + delta_prefix + " " + reference + " " + query
                os.system(line_to_write)
                
                # run delta-filter
                prefiltered_delta_file = delta_prefix + ".delta"
                filtered_delta_file = delta_prefix + ".delta-filter.options_rql0000.delta"
                
                line_to_write = "delta-filter -r -q -l 10000 " + prefiltered_delta_file + " > " + filtered_delta_file
                os.system(line_to_write)
                
                # parse through the delta to make a .coords file for graphing in R
                arg_coords_file = filename_extension_remover(filtered_delta_file) + ".coords"
                
                with open(filtered_delta_file, "r") as delta_fh, open(arg_coords_file, "w") as coords_fh:
                    line = delta_fh.readline()
                    # first line lists the files nucmer was run on
                    # assume that these files are still in the same folder
                    
                    # transfer the line to the coords file
                    #coords_fh.write(line)
                    
                    # second line should just read "NUCMER"
                    line = delta_fh.readline()
                    #coords_fh.write(line)
                    #coords_fh.write("\n")
                    coords_fh.write("rs\tre\tqs\tqe\terror\trid\tqid\tstrand\n")
                    
                    ref_start = 0
                    qry_start = 0
                    
                    ref_end = 0
                    qry_end = 0
                    strandedness = "+"
                    
                    ref_seq_id = ""
                    qry_seq_id = ""
                    
                    line = delta_fh.readline()
                    line_number = 3
                    
                    while line:
                        # identify the header lines
                        # grab the sequence strings
                        
                        if line[0] == ">":
                            line_items = line[1:].split()
                            # line_items[0] is the reference seq id
                            # line_items[1] is the query seq id
                            # line_items[2] is the length of ref
                            # line_items[3] is the length of qry
                            
                            ref_seq_id = line_items[0]
                            qry_seq_id = line_items[1]
                            
                        # if NOT a header line (the ones starting with >)
                        else:
                            line_items = line.split(" ")
                            line_to_write = ""
                            # if it's an summary line, it should have 7 items
                            if len(line_items) > 1:
                                ref_start = int(line_items[0])
                                ref_end = int(line_items[1])
                                
                                qry_start = int(line_items[2])
                                qry_end = int(line_items[3])
                                
                                errors = int(line_items[4])
                                
                                # determine the strandedness
                                if qry_start > qry_end:
                                    strandedness = "-"
                                else:
                                    strandedness = "+"
                                
                                line_to_write = str(ref_start) + "\t" + str(ref_end) + "\t" + str(qry_start) + "\t" + str(qry_end) + "\t" + str(errors) + "\t" + ref_seq_id + "\t" + qry_seq_id + "\t" + strandedness + "\n"
                                coords_fh.write(line_to_write)
                                
                            # it's an alignment line
                            else:
                                pass
                                
                        line = delta_fh.readline()

if __name__ == "__main__":
    main(sys.argv)