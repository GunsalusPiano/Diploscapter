####################################################################################################
# PCA_bedgraph_converter.py
####################################################################################################
#
# transform the PCA bedgraph file to another file compatible with the rest of the data tables
# bedgraph file format
# CHROM (ALL CAPS)      start (exclusive)       end (inclusive)      value
# The bedgraph file should already be sorted by CHROM and start
# make sure the bedgraph file has the same resolution as window_size

import sys, getopt
from Bio import SeqIO
from statistics import mean

def PCA_bedgraph_converter(current_seq, PCA_bedgraph_filename, window_size, increment, xaxis_label, yaxis_label, output_file_prefix):
    
    # construct PCA_value_array
    # stores all PCA values across every position in current_seq
    # this is to facilitate the conversion to ANY windows_size
    # make it an 1-index array
    
    length_of_seq = len(current_seq)
    
    PCA_value_array = [0]*(length_of_seq+1)
    
    # populate the array with values from PCA_bedgraph_file
    with open(PCA_bedgraph_filename, "r") as PCA_bedgraph_fh:
        line = PCA_bedgraph_fh.readline()
        
        while line:
            line_items = line.split()
            scaffold = line_items[0].lower()
                
            # only consider lines annotating the current seq (scaffold)
            # ignore all the rest
            if scaffold == current_seq.id.lower():
                start = int(line_items[1])
                end = int(line_items[2])
                value = float(line_items[3])
                
                # we want to turn the 0s in positions start to end to turn into 1s
                # start to end INCLUSIVE (python range() is exclusive on the right end)
                # bedgraph (I think) is exclusive on the left end
                for x in range(start+1, end+1):
                    PCA_value_array[x] = value
            else:
                pass
            
            line = PCA_bedgraph_fh.readline()
            
    
    # turn PCA_value_array back to 0-index
    PCA_value_array.pop(0)
    
    start = 0
    end = 0
    
    if window_size > len(current_seq):
        end = len(current_seq)
    else:
        end = window_size
    
    temp_output_filename = output_file_prefix + "." + current_seq.id + "." + yaxis_label + ".window" + str(window_size) + ".increment" + str(increment) + ".txt"
    
    with open(temp_output_filename, "w") as temp_output_fh:
        line_to_write = xaxis_label + "\t" + yaxis_label + "\n"
        temp_output_fh.write(line_to_write)
        
        while start < length_of_seq:
            
            # take the average from start to end
            
            value = mean(PCA_value_array[start:end])
            line_to_write = str((start + end)//2) + "\t" + str(value) + "\n"
            temp_output_fh.write(line_to_write)
                    
            start = start + increment
                    
            if end + increment > len(current_seq):
                end = len(current_seq)
            else:
                end = end + increment



def main(argv):
    arg_fasta = ""
    arg_bedgraph = ""
    arg_window = 100000
    arg_increment = 100000
    arg_xaxis = ""
    arg_yaxis = ""
    arg_prefix = ""
    
    arg_help = "\n\n{0}\n".format(argv[0]) + \
                "Turns a PCA bedgraph file (from hicExplorer) into a data table with fixed window and increment.\n" + \
                "\n\nUsage: " + \
                "{0} [options] -f <genome fasta> -b <PCA bedgraph> -w <window size> -i <increment> -x <x axis label> -y <y axis label> -p <output file prefix>".format(argv[0]) + \
                "\n\nWindow size and increment default to 100000.\n" + \
                "all other arguments required. \n\n"
    
    opts, args = getopt.getopt(argv[1:], "hf:b:w:i:x:y:p:", ["fasta=", "bedgraph=", "window=", "increment=", "xaxis=", "yaxis=", "prefix="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta = arg
        elif opt in ("-b", "--bedgraph"):
            arg_bedgraph = arg
        elif opt in ("-w", "--window"):
            arg_window = int(arg)
        elif opt in ("-i", "--increment"):
            arg_increment = int(arg)
        elif opt in ("-x", "--xaxis"):
            arg_xaxis = arg
        elif opt in ("-y", "--yaxis"):
            arg_yaxis = arg            
        elif opt in ("-p", "--prefix"):
            arg_prefix = arg
    line_to_write = "COMMAND: " + " ".join(argv)
    print(line_to_write)
	
    for current_seq in SeqIO.parse(arg_fasta, "fasta"):
        PCA_bedgraph_converter(current_seq, arg_bedgraph, arg_window, arg_increment, arg_xaxis, arg_yaxis, arg_prefix)

if __name__ == "__main__":
    main(sys.argv)