################################################################
# fraction_GC_along_length.py
# Calculates the %GC content along the length of sequences in
# a fasta file.
# Allows the user to adjust the window size and increment
# 
################################################################


from Bio import SeqIO
import sys, getopt

# count_percent_GC(seq, window_size, increment, xaxis_label, yaxis_label, file_prefix)
# Generates one file per contig. seq is a SeqRecord object.
# Calculates the GC content along the length of the seq. Outputs a table with two columns with headers:
# xaxis_label (usually "coordinate" or "scaffold_position" etc.) and yaxis_label (usually "fraction_GC")
def count_percent_GC(seq, window_size, increment, xaxis_label, yaxis_label, file_prefix):
    percent_GC_list = []
    
    length_of_seq = len(seq)
    
    start = 0
    end = 0
    
    if window_size > length_of_seq:
        end = length_of_seq
    else:
        end = window_size
        
    temp_output_filename = file_prefix + "." + seq.id + "." + yaxis_label + ".window" + str(window_size) + ".increment" + str(increment) + ".txt"
    
    # write the percent GC to temp_output_filename
    
    with open(temp_output_filename, "w") as temp_output_fh:
        string_to_write = xaxis_label + "\t" + yaxis_label + "\n"
        temp_output_fh.write(string_to_write)
        
        while start < length_of_seq:
            number_of_G = seq.seq.upper().count("G", start, end)
            number_of_C = seq.seq.upper().count("C", start, end)
        
            number_of_N = seq.seq.upper().count("N", start, end)
            
            temp_output_fh.write(str((start + end)//2) + "\t")
            
            # if the current window only contains Ns
            if end - start - number_of_N == 0:
                temp_output_fh.write("\n")
            
            # if the current window actually has some sequence data
            else:
                percentage_GC = (number_of_G + number_of_C)/(end - start - number_of_N)
                # write to file this value
                temp_output_fh.write(str(percentage_GC) + "\n")           
        
            start = start + increment
        
            if end + increment > length_of_seq:
                end = length_of_seq
            else:
                end = end + increment

def main(argv):
    arg_fasta = ""
    arg_prefix = "tmp"
    arg_x = "contig_position"
    arg_y = "fraction_GC"
    arg_window = 100000
    arg_increment = 100000
    
    arg_help = "\n\n{0}\n".format(argv[0]) + \
                "Calculates the fraction GC along the length of sequences in a fasta file.\n" + \
                "Every contig produces one output file.\n\nUsage: " + \
                "{0} [options] -f <genome fasta> -p <output file prefix> -x <x-axis label> -y <y-axis label> -w <window size> -i <increment>".format(argv[0]) + \
                "\nrequired argument -f, or the fasta file containing sequences." + \
                "\nwindow size -w and increment -i default to 100000.\n\n"
    
    opts, args = getopt.getopt(argv[1:], "hf:p:x:y:w:i:", ["fasta=", "prefix=", "xaxis=", "yaxis=", "windowSize=", "increment="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta = arg
        elif opt in ("-p", "--prefix"):
            arg_prefix = arg
        elif opt in ("-x", "--xaxis"):
            arg_x = arg
        elif opt in ("-y", "--yaxis"):
            arg_y = arg
        elif opt in ("-w", "--windowSize"):
            arg_window = int(arg)
        elif opt in ("-i", "--increment"):
            arg_increment = int(arg)
    line_to_write = "COMMAND: " + " ".join(argv)
    print(line_to_write)
    
    for seq in SeqIO.parse(arg_fasta, "fasta"):
        count_percent_GC(seq, arg_window, arg_increment, arg_x, arg_y, arg_prefix)


if __name__ == "__main__":
	main(sys.argv)