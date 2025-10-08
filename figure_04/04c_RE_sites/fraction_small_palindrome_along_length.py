################################################################
# fraction_small_palindrome_along_length.py
# Calculates the occupancy of a small palindromic sequence 
# along the length of contigs in a fasta file.
# Allows the user to adjust the window size and increment
# 
################################################################


from Bio import SeqIO
import sys, getopt


# intended to count fraction_of_NlaIII_sites coverage for PoreC
# calculates the COVERAGE of NlaIII sites (1 site = 4 bps)
def count_fraction_of_small_seq_sites_coverage(seq, small_seq_string, window_size, increment, xaxis_label, yaxis_label, file_prefix):
    
    length_of_seq = len(seq)
    
    start = 0
    end = 0
    
    if window_size > length_of_seq:
        end = length_of_seq
    else:
        end = window_size
  
    temp_output_filename = file_prefix + "." + seq.id + "." + yaxis_label + ".window" + str(window_size) + ".increment" + str(increment) + ".txt"
    
    small_seq_string = small_seq_string.upper()
        
    with open(temp_output_filename, "w") as temp_output_fh:
    
        string_to_write = xaxis_label + "\t" + yaxis_label + "\n"
        temp_output_fh.write(string_to_write)
        
        # go through all the windows
        while start < length_of_seq:
            # NlaIII RE site coverage is the number of NlaIII sites x 4 (bps per site)
            small_seq_site_cov = seq.seq.upper().count(small_seq_string, start, end)*len(small_seq_string)
            
            # account for Ns (scaffolding gaps)
            number_of_N = seq.seq.upper().count("N", start, end)
            
            temp_output_fh.write(str((start + end)//2) + "\t")
            
            # account for if the window only contains Ns
            if end-start-number_of_N == 0:
                temp_output_fh.write("\n")
            else:
                # fraction coverage by NlaIII site
                fraction_small_seq_coverage = float(small_seq_site_cov/(end-start-number_of_N))
                #fraction_small_seq_coverage = round(fraction_small_seq_coverage, 8)
                temp_output_fh.write(str(fraction_small_seq_coverage)+"\n")
            
            start = start + increment
            
            if end + increment > length_of_seq:
                end = length_of_seq
            else:
                end = end + increment


def main(argv):
    arg_fasta = ""
    arg_seq = ""
    arg_prefix = "tmp"
    arg_x = "contig_position"
    arg_y = ""
    arg_window = 100000
    arg_increment = 100000
    
    arg_help = "\n\n{0}\n".format(argv[0]) + \
                "Calculates the fraction occupancy of a small palindromic sequence along the length of sequences in a fasta file.\n" + \
                "Every contig produces one output file.\n\nUsage: " + \
                "{0} [options] -f <genome fasta> -s <short palindromic sequence> -p <output file prefix> -x <x-axis label> -y <y-axis label> -w <window size> -i <increment>".format(argv[0]) + \
                "\n\nrequired arguments:" + \
                "\n-f, or the fasta file containing sequences." + \
                "\n-s, or the short palindromic sequence to look for." + \
                "\nwindow size -w and increment -i default to 100000.\n\n"
    
    opts, args = getopt.getopt(argv[1:], "hf:s:p:x:y:w:i:", ["fasta=", "shortPalindrome=", "prefix=", "xaxis=", "yaxis=", "windowSize=", "increment="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta = arg
        elif opt in ("-s", "--shortPalindrome"):
            arg_seq = arg
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
        count_fraction_of_small_seq_sites_coverage(seq, arg_seq, arg_window, arg_increment, arg_x, arg_y, arg_prefix)


if __name__ == "__main__":
	main(sys.argv)