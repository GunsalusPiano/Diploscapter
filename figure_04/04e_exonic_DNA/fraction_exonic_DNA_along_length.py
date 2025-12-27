################################################################
# fraction_exonic_DNA_along_length.py
# Calculates the exonic DNA along the length of sequences in
# a fasta file.
# Allows the user to adjust the window size and increment
# 
################################################################

from Bio import SeqIO
import sys, getopt

# count the fraction exonic DNA using an augustus.hints.gff output
# probably works best with gff with only t1 transcripts!
# but the code probably won't break if there are multiple transcripts per gene
def count_fraction_exonic_DNA(seq, augustus_gff_output, window_size, increment, xaxis_label, yaxis_label, file_prefix):
    # exon_position_array
    # array with length of sequence number of items + 1, prepresenting the nucleotides in the seq
    # exon positions = 1
    # everywhere else = 0
    
    length_of_seq = len(seq)
    
    exon_position_array = [0]*(length_of_seq+1)
    
    with open(augustus_gff_output, "r") as augustus_gff_output_fh:
        
        # read the augustus gff file line by line
        line = augustus_gff_output_fh.readline()
        
        while line:
            # ignore lines starting with #
            if line[0] == "#":
                pass
            else:
                line_items = line.split()
                feature = line_items[2]
                
                # only deal with lines that read are "exon" annotations
                if feature == "exon":
                    scaffold = line_items[0]
                    
                    # only deal with lines that pertain to the current seq
                    if scaffold == seq.id:
                        # capture the exon start and exon end positions and store in the exon_position_array
                        # have position go from exon start (line_items[3]) to exon end (line_items[4])
                        position = int(line_items[3])
                        
                        while position <= int(line_items[4]):
                            exon_position_array[position] = 1
                            position = position + 1
                    else:
                        pass
                else:
                    pass
            line = augustus_gff_output_fh.readline()
    
    # now deal with the exonic DNA content per window_size
    # make the exon_position_array 0-based 
    
    exon_position_array.pop(0)
        
    start = 0
    end = 0
    
    if window_size > length_of_seq:
        end = length_of_seq
    else:
        end = window_size       
    
    # record exonic DNA occupancy in a file
    temp_output_filename = file_prefix + "." + seq.id + "." + yaxis_label + ".window" + str(window_size) + ".increment" + str(increment) + ".txt"
    
    with open(temp_output_filename, "w") as temp_output_fh:
        
        string_to_write = xaxis_label + "\t" + yaxis_label + "\n"
        temp_output_fh.write(string_to_write)
        
        while start < length_of_seq:
            # count the number of exonic nucleotides in the range
            number_of_exon_bps = sum(exon_position_array[start:end])
        
            # count the number of nucleotides that are "N"
            number_of_N = seq.seq.upper().count("N", start, end)
            
            # start writing to file: x-axis value
            temp_output_fh.write(str((start + end)//2) + "\t")
            
            # account for if the window contains all Ns
            if end - start - number_of_N == 0:
                temp_output_fh.write("\n")
            else:
                fraction_exonic_DNA = number_of_exon_bps/(end - start - number_of_N)
                temp_output_fh.write(str(fraction_exonic_DNA)+"\n")
            
            start = start + increment
        
            if end + increment > length_of_seq:
                end = length_of_seq
            else:
                end = end + increment


def main(argv):
    arg_fasta = ""
    arg_gff = ""
    arg_prefix = "tmp"
    arg_x = "contig_position"
    arg_y = "fraction_exonic_DNA"
    arg_window = 100000
    arg_increment = 100000
    
    arg_help = "\n\n{0}\n".format(argv[0]) + \
                "Calculates the fraction occupancy of exonic DNA along the length of sequences in a fasta file.\n" + \
                "Every contig produces one output file.\n\nUsage: " + \
                "{0} [options] -f <genome fasta> -g <gff containing gene annotations> -p <output file prefix> -x <x-axis label> -y <y-axis label> -w <window size> -i <increment>".format(argv[0]) + \
                "\n\nrequired arguments:" + \
                "\n-f, or the fasta file containing sequences." + \
                "\n-g, or the gff file containing gene annotations." + \
                "\nwindow size -w and increment -i default to 100000.\n\n"
    
    opts, args = getopt.getopt(argv[1:], "hf:g:p:x:y:w:i:", ["fasta=", "gff=", "prefix=", "xaxis=", "yaxis=", "windowSize=", "increment="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-f", "--fasta"):
            arg_fasta = arg
        elif opt in ("-g", "--gff"):
            arg_gff = arg
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
        count_fraction_exonic_DNA(seq, arg_gff, arg_window, arg_increment, arg_x, arg_y, arg_prefix)


if __name__ == "__main__":
	main(sys.argv)