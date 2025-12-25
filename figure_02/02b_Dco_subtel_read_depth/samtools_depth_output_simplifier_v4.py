# Python script "samtools_depth_output_simplifier_v4.py"
# Simplifies the samtools depth output to a window size of your choosing
# can specify the scaffold and the range


from Bio import SeqIO
import sys, os, getopt

def file_extension_remover(filename):
    list_of_filename_parts = filename.split(".")
    new_filename = ""
    
    for index, filename_part in enumerate(list_of_filename_parts):
        if index < len(list_of_filename_parts) - 1:
            new_filename = new_filename + filename_part + "."
        else:
            pass
    
    return new_filename    


def main(argv):
    arg_depth = ""
    arg_chrom = ""
    arg_window = 10000
    arg_help = "{0} [options] -d <samtools depth output> -w <window size> -c <chromosome and range to show> \n".format(argv[0]) + \
                "-c option: a chromosome ranges file, if provided, will generate a second set of files that has the depth information contained in those ranges. " + \
                "The chromosome ranges file format is <chromosome name> <start coordinate> <end coordinate> separated by tabs, one range per line.\n" + \
                "Window size defaults to 10000 unless specified.\n"
    
    # get the options
    opts, args = getopt.getopt(argv[1:], "hd:w:c:", ["help", "depth=", "window=", "chrom="])
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-d", "--depth"):
            arg_depth = arg
        elif opt in ("-w", "--window"):
            arg_window = int(arg)
        elif opt in ("-c", "--chrom"):
            arg_chrom = arg
    
    #############################################################################################################################################################
    # generate the simplified samtools depth file                                                                                                               #
    #############################################################################################################################################################
    
    simplified_samtools_depth_output = file_extension_remover(arg_depth) + "window" + str(arg_window) + ".txt"
    
    if os.path.exists(simplified_samtools_depth_output) == True:
        pass
    else:
        with open(arg_depth, "r") as samtools_depth_output_filehandle, open(simplified_samtools_depth_output, "w") as samtools_depth_simplified_filehandle:
            contig_name = ""
            current_position = 0
            sum_coverage = 0
            plot_x = 0
            plot_y = 0
            
            samtools_depth_simplified_filehandle.write("#scaffold\tcoordinate\tread_depth\n")
            
            line = samtools_depth_output_filehandle.readline()
            
            while line:
                line_items = line.split("\t")
                if contig_name != line_items[0]:
                    if current_position % arg_window != 0:
                        plot_x = current_position - (current_position % arg_window)/2
                        plot_y = sum_coverage/(current_position % arg_window)
                        
                        string_to_write = contig_name + "\t" + str(plot_x) + "\t" + str(plot_y) + "\n"
                        samtools_depth_simplified_filehandle.write(string_to_write)
                    else:
                        pass
                    contig_name = line_items[0]
                    current_position = int(line_items[1])
                    sum_coverage = int(line_items[2].rstrip())
                
                else:
                    if current_position % arg_window == 0:
                        plot_x = current_position - arg_window/2
                        plot_y = sum_coverage/arg_window
                        
                        string_to_write = contig_name + "\t" + str(plot_x) + "\t" + str(plot_y) + "\n"
                        samtools_depth_simplified_filehandle.write(string_to_write)
                        sum_coverage = 0
                    else:
                        pass
                    current_position = int(line_items[1])
                    sum_coverage = sum_coverage + int(line_items[2].rstrip())
                
                line = samtools_depth_output_filehandle.readline()
            
            plot_x = current_position - (current_position % arg_window)/2
            print(current_position)
            print(arg_window)
            plot_y = sum_coverage/(current_position % arg_window)
            string_to_write = contig_name + "\t" + str(plot_x) + "\t" + str(plot_y) + "\n"
            samtools_depth_simplified_filehandle.write(string_to_write)
    
    
    
    #############################################################################################################################################################
    # now generate simplified depth files for the ranges specified in the -c option
    #############################################################################################################################################################
    
    if arg_chrom != "":
        with open(arg_chrom, "r") as chromosome_ranges_fh:
            line = chromosome_ranges_fh.readline()
            
            while line:
                # generate a new file per line, describing all the coordinates with values that are within the chromosome ranges in that line
                contig_name = line.split("\t")[0]
                start_coordinate = int(line.split("\t")[1])
                end_coordinate = int(line.split("\t")[2])
                
                new_depth_file = file_extension_remover(simplified_samtools_depth_output) + contig_name + "." + str(start_coordinate) + "-" + str(end_coordinate) + ".txt"
                
                with open(simplified_samtools_depth_output, "r") as simplified_samtools_depth_fh, open(new_depth_file, "w") as new_depth_fh:
                    simplified_samtools_depth_output_line = simplified_samtools_depth_fh.readline()
                    
                    while simplified_samtools_depth_output_line:
                        if simplified_samtools_depth_output_line[0] == "#":
                            new_depth_fh.write(simplified_samtools_depth_output_line)
                        elif simplified_samtools_depth_output_line.split("\t")[0] == contig_name and float(simplified_samtools_depth_output_line.split("\t")[1]) >= start_coordinate and float(simplified_samtools_depth_output_line.split("\t")[1]) <= end_coordinate:
                            new_depth_fh.write(simplified_samtools_depth_output_line)
                        else:
                            pass
                        simplified_samtools_depth_output_line = simplified_samtools_depth_fh.readline()
                        
                line = chromosome_ranges_fh.readline()
    else:
        pass

    print("\nscript finished!\n\n")



if __name__ == "__main__":
    main(sys.argv)	
	
		
	
	