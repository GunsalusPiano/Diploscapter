#parse through the output of samtools depth
#The samtools depth output contains three columns
#Scaffold/contig name, position, and depth of bam/sam coverage
#This script extracts the coverage information (3rd column) for all scaffolds/contigs


print("\nParser for samtools depth output.\n\n")


#Enter filename for samtools depth output file here.

samtools_depth_output = input("samtools depth output filename:\n")

#Enter filename for output file for this script here

output_filename = samtools_depth_output + ".tabulated.txt"

cumulative_depth = 0
#array to store depth values
depth_values_array = []

with open(samtools_depth_output, "r") as filehandle:
	#print("in with block")
	while True:
		#print("in 1st while block")
		line = filehandle.readline()
		#print("this is the current line: ", line, "\n")
		if not line:
			break		
		line_items = line.split("\t")
		current_depth_value = int(line_items[2])
		while len(depth_values_array) < current_depth_value+1:
			depth_values_array.append(0)
		depth_values_array[current_depth_value]=depth_values_array[current_depth_value]+1
		cumulative_depth=cumulative_depth+current_depth_value



with open(output_filename, "w") as filehandle:
	filehandle.write("depth\tnum_of_positions\n")
	for index, depth in enumerate(depth_values_array):
		string_to_write = str(index) + "\t" + str(depth) + "\n"
		filehandle.write(string_to_write)

	
	
	