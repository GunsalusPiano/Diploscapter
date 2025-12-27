# Diploscapter_telomere_gggenomes_visualization.py

# organizes the data from *subtel.sorted.coords
# and *subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop.txt
#
# into something graphable

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

class Read:
    # read_length: INT read length in bp
    # max_dropoff: INT the position on the read where occupancy of telomeric repeats is the most different to the right vs to the left
    # orientation: INT + for pointing to the right, - for pointing to the left
    def __init__(self, read_name, read_length, max_dropoff, orientation):
        self.read_name = read_name
        self.read_length = read_length
        self.max_dropoff = max_dropoff
        self.orientation = orientation
        self.read_end = 0
        
        self.subtel_list = [] # each element will be [subtel_name, start, end]
        
        
        if orientation > 0:
            self.read_end = read_length
        else:
            self.read_end = 1
        
    def reorient(self):
        if self.orientation > 0:
            self.subtel_list.sort(key=lambda x: x[1])
        else:
            new_list = []
            for feature in self.subtel_list:
                feature[1] = self.read_length - feature[1] # start position
                feature[2] = self.read_length - feature[2] # end position
            
            self.subtel_list.sort(key=lambda x: x[1])
            
            self.max_dropoff = self.read_length - self.max_dropoff
            
            self.orientation = -self.orientation
            self.read_end = self.read_length
    
    # subtel_name is STR the name of the subtelomeric repeat
    # start and end are INT positions of start and end of a subtelomeric repeat
    def add_start_and_end(self, subtel_name, start, end):
        self.subtel_list.append([subtel_name, start, end])
        self.subtel_list.sort(key=lambda x: x[1])
        
    # nudge features over by bp
    def move_features(self, bp):
        new_list = []
        
        for item in self.subtel_list:
            item[1] = item[1] + bp # start position shift
            item[2] = item[2] + bp # end position shift
        
        self.max_dropoff = self.max_dropoff + bp
        self.read_end = self.read_end + bp
        
    def set_second_last_subtel_end_to(self, position):
        bp = position - self.subtel_list[-2][2]
        self.move_features(bp)
    
    def set_last_subtel_start_to(self, position):
        bp = position - self.subtel_list[-1][1]
        self.move_features(bp)
    
    # returns a string that can be printed to file of the data within this Object
    def table_for_gggenomes(self):
        # gggenomes columns are:
        # seq_id	start	end	strand	feat_id	width	name	geom_id	type
        lines_to_write = ""
        orientation_symbol = ""
        
        if self.orientation >= 0:
            orientation_symbol = "+"
        else:
            orientation_symbol = "-"
            
        for i in range(len(self.subtel_list)):
            lines_to_write = lines_to_write + f"{self.read_name}\t{self.subtel_list[i][1]}\t{self.subtel_list[i][2]}\t{orientation_symbol}\t{self.subtel_list[i][0]}\t{self.subtel_list[i][2]-self.subtel_list[i][1]+1}\t{self.subtel_list[i][0]}\t{self.subtel_list[i][0]}\tCDS\n"
        
        lines_to_write = lines_to_write + f"{self.read_name}\t{self.max_dropoff}\t{self.read_end}\t{orientation_symbol}\ttelomere\t{self.read_end-self.max_dropoff+1}\ttelomere\ttelomere\tCDS\n"
        
        return lines_to_write
    
    # returns the coordinates from:
    #   1 + (the end of the last subtelomere)
    #   max_dropoff value (start of the telomere)
    #   end of the read (read_end)
    def return_coord_end_of_last_subtel_to_end(self):
        if self.orientation > 1:
            return self.subtel_list[-1][2] + 1, self.max_dropoff, self.read_end
        else:
            return self.subtel_list[0][1] + 1, self.max_dropoff, self.read_end

# delete trailing dashes
def delete_trailing_dashes(alignment_string):
    # reverse the alignment_string then traverse until it is not a dash
    
    for index, letter in enumerate(alignment_string[::-1]):
        if letter == "-":
            pass
        else:
            return alignment_string[0:len(alignment_string)-index]
            break



def main():
    # Dpachys data files
    # Dpachys_porechopped.telomereScore20.Dpa1L2L_subtel.sorted.coords, Dpachys_porechopped.telomereScore20.Dpa1L2L_subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop.txt,
    # Dpachys_porechopped.telomereScore20.Dpa2R_subtel.sorted.coords, and Dpachys_porechopped.telomereScore20.Dpa2R_subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop.txt
    # have some manually curated repeat boundaries missed by nucmer
    subtel_nucmer_files = ["Dpachys_porechopped.telomereScore20.Dpa1L2L_subtel.sorted.coords",
                           "Dpachys_porechopped.telomereScore20.Dpa1R_subtel.sorted.coords",
                           "Dpachys_porechopped.telomereScore20.Dpa2R_subtel.sorted.coords"]
                           
    max_dropoff_files = ["Dpachys_porechopped.telomereScore20.Dpa1L2L_subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop.txt",
                         "Dpachys_porechopped.telomereScore20.Dpa1R_subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop.txt",
                         "Dpachys_porechopped.telomereScore20.Dpa2R_subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop.txt"]
    
    sequencing_reads = "Dpachys_porechopped.telomereScore20.fasta"
    max_leeway = 50
    
    
    orig_orientation_gggenomes_genes_file = "Dpa.genes.txt"
    orig_orientation_gggenomes_seq_file = "Dpa.seq.txt"
    centered_gggenomes_genes_file = "Dpa.centered.genes.txt"
    centered_gggenomes_seq_file = "Dpa.centered.seq.txt"
    zoomed_in_bases = "Dpa.zoomed_in_bases.gapped.txt"
    zoomed_in_bases_no_gaps = "Dpa.zoomed_in_bases.no_gaps.txt"
    
    all_subtels = "Dpa_subtelomeres.fasta"
    
    
    # extract read sequences
    reads_dict = SeqIO.to_dict(SeqIO.parse(sequencing_reads, "fasta"))
    subtel_dict = SeqIO.to_dict(SeqIO.parse(all_subtels, "fasta"))
    
    # go through the occupancy drop files
    # record read_name and max_dropoff values
    # record the subtelomere sequences, and do a Clustal alignment, then output a "zoomed in" view of the subtelomeres
    
    max_dropoff_dict = {}
    
    # write subtel sequences to a file so Clustal can be performed
    # store sequences after the subtel in a list
    subtel_file = "tmp.subtel.fasta"
    subtel_aligned_file = "tmp.subtel.aligned.fasta"
    
    
    # save 'zoomed_in' file as well with subtel sequences aligned + tel sequences
    with open(zoomed_in_bases, "w") as zoomed_fh, open(zoomed_in_bases_no_gaps, "w") as gapped_fh:
        for file in max_dropoff_files:
            # make a fasta file to collect the reads being analyzed
            fasta_file_to_collect_the_reads = file + ".reads.fasta"
            
            # define lines in max_dropoff files
            # (will sort the reads based on subtel length)
            
            max_dropoff_lines = []
                
            with open(file, "r") as fh:
                
                line = fh.readline() # this is the header
                line = fh.readline()
                while line:
                    if line[0] != "#":
                        subtel_start = int(line.split("\t")[2])
                        subtel_end = int(line.split("\t")[3])
                        last_subtel_repeat_length = 0
                        if subtel_end > subtel_start:
                            last_subtel_repeat_length = subtel_end - subtel_start + 1
                        else:
                            last_subtel_repeat_length = subtel_start - subtel_end + 1
                            
                        max_dropoff_lines.append([line, last_subtel_repeat_length])
                    else:
                        pass
                    
                    line = fh.readline()
            

            # save sequences after subtel_end as a dict of strings
            tel_dict = {}
            
            # sort the lines by the subtel length from shortest to longest
            max_dropoff_lines.sort(key = lambda x: x[1])
            
            # get ready to write to the subtel_file
            
            with open(subtel_file, "w") as subtel_fh, open(fasta_file_to_collect_the_reads, "w") as reads_fh:
                for line in max_dropoff_lines:
                    max_dropoff = 0
                    read_name = line[0].split("\t")[0]
                    subtel_start = int(line[0].split("\t")[2])
                    subtel_end = int(line[0].split("\t")[3])
                    read_length = int(line[0].split("\t")[1])
                    position_difference = int(line[0].split("\t")[8])
                    
                    if position_difference < max_leeway:
                        # CAREFUL!
                        # nucmer uses 1-based positions, and ranges are inclusive both sides
                        # BioPython uses 0-based positions, and ranges are inclusive on the left, exclusive on the right
                        #              ACGTACGTA
                        # nucmer       123456789   [1:6]   ACGTAC
                        # BioPython    012345678   [0:6]   ACGTAC   
                        if subtel_end - subtel_start > 0:
                            max_dropoff = int(line[0].split("\t")[5])
                            subtel_seq_text = str(reads_dict[read_name].seq[0:subtel_end]).lower() + str(reads_dict[read_name].seq[subtel_end:read_length]).upper()
                            subtel_seqrecord = SeqRecord(Seq(subtel_seq_text[subtel_start-1:subtel_end]),
                                                         id=read_name,
                                                         name="",
                                                         description="")
                            SeqIO.write(subtel_seqrecord, subtel_fh, "fasta")
                            tel_dict[read_name] = subtel_seq_text[subtel_end:read_length]
                            
                            SeqIO.write(SeqRecord(reads_dict[read_name].seq,
                                                  id=read_name,
                                                  name="",
                                                  description=""),
                                        reads_fh, "fasta")
                        else:
                            max_dropoff = int(line[0].split("\t")[7])
                            subtel_seq_text = str(reads_dict[read_name].seq[0:subtel_end-1]).upper() + str(reads_dict[read_name].seq[subtel_end-1:read_length]).lower()
                            subtel_seqrecord = SeqRecord(Seq(subtel_seq_text[subtel_end-1:subtel_start]).reverse_complement(),
                                                         id=read_name,
                                                         name="",
                                                         description="")
                            SeqIO.write(subtel_seqrecord, subtel_fh, "fasta")
                            tel_dict[read_name] = str(Seq(subtel_seq_text[0:subtel_end-1]).reverse_complement())
                            
                            SeqIO.write(SeqRecord(reads_dict[read_name].seq.reverse_complement(),
                                                  id=read_name,
                                                  name="",
                                                  description=""),
                                        reads_fh, "fasta")
                        
                        max_dropoff_dict[read_name] = max_dropoff
                    else:
                        pass
                # nxDipCoro1_hifi.telomereScore20.Dco1L_subtel.sorted.multiple_hits_eliminated.scanningWindow20.telOccupancyDrop
                subtel_name = file.split(".")[2]
                subtel_seqrecord = SeqRecord(Seq(str(subtel_dict[subtel_name].seq).lower()),
                                             id=subtel_name,
                                             name="",
                                             description="")
                SeqIO.write(subtel_seqrecord, subtel_fh, "fasta")
                tel_dict[subtel_name] = ""
            # now do smoe clustal on the extracted subtel sequences
            lines_to_write = f"clustalo -i {subtel_file} -o {subtel_aligned_file} --seqtype DNA --output-order=input-order --force"
            os.system(lines_to_write)
            
            # now record the subtels and telomeres into zoomed_fh
            zoomed_fh.write(f"{file}\n")
            gapped_fh.write(f"{file}\n")
            for aligned_subtel_seqrecord in SeqIO.parse(subtel_aligned_file, "fasta"):
                read_name = aligned_subtel_seqrecord.id
                lines_to_write = f"{read_name}\t{str(aligned_subtel_seqrecord.seq)}{tel_dict[read_name]}\n"
                zoomed_fh.write(lines_to_write)
                
                lines_to_write = f"{read_name}\t{delete_trailing_dashes(str(aligned_subtel_seqrecord.seq))}{tel_dict[read_name]}\n"
                gapped_fh.write(lines_to_write)
            zoomed_fh.write("\n\n\n\n\n\n\n\n\n\n")
            gapped_fh.write("\n\n\n\n\n\n\n\n\n\n")

    
    # read_information_dict is a collection of Read objects
    # key = read_name, value = Read object
    read_information_dict = {}
    read_lengths_list = []
    for file in subtel_nucmer_files: # read through the nucmer coords files
        with open(file, "r") as fh:
            # four lines of header
            line = fh.readline()
            line = fh.readline()
            line = fh.readline()
            line = fh.readline()
            
            line = fh.readline()
            while line:
                read_name = line.split("\t")[10].rstrip()
                subtel_name = line.split("\t")[9]
                subtel_start = int(line.split("\t")[2])
                subtel_end = int(line.split("\t")[3])
                
                # consider only if read_name was already added to max_dropoff_dict
                if read_name in max_dropoff_dict:
                    # if read not in the read_information_dict, add new Read object
                    if read_name not in read_information_dict:
                        read_length = int(line.split("\t")[8])
                        read_lengths_list.append(read_length)
                        read_information_dict[read_name] = Read(read_name, read_length, max_dropoff_dict[read_name], subtel_end-subtel_start)
                    else:
                        pass
                    # Add new subtel feature
                    read_information_dict[read_name].add_start_and_end(subtel_name, subtel_start, subtel_end)
                
                line = fh.readline()
    
    # print the collected subtel and telomere info to file
    with open(orig_orientation_gggenomes_genes_file, "w") as orig_genes_fh, \
         open(orig_orientation_gggenomes_seq_file, "w") as orig_seq_fh, \
         open(centered_gggenomes_genes_file, "w") as centered_genes_fh, \
         open(centered_gggenomes_seq_file, "w") as centered_seq_fh:
        genes_header = "seq_id\tstart\tend\tstrand\tfeat_id\twidth\tname\tgeom_id\ttype\n"
        orig_genes_fh.write(genes_header)
        centered_genes_fh.write(genes_header)
        
        seq_header = "seq_id\tlength\n"
        orig_seq_fh.write(seq_header)
        centered_seq_fh.write(seq_header)
        
        for read_name in max_dropoff_dict:
            orig_genes_fh.write(read_information_dict[read_name].table_for_gggenomes())
            orig_seq_fh.write(f"{read_name}\t{read_information_dict[read_name].read_length}\n")
            
            read_information_dict[read_name].reorient()
            read_information_dict[read_name].set_last_subtel_start_to(max(read_lengths_list))
            centered_genes_fh.write(read_information_dict[read_name].table_for_gggenomes())
            
            centered_seq_fh.write(read_name + "\t" + str(read_information_dict[read_name].read_end) + "\n")
    
    os.system("rm tmp.*")

if __name__ == "__main__":
    main()