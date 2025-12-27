# chromosome_wide_MSA_of_BUSCOs_2024-01-02_update.py
# attempt to align many BUSCOs from DpaA, DpaB, DcoA, DcoB and C. elegans to figure out the evolutionary relationships between them



from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os


def Cbriggsae_FASTA_eliminate_multiple_isoforms(multiple_isoform_fasta, single_isoform_fasta):
	with open(multiple_isoform_fasta, "r") as multiple_isoform_fasta_fh, open(single_isoform_fasta, "w") as single_isoform_fasta_fh:
		for current_seq in SeqIO.parse(multiple_isoform_fasta_fh, "fasta"):
			if current_seq.id[-1].isdigit() or current_seq.id[-1] == "a":
				SeqIO.write(current_seq, single_isoform_fasta_fh, "fasta")
			
			else:
				pass

def Celegans_FASTA_eliminate_multiple_isoforms(multiple_isoform_fasta, single_isoform_fasta):
	with open(multiple_isoform_fasta, "r") as multiple_isoform_fasta_fh, open(single_isoform_fasta, "w") as single_isoform_fasta_fh:
		for current_seq in SeqIO.parse(multiple_isoform_fasta_fh, "fasta"):
			if current_seq.id[-1].isdigit() or current_seq.id[-1] == "a":
				SeqIO.write(current_seq, single_isoform_fasta_fh, "fasta")
			
			else:
				pass

def eliminate_non_complete_BUSCOs(BUSCO_full_table_filename, list_of_complete_BUSCOs):
	with open(BUSCO_full_table_filename, "r") as BUSCO_full_table_fh:
		line = BUSCO_full_table_fh.readline()
		
		while line:
			if line[0] == "#":
				pass
			else:
				line_items = line.split("\t")
				BUSCO_id = line_items[0]
				BUSCO_status = line_items[1]
				
				# keep the BUSCO only if it's already in the list AND "Complete" in C. briggsae
				# 4 conditions:
				# BUSCO_status is "Complete" or otherwise   x   BUSCO_id is in list_of_complete_BUSCOs or otherwise
				# only when it's already in the list_of_complete_BUSCOs but BUSCO_status != "Complete" do we need to do something
				# ie. remove it from the list.
				
				if BUSCO_status != "Complete" and BUSCO_id in list_of_complete_BUSCOs:
					list_of_complete_BUSCOs.remove(BUSCO_id)
				else:
					pass
			line = BUSCO_full_table_fh.readline()
	return list_of_complete_BUSCOs

# returns a dictionary from which you can access SeqRecord by looking up BUSCO_id
def make_BUSCO_id_to_seqRecord_dictionary(single_isoform_fasta, BUSCO_full_table_filename, list_of_complete_BUSCOs):
	gene_name_to_BUSCO_id_dictionary = {}
	with open(BUSCO_full_table_filename, "r") as BUSCO_full_table_fh:
		line = BUSCO_full_table_fh.readline()
		
		while line:
			if line[0] == "#":
				pass
			else:
				# retrieve gene names corresponding to BUSCO ids
				line_items = line.split("\t")
				BUSCO_id = line_items[0]
				
				if BUSCO_id in list_of_complete_BUSCOs:
					gene_name = line_items[2]
					gene_name_to_BUSCO_id_dictionary[gene_name] = BUSCO_id
				else:
					pass
			
			line = BUSCO_full_table_fh.readline()
		
	# step through the FASTA and retrieve sequences that are in list_of_complete_BUSCOs
	
	BUSCO_id_to_seqRecord_dictionary = {}
	for current_seq in SeqIO.parse(single_isoform_fasta, "fasta"):
		if current_seq.id in gene_name_to_BUSCO_id_dictionary:
			BUSCO_id_to_seqRecord_dictionary[gene_name_to_BUSCO_id_dictionary[current_seq.id]] = current_seq
		else:
			pass
	
	return BUSCO_id_to_seqRecord_dictionary
	

def turn_protein_FASTA_into_CDS_FASTA(protein_seqRecord, transcript_dictionary):
	final_aligned_transcript_sequence = ""
	
	gene_name = protein_seqRecord.description.split()[1]
	
	unaligned_transcript_sequence = transcript_dictionary[gene_name]
	
	#aa_counter = 0
	nt_counter = 0
	
	for amino_acid in protein_seqRecord.seq:
		if amino_acid != "-":
			final_aligned_transcript_sequence = final_aligned_transcript_sequence + unaligned_transcript_sequence[nt_counter:nt_counter+3]
			nt_counter = nt_counter + 3
		else:
			final_aligned_transcript_sequence = final_aligned_transcript_sequence + "---"
	
	record = SeqRecord(seq = final_aligned_transcript_sequence,
						id = protein_seqRecord.id,
						name = "",
						description = protein_seqRecord.description)
	return record
	

def repopulate_list_of_complete_BUSCOs(list_of_complete_BUSCOs_file):
    list_of_complete_BUSCOs = []
    with open(list_of_complete_BUSCOs_file, "r") as list_of_complete_BUSCOs_fh:
        line = list_of_complete_BUSCOs_fh.readline()
        
        while line:
            if line[0] != "#":
                list_of_complete_BUSCOs.append(line.rstrip())
            else:
                pass
            line = list_of_complete_BUSCOs_fh.readline()
    return list_of_complete_BUSCOs
    
def file_extension_remover(filename):
    new_filename = ""
    filename_split = filename.split(".")
    for split_part in range(len(filename_split)-1):
        new_filename = new_filename + filename_split[split_part] + "."
    return new_filename[0:-1]

def main():
    
    # We will not use the C. briggsae BUSCOs for this version
    # eliminate the multiple isoforms per gene
    
    multiple_isoform_fasta = "c_elegans.PRJNA13758.WS285.protein.fa"
    single_isoform_fasta = "Cel_proteins.faa"
    
    Celegans_FASTA_eliminate_multiple_isoforms(multiple_isoform_fasta, single_isoform_fasta)
    
    multiple_isoform_fasta = "c_elegans.PRJNA13758.WS285.CDS_transcripts.fa"
    single_isoform_fasta = "Cel_transcripts.fna"
    
    Celegans_FASTA_eliminate_multiple_isoforms(multiple_isoform_fasta, single_isoform_fasta)
    
    
    # did 2 BUSCO runs - one with Celegans data and one with Cbriggsae data
    # command = busco -i Cel_proteins.faa -l nematoda_odb10 -o Cel_proteins_busco -m prot -c 14
    # command = busco -i Cbr_proteins.faa -l nematoda_odb10 -o Cbr_proteins_busco -m prot -c 14
    
    
    # walk through "full_table.tsv" for every species and only keep entries that are Busco complete (ie. single full copy)
    
    list_of_complete_BUSCOs = []
    
    # first run through the list: Let's look through the C. elegans list of complete BUSCOs
    
    BUSCO_full_table_filename = "Cel_full_table.tsv"
    
    with open(BUSCO_full_table_filename, "r") as BUSCO_full_table_fh:
        # If the line begins with "#", ignore
        # else, record the busco ID only if the second column [1] reads "Complete"
        line = BUSCO_full_table_fh.readline()
        
        while line:
            if line[0] == "#":
                pass
            else:
                line_items = line.split("\t")
                BUSCO_id = line_items[0]
                BUSCO_status = line_items[1]
                if BUSCO_status == "Complete":
                    list_of_complete_BUSCOs.append(BUSCO_id)
                else:
                    pass
            line = BUSCO_full_table_fh.readline()
    
    # Now, list_of_complete_BUSCOs should only contain BUSCOs listed as "complete" in C. elegans
    # Go through this process again, but with C. briggsae
    # This time, just consider those BUSCOs already in the list, and consider if it's "Complete" in Diploscapter chromosomes
    # Starting with DpaA
    BUSCO_full_table_filename = "DpaA_full_table.tsv"
    list_of_complete_BUSCOs = eliminate_non_complete_BUSCOs(BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # and DpaB
    BUSCO_full_table_filename = "DpaB_full_table.tsv"
    list_of_complete_BUSCOs = eliminate_non_complete_BUSCOs(BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # and DcoA
    BUSCO_full_table_filename = "DcoA_full_table.tsv"
    list_of_complete_BUSCOs = eliminate_non_complete_BUSCOs(BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # and DcoB
    BUSCO_full_table_filename = "DcoB_full_table.tsv"
    list_of_complete_BUSCOs = eliminate_non_complete_BUSCOs(BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    
    # now list_of_complete_BUSCOs should contain only "Complete" BUSCOs from all 4 species (6 haploid genomes)
    # save this list to file.
    with open("list_of_complete_BUSCOs.txt", "w") as complete_BUSCO_fh:
        complete_BUSCO_fh.write("# number of complete BUSCOs in Cel, Cbr, DpaA, DpaB, DcoA, DcoB: ")
        complete_BUSCO_fh.write(str(len(list_of_complete_BUSCOs)) + "\n")
        for BUSCO in list_of_complete_BUSCOs:
            complete_BUSCO_fh.write(BUSCO + "\n")
    
    
    # Go back to the .tsv files, retrieve the actual sequences and put in 6 dictionaries.
    # for these dictionaries, you can access the SeqRecord object by looking up BUSCO_id
    # C. elegans
    single_isoform_fasta = "Cel_proteins.faa"
    BUSCO_full_table_filename = "Cel_full_table.tsv"
    Cel_complete_BUSCO_proteins_dictionary = make_BUSCO_id_to_seqRecord_dictionary(single_isoform_fasta, BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # DpaA
    single_isoform_fasta = "DpaA_proteins.faa"
    BUSCO_full_table_filename = "DpaA_full_table.tsv"
    DpaA_complete_BUSCO_proteins_dictionary = make_BUSCO_id_to_seqRecord_dictionary(single_isoform_fasta, BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # DpaB
    single_isoform_fasta = "DpaB_proteins.faa"
    BUSCO_full_table_filename = "DpaB_full_table.tsv"
    DpaB_complete_BUSCO_proteins_dictionary = make_BUSCO_id_to_seqRecord_dictionary(single_isoform_fasta, BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # DcoA
    single_isoform_fasta = "DcoA_proteins.faa"
    BUSCO_full_table_filename = "DcoA_full_table.tsv"
    DcoA_complete_BUSCO_proteins_dictionary = make_BUSCO_id_to_seqRecord_dictionary(single_isoform_fasta, BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    # DcoB
    single_isoform_fasta = "DcoB_proteins.faa"
    BUSCO_full_table_filename = "DcoB_full_table.tsv"
    DcoB_complete_BUSCO_proteins_dictionary = make_BUSCO_id_to_seqRecord_dictionary(single_isoform_fasta, BUSCO_full_table_filename, list_of_complete_BUSCOs)
    
    
    
    # save BUSCO orthologues in individual FASTA files so they can be MUSCLE-aligned later
    
    for BUSCO_id in list_of_complete_BUSCOs:
        FASTA_filename = "unaligned." + BUSCO_id + ".AA.orthologues.faa"
        
        with open(FASTA_filename, "w") as FASTA_fh:
            # write in this order: Cel, Cbr, DpaA, DpaB, DcoA, DcoB
            
            # write the Cel orthologue
            # remove the terminal "*"
            sequence = Cel_complete_BUSCO_proteins_dictionary[BUSCO_id].seq
            if sequence[-1] == "*":
                sequence = sequence[0:-1]
            else:
                pass
				
            record = SeqRecord(
                        seq = sequence,
                        id = "Cel",
                        name = "",
                        description = Cel_complete_BUSCO_proteins_dictionary[BUSCO_id].id)
            SeqIO.write(record, FASTA_fh, "fasta")
            
            # write the DpaA orthologue, but only if the gene name number is > g.r.17444.t1 (heterozygous region)
            # or > g.f.17673.t1
            if DpaA_complete_BUSCO_proteins_dictionary[BUSCO_id].id[2] == "f" and \
                int(DpaA_complete_BUSCO_proteins_dictionary[BUSCO_id].id[4:-3]) > 17673 or \
                DpaA_complete_BUSCO_proteins_dictionary[BUSCO_id].id[2] == "r" and \
                int(DpaA_complete_BUSCO_proteins_dictionary[BUSCO_id].id[4:-3]) > 17444:
                
                # remove the terminal "*"
                sequence = DpaA_complete_BUSCO_proteins_dictionary[BUSCO_id].seq
                if sequence[-1] == "*":
                    sequence = sequence[0:-1]
                else:
                    pass
                record = SeqRecord(
                        seq = sequence,
                        id = "DpaA",
                        name = "",
                        description = DpaA_complete_BUSCO_proteins_dictionary[BUSCO_id].id)
                SeqIO.write(record, FASTA_fh, "fasta")
            else:
                pass
            
            # write the DpaB orthologue
            # remove the terminal "*"
            sequence = DpaB_complete_BUSCO_proteins_dictionary[BUSCO_id].seq
            if sequence[-1] == "*":
                sequence = sequence[0:-1]
            else:
                pass
            
            record = SeqRecord(
                        seq = sequence,
                        id = "DpaB",
                        name = "",
                        description = DpaB_complete_BUSCO_proteins_dictionary[BUSCO_id].id)
            SeqIO.write(record, FASTA_fh, "fasta")
            
            # write the DcoA orthologue
            # remove the terminal "*"
            sequence = DcoA_complete_BUSCO_proteins_dictionary[BUSCO_id].seq
            if sequence[-1] == "*":
                sequence = sequence[0:-1]
            else:
                pass
            record = SeqRecord(
                        seq = sequence,
                        id = "DcoA",
                        name = "",
                        description = DcoA_complete_BUSCO_proteins_dictionary[BUSCO_id].id)
            SeqIO.write(record, FASTA_fh, "fasta")
            
            # write the DcoB orthologue
            # remove the terminal "*"
            sequence = DcoB_complete_BUSCO_proteins_dictionary[BUSCO_id].seq
            if sequence[-1] == "*":
                sequence = sequence[0:-1]
            else:
                pass
            record = SeqRecord(
                        seq = sequence,
                        id = "DcoB",
                        name = "",
                        description = DcoB_complete_BUSCO_proteins_dictionary[BUSCO_id].id)
            SeqIO.write(record, FASTA_fh, "fasta")
    
    
    # attempt to muscle-align through all the multi-fastas
    # muscle v5.1.0 (https://www.drive5.com/muscle/)
    # make sure muscle is in the system PATH
    
    
    for BUSCO_id in list_of_complete_BUSCOs:
        unaligned_FASTA_filename = "unaligned." + BUSCO_id + ".AA.orthologues.faa"
        aligned_FASTA_filename = "aligned." + BUSCO_id + ".AA.orthologues.faa"
        
        line_to_write = "muscle -align " + unaligned_FASTA_filename + " -output " + aligned_FASTA_filename
        
        os.system(line_to_write)
    
    # now turn amino acid alignment to nucleotide alignment
    # generate transcript dictionaries
    # key = gene name in respective genomes/proteomes
    # value = string representing the nucleotide sequences of the transcripts
    
    # C. elegans transcripts dictionary
    transcripts_filename = "Cel_transcripts.fna"
    Cel_transcripts_dict = {}	
    for seqRecord in SeqIO.parse(transcripts_filename, "fasta"):
        Cel_transcripts_dict[seqRecord.id] = seqRecord.seq
    
    
    # D. pachys transcripts dictionary
    transcripts_filename = "Dpa_transcripts.fna"
    Dpa_transcripts_dict ={}
    for seqRecord in SeqIO.parse(transcripts_filename, "fasta"):
        Dpa_transcripts_dict[seqRecord.id] = seqRecord.seq
    
    
    # D. coronatus transcripts dictionary
    transcripts_filename = "Dco_transcripts.fna"
    Dco_transcripts_dict ={}
    for seqRecord in SeqIO.parse(transcripts_filename, "fasta"):
        Dco_transcripts_dict[seqRecord.id] = seqRecord.seq
    
    list_of_complete_BUSCOs = repopulate_list_of_complete_BUSCOs("list_of_complete_BUSCOs.txt")
       
    for BUSCO_id in list_of_complete_BUSCOs:
        protein_aligned_FASTA_filename = "aligned." + BUSCO_id + ".AA.orthologues.faa"
        
        # open a file for nucleotide alignments
        nucleotide_aligned_FASTA_filename = "aligned." + BUSCO_id + ".NT.orthologues.fna"
        
        with open(nucleotide_aligned_FASTA_filename, "w") as nucleotide_aligned_FASTA_fh:
            for protein_seqRecord in SeqIO.parse(protein_aligned_FASTA_filename, "fasta"):
                if protein_seqRecord.id == "Cel":
                    SeqIO.write(turn_protein_FASTA_into_CDS_FASTA(protein_seqRecord, Cel_transcripts_dict), nucleotide_aligned_FASTA_fh, "fasta")
                elif protein_seqRecord.id == "Cbr":
                    SeqIO.write(turn_protein_FASTA_into_CDS_FASTA(protein_seqRecord, Cbr_transcripts_dict), nucleotide_aligned_FASTA_fh, "fasta")
                elif protein_seqRecord.id == "DpaA" or protein_seqRecord.id == "DpaB":
                    SeqIO.write(turn_protein_FASTA_into_CDS_FASTA(protein_seqRecord, Dpa_transcripts_dict), nucleotide_aligned_FASTA_fh, "fasta")
                elif protein_seqRecord.id == "DcoA" or protein_seqRecord.id == "DcoB":
                    SeqIO.write(turn_protein_FASTA_into_CDS_FASTA(protein_seqRecord, Dco_transcripts_dict), nucleotide_aligned_FASTA_fh, "fasta")
                else:
                    pass
    
    # for organizational purposes, move the unaligned AA sequences to another folder
    line_to_write = "mkdir unaligned_BUSCO_proteins"
    os.system(line_to_write)
    
    line_to_write = "mv unaligned.* unaligned_BUSCO_proteins"
    os.system(line_to_write)
    
    # for organizational purposes, move the aligned AA sequences to another folder
    line_to_write = "mkdir aligned_BUSCO_proteins"
    os.system(line_to_write)
    
    line_to_write = "mv *AA.orthologues.faa aligned_BUSCO_proteins"
    os.system(line_to_write)
    
    
    ##################################################################################################################################################
    # trim alignments and concatenate the alignments from the Dpa homozygous region and alignments from the Dpa het region 
    ##################################################################################################################################################
    # at the same time, separate alignments corresponding to the homozygous and heterozygous regions in DpaA
    # This is so we can concatenate the alignments later on
    list_of_complete_quality_BUSCOs_Dpa_hom = []
    list_of_complete_quality_BUSCOs_Dpa_het = []
    
    # re-populate list_of_complete_BUSCOs (earlier block of code isn't re-run)
    list_of_complete_BUSCOs = repopulate_list_of_complete_BUSCOs("list_of_complete_BUSCOs.txt")
    
    # Use trimAl to trim the alignments and eliminate all positions with gaps
    for BUSCO_id in list_of_complete_BUSCOs:
        # name of each alignment file follows this format
        nucleotide_aligned_FASTA_filename = "aligned." + BUSCO_id + ".NT.orthologues.fna"
        
        # test if DpaA is in the alignment - if it is, the BUSCO is in the Dpa heterozygous region
        has_DpaA = False
        for current_seq in SeqIO.parse(nucleotide_aligned_FASTA_filename, "fasta"):
            if current_seq.id == "DpaA":
                has_DpaA = True
            else:
                pass
        
        # trimAl output filename below
        trimAl_output_filename = "trimAl.aligned." + BUSCO_id + ".NT.orthologues.fna"
        
        if has_DpaA:
            list_of_complete_quality_BUSCOs_Dpa_het.append(trimAl_output_filename)
        else:
            list_of_complete_quality_BUSCOs_Dpa_hom.append(trimAl_output_filename)

        # run trimAl
        # the -gt 1 option from trimAl is to remove any column with any gaps
        string_to_write = "trimal -in " + nucleotide_aligned_FASTA_filename + " -out " + trimAl_output_filename + " -gt 1"
        os.system(string_to_write)
        
    # write the het/hom BUSCO alignment lists to file
    filename = "02.info.Dpa_het_alignments.txt"
    with open(filename, "w") as fh:
        for alignment_filename in list_of_complete_quality_BUSCOs_Dpa_het:
            fh.write(alignment_filename + "\n")
    filename = "01.info.Dpa_hom_alignments.txt"
    with open(filename, "w") as fh:
        for alignment_filename in list_of_complete_quality_BUSCOs_Dpa_hom:
            fh.write(alignment_filename + "\n")
    
    # concatenate the alignments.
    # deal with the HOMOZYGOUS region first
    # we will also record down the total alignment length for RAxML
    
    # first deal with the sequence (text) concatenation
    Cel_sequence = ""
    DpaB_sequence = ""
    DcoA_sequence = ""
    DcoB_sequence = ""
    
    total_aligned_seq_length = 0
    
    for alignment_file in list_of_complete_quality_BUSCOs_Dpa_hom:
        for seqRecord in SeqIO.parse(alignment_file, "fasta"):
            if seqRecord.id == "Cel":
                Cel_sequence = Cel_sequence + seqRecord.seq
            elif seqRecord.id == "DpaB":
                DpaB_sequence = DpaB_sequence + seqRecord.seq
            elif seqRecord.id == "DcoA":
                DcoA_sequence = DcoA_sequence + seqRecord.seq
            elif seqRecord.id == "DcoB":
                DcoB_sequence = DcoB_sequence + seqRecord.seq
            else:
                pass
    
    total_aligned_seq_length = len(Cel_sequence)
    
    # turn the text into SeqRecord objects
    Cel_seqRecord = SeqRecord(seq = Cel_sequence,
                                id = "Cel",
                                name = "",
                                description = "")
	
    DpaB_seqRecord = SeqRecord(seq = DpaB_sequence,
                                id = "DpaB",
                                name = "",
                                description = "")
								
    DcoA_seqRecord = SeqRecord(seq = DcoA_sequence,
                                id = "DcoA",
                                name = "",
                                description = "")
	
    DcoB_seqRecord = SeqRecord(seq = DcoB_sequence,
                                id = "DcoB",
                                name = "",
                                description = "")
    
    # write the SeqRecords to file
    filename = "01.concatenated.alignments_mapping_to_Dpa_hom_region.fna"
    with open(filename, "w") as fh:
        SeqIO.write(Cel_seqRecord, fh, "fasta")
        SeqIO.write(DpaB_seqRecord, fh, "fasta")
        SeqIO.write(DcoA_seqRecord, fh, "fasta")
        SeqIO.write(DcoB_seqRecord, fh, "fasta")
	
    # record down the alignment length
    filename = "01.concatenated.alignments_mapping_to_Dpa_hom_region.alignmentlength"
    with open(filename, "w") as fh:
        line_to_write = str(total_aligned_seq_length)
        fh.write(line_to_write)
    
    # write the partition file for modeltest-ng (later)
    filename = file_extension_remover(filename) + ".modeltest-ng.partitions"
    with open(filename, "w") as fh:
        line_to_write = "DNA, CDS_1 = 1-" + str(total_aligned_seq_length) + "\\3\n"
        fh.write(line_to_write)
        line_to_write = "DNA, CDS_2 = 2-" + str(total_aligned_seq_length) + "\\3\n"
        fh.write(line_to_write)
        line_to_write = "DNA, CDS_3 = 3-" + str(total_aligned_seq_length) + "\\3\n"
        fh.write(line_to_write)
    
    # deal with the heterozygous region now
    # we will also record down the total alignment length for RAxML
    # deal with the sequence (text) concatenation
    Cel_sequence = ""
    DpaA_sequence = ""
    DpaB_sequence = ""
    DcoA_sequence = ""
    DcoB_sequence = ""
    
    total_aligned_seq_length = 0
    
    for alignment_file in list_of_complete_quality_BUSCOs_Dpa_het:
        for seqRecord in SeqIO.parse(alignment_file, "fasta"):
            if seqRecord.id == "Cel":
                Cel_sequence = Cel_sequence + seqRecord.seq
            elif seqRecord.id == "DpaA":
                DpaA_sequence = DpaA_sequence + seqRecord.seq
            elif seqRecord.id == "DpaB":
                DpaB_sequence = DpaB_sequence + seqRecord.seq
            elif seqRecord.id == "DcoA":
                DcoA_sequence = DcoA_sequence + seqRecord.seq
            elif seqRecord.id == "DcoB":
                DcoB_sequence = DcoB_sequence + seqRecord.seq
            else:
                pass
    
    total_aligned_seq_length = len(Cel_sequence)
    
    # turn the text into SeqRecord objects
    Cel_seqRecord = SeqRecord(seq = Cel_sequence,
                                id = "Cel",
                                name = "",
                                description = "")
	
    DpaA_seqRecord = SeqRecord(seq = DpaA_sequence,
                                id = "DpaA",
                                name = "",
                                description = "")
    
    DpaB_seqRecord = SeqRecord(seq = DpaB_sequence,
                                id = "DpaB",
                                name = "",
                                description = "")
								
    DcoA_seqRecord = SeqRecord(seq = DcoA_sequence,
                                id = "DcoA",
                                name = "",
                                description = "")
	
    DcoB_seqRecord = SeqRecord(seq = DcoB_sequence,
                                id = "DcoB",
                                name = "",
                                description = "")
    
    # write the SeqRecords to file
    filename = "02.concatenated.alignments_mapping_to_Dpa_het_region.fna"
    with open(filename, "w") as fh:
        SeqIO.write(Cel_seqRecord, fh, "fasta")
        SeqIO.write(DpaA_seqRecord, fh, "fasta")
        SeqIO.write(DpaB_seqRecord, fh, "fasta")
        SeqIO.write(DcoA_seqRecord, fh, "fasta")
        SeqIO.write(DcoB_seqRecord, fh, "fasta")
	
    # record down the alignment length
    filename = "02.concatenated.alignments_mapping_to_Dpa_het_region.alignmentlength"
    with open(filename, "w") as fh:
        line_to_write = str(total_aligned_seq_length)
        fh.write(line_to_write)
        
    # write the partition file for modeltest-ng (later)
    filename = file_extension_remover(filename) + ".modeltest-ng.partitions"
    with open(filename, "w") as fh:
        line_to_write = "DNA, CDS_1 = 1-" + str(total_aligned_seq_length) + "\\3\n"
        fh.write(line_to_write)
        line_to_write = "DNA, CDS_2 = 2-" + str(total_aligned_seq_length) + "\\3\n"
        fh.write(line_to_write)
        line_to_write = "DNA, CDS_3 = 3-" + str(total_aligned_seq_length) + "\\3\n"
        fh.write(line_to_write)

    # Do a final bit of organization.
    # Move all untrimmed alignments to a new folder
    line_to_write = "mkdir untrimmed_aligned_BUSCO_CDSs"
    os.system(line_to_write)
    line_to_write = "mv aligned.* untrimmed_aligned_BUSCO_CDSs"
    os.system(line_to_write)
    
    # Move all trimmed alignments to a new folder
    line_to_write = "mkdir trimmed_aligned_BUSCO_CDSs"
    os.system(line_to_write)
    line_to_write = "mv trimAl.* trimmed_aligned_BUSCO_CDSs"
    os.system(line_to_write)
    
    # Now the concatenated alignments should be ready for use in RAxML-ng
    # They are:
    #          concatenated.alignments_mapping_to_Dpa_hom_region.fna
    #          concatenated.alignments_mapping_to_Dpa_het_region.fna
        
    ##################################################################################################################################################
    # The alignments are now exactly of lengths %3 = 0, because these are CDSs aligned using their translations. This means each of the three
    # positions of the codons can have a different substitution model.
    # use modeltest-ng to find the most suitable substitution model.
    ##################################################################################################################################################
    
    # region corresponding to the DpaA homozygous
    concatenated_alignment = "01.concatenated.alignments_mapping_to_Dpa_hom_region.fna"
    modeltest_partitions = file_extension_remover(concatenated_alignment) + ".modeltest-ng.partitions"
    modeltest_output = file_extension_remover(concatenated_alignment) + ".modeltest-ng.results"
    
    line_to_write = "modeltest-ng -i " + concatenated_alignment + " -q " + modeltest_partitions + " -o " + modeltest_output + " --force"
    os.system(line_to_write)
    
    # Hom region partitions should use:
    #   CDS_1: GTR+I
    #   CDS_2: GTR+G4
    #   CDS_3: GTR+I
    
    
    # region corresponding to the DpaA het
    concatenated_alignment = "02.concatenated.alignments_mapping_to_Dpa_het_region.fna"
    modeltest_partitions = file_extension_remover(concatenated_alignment) + ".modeltest-ng.partitions"
    modeltest_output = file_extension_remover(concatenated_alignment) + ".modeltest-ng.results"
    
    line_to_write = "modeltest-ng -i " + concatenated_alignment + " -q " + modeltest_partitions + " -o " + modeltest_output + " --force"
    os.system(line_to_write)
    
    # Het region partions should use:
    #   CDS_1: GTR+I
    #   CDS_2: GTR+I+G4
    #   CDS_3: GTR+I+G4
    
    ##############################################################################
    # write the partition files for RAxML-ng and run RAxML-NG (trial01)
    ##############################################################################
    
    # hom region
    concatenated_alignment = "01.concatenated.alignments_mapping_to_Dpa_hom_region.fna"
    alignment_length_file = file_extension_remover(concatenated_alignment) + ".alignmentlength"
    raxmlng_partitions = file_extension_remover(concatenated_alignment) + ".raxml-ng.partitions"
    with open(raxmlng_partitions, "w") as raxmlng_partitions_fh, open(alignment_length_file, "r") as alignment_length_fh:
        alignment_length = alignment_length_fh.readline().rstrip()
        line_to_write = "GTR+I, CDS_1=1-" + alignment_length + "\\3\n"
        raxmlng_partitions_fh.write(line_to_write)
        line_to_write = "GTR+G4, CDS_2=2-" + alignment_length + "\\3\n"
        raxmlng_partitions_fh.write(line_to_write)
        line_to_write = "GTR+I, CDS_3=3-" + alignment_length + "\\3\n"
        raxmlng_partitions_fh.write(line_to_write)
    
    # run raxml-ng
    trial_name = "01.trial01.Dpa_hom"
    line_to_write = "raxml-ng --all --msa " + concatenated_alignment + " --model " + raxmlng_partitions + " --prefix " + trial_name
    os.system(line_to_write)
    
    # het region
    concatenated_alignment = "02.concatenated.alignments_mapping_to_Dpa_het_region.fna"
    alignment_length_file = file_extension_remover(concatenated_alignment) + ".alignmentlength"
    raxmlng_partitions = file_extension_remover(concatenated_alignment) + ".raxml-ng.partitions"
    with open(raxmlng_partitions, "w") as raxmlng_partitions_fh, open(alignment_length_file, "r") as alignment_length_fh:
        alignment_length = alignment_length_fh.readline().rstrip()
        line_to_write = "GTR+I, CDS_1=1-" + alignment_length + "\\3\n"
        raxmlng_partitions_fh.write(line_to_write)
        line_to_write = "GTR+I+G4, CDS_2=2-" + alignment_length + "\\3\n"
        raxmlng_partitions_fh.write(line_to_write)
        line_to_write = "GTR+I+G4, CDS_3=3-" + alignment_length + "\\3\n"
        raxmlng_partitions_fh.write(line_to_write)
    
    # run raxml-ng
    trial_name = "02.trial02.Dpa_het"
    line_to_write = "raxml-ng --all --msa " + concatenated_alignment + " --model " + raxmlng_partitions + " --prefix " + trial_name
    os.system(line_to_write)
    
    

if __name__ == "__main__":
    main()