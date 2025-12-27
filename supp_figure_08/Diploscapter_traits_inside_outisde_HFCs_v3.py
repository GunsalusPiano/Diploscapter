################################################
# Diploscapter_traits_inside_outisde_HFCs.py
# 
# Takes a chromosome traits file (eg. GC content)
# and a file defining the HFCs, outputs two files
# - one with traits inside the HFCs and one out-
# side. This is so one could do statistics on 
# these traits.
################################################

import sys
import getopt

version_number = "0.0.1"

help_text = """Diploscapter_traits_inside_outisde_HFCs.py
                
Version: {0}   Contact: Dr. George Chung (gc95@nyu.edu)
                
Usage:   python {1} [options]
                
Options:
  -f    STR    tab-delimited file describing the feature name (1st column), start coordinate (2nd column) and end coordinate (3rd column)
  -d    STR    data file describing the coordinate (1st column) and fraction of a feature (eg. %GC, 2nd column)
  -o    STR    output file with a new column added separating data INSIDE and OUTSIDE the features.\n
"""

def filename_extension_remover(filename):
    filename_parts = filename.split(".")
    filename_parts.pop(-1)
    
    return ".".join(filename_parts)

# HFC coordinates to 
def get_features_start_and_end_coords(features_file):
    # features_list is a nested list containing 3-member lists:
    # [[feature_name1, start1, end1], [feature_name2, start2, end2], ...etc]
    features_list = []
        
    
    with open(features_file, "r") as features_fh:
        # first line is the header
        features_fh.readline()
        
        # grab the second line data
        line = features_fh.readline()
        
        while line:
            line_items = line.rstrip().split("\t")
            # the line structure is [HFC_name]  [start]  [end]
            features_list.append([line_items[0], int(line_items[1]), int(line_items[2])])
            line = features_fh.readline()
    return features_list

def add_feature_column(data_file, features_list, non_feature_name, output_file):
    with open(data_file, "r") as data_fh, open(output_file, "w") as output_fh:
        # transfer the first line to the new file
        line = data_fh.readline().rstrip() + "\tfeature\n"
        output_fh.write(line)
        
        line = data_fh.readline()
        
        while line:
            line_items = line.split("\t")
            current_line_coord = int(line_items[0])
            
            feature_name = non_feature_name
            
            # feature (eg. HFC i)
            # check to see if the current coordinate is a feature
            # features_list is a nested list containing 3-member lists:
            # [[feature_name1, start1, end1], [feature_name2, start2, end2], ...etc]
            for feature in features_list:
                start = feature[1]
                end = feature[2]
                
                # if the current coordinate falls within the range of a feature
                # update the feature name
                # change the 
                if current_line_coord >= start and current_line_coord <= end:
                    feature_name = feature[0]
                    break
                else:
                    pass
                    
            line_to_write = line.rstrip() + "\t" + feature_name + "\n"
            output_fh.write(line_to_write)
            
            line = data_fh.readline()

def main(argv):
    '''# Start of the script
    feature_file_HFC_dict = {"Dcor.08.SUPER_1.repeats.window100000.increment100000.txt": "DcoB_HFCs.txt",
                             "Dcor.08.SUPER_2.repeats.window100000.increment100000.txt": "DcoA_HFCs.txt",
                             "Dcor.SUPER_1.fraction_Arima_RE.window100000.increment100000.txt": "DcoB_HFCs.txt",
                             "Dcor.SUPER_1.fraction_GC.window100000.increment100000.txt": "DcoB_HFCs.txt",
                             "Dcor.SUPER_1.fraction_exonic.window100000.increment100000.txt": "DcoB_HFCs.txt",
                             "Dcor.SUPER_2.fraction_Arima_RE.window100000.increment100000.txt": "DcoA_HFCs.txt",
                             "Dcor.SUPER_2.fraction_GC.window100000.increment100000.txt": "DcoA_HFCs.txt",
                             "Dcor.SUPER_2.fraction_exonic.window100000.increment100000.txt": "DcoA_HFCs.txt",
                             "Dpac.08.scaffold_1.repeats.window100000.increment100000.txt": "DpaA_HFCs.txt",
                             "Dpac.08.scaffold_2.repeats.window100000.increment100000.txt": "DpaB_HFCs.txt",
                             "Dpac.scaffold_1.fraction_GC.window100000.increment100000.txt": "DpaA_HFCs.txt",
                             "Dpac.scaffold_1.fraction_NlaIII.window100000.increment100000.txt": "DpaA_HFCs.txt",
                             "Dpac.scaffold_1.fraction_exonic.window100000.increment100000.txt": "DpaA_HFCs.txt",
                             "Dpac.scaffold_2.fraction_GC.window100000.increment100000.txt": "DpaB_HFCs.txt",
                             "Dpac.scaffold_2.fraction_NlaIII.window100000.increment100000.txt": "DpaB_HFCs.txt",
                             "Dpac.scaffold_2.fraction_exonic.window100000.increment100000.txt": "DpaB_HFCs.txt"}
    
    for feature_file in feature_file_HFC_dict:
        feature_list = get_features_start_and_end_coords(feature_file_HFC_dict[feature_file])
        non_feature_name = "NA"
        output_file = filename_extension_remover(feature_file) + ".by_HFCs.txt"
        add_feature_column(feature_file, feature_listsas, non_feature_name, output_file)'''
    
    # build a big file containing all the data
    RE_files = {"DcoA":"Dcor.SUPER_2.fraction_Arima_RE.window100000.increment100000.txt",
                "DcoB": "Dcor.SUPER_1.fraction_Arima_RE.window100000.increment100000.txt",
                "DpaA": "Dpac.scaffold_1.fraction_NlaIII.window100000.increment100000.txt",
                "DpaB": "Dpac.scaffold_2.fraction_NlaIII.window100000.increment100000.txt"}

    rep_files = {"DcoA": "Dcor.08.SUPER_2.repeats.window100000.increment100000.txt",
                 "DcoB": "Dcor.08.SUPER_1.repeats.window100000.increment100000.txt",
                 "DpaA": "Dpac.08.scaffold_1.repeats.window100000.increment100000.txt",
                 "DpaB": "Dpac.08.scaffold_2.repeats.window100000.increment100000.txt"}
                 
    exon_files = {"DcoA": "Dcor.SUPER_2.fraction_exonic.window100000.increment100000.txt",
                  "DcoB": "Dcor.SUPER_1.fraction_exonic.window100000.increment100000.txt",
                  "DpaA": "Dpac.scaffold_1.fraction_exonic.window100000.increment100000.txt",
                  "DpaB": "Dpac.scaffold_2.fraction_exonic.window100000.increment100000.txt"}
                  
    GC_files = {"DcoA": "Dcor.SUPER_2.fraction_GC.window100000.increment100000.txt",
                "DcoB": "Dcor.SUPER_1.fraction_GC.window100000.increment100000.txt",
                "DpaA": "Dpac.scaffold_1.fraction_GC.window100000.increment100000.txt",
                "DpaB": "Dpac.scaffold_2.fraction_GC.window100000.increment100000.txt"}
    
    HFC_files = {"DcoA": "DcoA_HFCs.txt",
                 "DcoB": "DcoB_HFCs.txt",
                 "DpaA": "DpaA_HFCs.txt",
                 "DpaB": "DpaB_HFCs.txt"}
                
    dip_chromosomes = ["DcoA", "DcoB", "DpaA", "DpaB"]
    
    all_data_file = "all_data.txt"
    
    with open(all_data_file, "w") as all_fh:
        line_to_write = "chromosome\tcoordinate\tvalue\tdata_type\tHFC\n"
        all_fh.write(line_to_write)
        
    for chromosome in dip_chromosomes:
        with open(RE_files[chromosome], "r") as RE_fh, open(rep_files[chromosome], "r") as rep_fh, open(exon_files[chromosome], "r") as exon_fh, \
        open(GC_files[chromosome], "r") as GC_fh, open(all_data_file, "a") as all_fh:
            # clear the first line of each file
            RE_fh.readline()
            rep_fh.readline()
            exon_fh.readline()
            GC_fh.readline()
            
            # proceed to extract data from each file and transfer to the new data file "all_data.txt"
            line = RE_fh.readline()
            
            # grab the features list (list of HFC coordinates)
            features_list = get_features_start_and_end_coords(HFC_files[chromosome])
            
            
            while line:
                coordinate = int(line.split("\t")[0])
                value = line.rstrip().split("\t")[1]
                HFC = "non-HFC"
                
                # feature is the 3-member list [feature, start, stop]
                for feature in features_list:
                    start = feature[1]
                    end = feature[2]
                    
                    if coordinate >= start and coordinate <= end:
                        HFC = feature[0]
                        break
                    else:
                        pass
                
                line_to_write = chromosome + "\t" + str(coordinate) + "\t" + value + "\t" + "RE_density\t" + HFC + "\n"
                all_fh.write(line_to_write)
                
                # move on to repeats
                line = rep_fh.readline()
                value = ""
                if len(line.split("\t")) == 4:
                    value = line.rstrip().split("\t")[3]
                else:
                    pass
                line_to_write = chromosome + "\t" + str(coordinate) + "\t" + value + "\t" + "repeat_density\t" + HFC + "\n"
                all_fh.write(line_to_write)
                
                # move on to exon density
                line = exon_fh.readline()
                value = ""
                if len(line.split("\t")) == 2:
                    value = line.rstrip().split("\t")[1]
                else:
                    pass
                line_to_write = chromosome + "\t" + str(coordinate) + "\t" + value + "\t" + "exon_density\t" + HFC + "\n"
                all_fh.write(line_to_write)
                
                # move to GC_density
                line = GC_fh.readline()
                value = ""
                if len(line.split("\t")) == 2:
                    value = line.rstrip().split("\t")[1]
                else:
                    pass
                line_to_write = chromosome + "\t" + str(coordinate) + "\t" + value + "\t" + "GC_content\t" + HFC + "\n"
                all_fh.write(line_to_write)
                
                line = RE_fh.readline()

            
    
if __name__ == "__main__":
    main(sys.argv)