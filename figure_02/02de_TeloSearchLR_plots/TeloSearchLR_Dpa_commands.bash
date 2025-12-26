# Use TeloSearchLR to find the telomeric repeat motif in the Porechopped Dpachys ONT reads (Dpachys_porechopped.fasta)
TeloSearchLR.py -f Dpachys_porechopped.fasta -k 4 -K 20 -t 1000 -n 1000 -m 1 -M 100 -c 10

# We find a mirrored occupancy for the pattern TAAGGGTAAGGC and its reverse complement, as the 38th most frequently found 
# pattern in the terminal 1000 bp of the reads. Upon closer inspection, this pattern appears to be the result of alternating 
# between the TAAGGG and TAAGGC repeat motifs.

# Fetch the occupancy data of TAAGGG and TAAGGGTAAGGC, using the TideHunter table generated in the last step
TeloSearchLR.py -s TAAGGGTAAGGC -T TeloSearchLR.k4.K20.Dpachys_porechopped.TideHunterTable.txt -f Dpachys_porechopped.fasta -n 1000
TeloSearchLR.py -s TAAGGG -T TeloSearchLR.k4.K20.Dpachys_porechopped.TideHunterTable.txt -f Dpachys_porechopped.fasta -n 1000

# Do the same with Dcoronatus HiFi reads (nxDipCoro1_hifi.fasta)
TeloSearchLR.py -f nxDipCoro1_hifi.fasta -k 4 -K 20 -t 1000 -n 1000 -m 1 -M 100 -c 10
TeloSearchLR.py -s TAAGGGTAAGGC -T TeloSearchLR.k4.K20.nxDipCoro1_hifi.TideHunterTable.txt -f nxDipCoro1_hifi.fasta -n 1000
TeloSearchLR.py -s TAAGGG -T TeloSearchLR.k4.K20.nxDipCoro1_hifi.TideHunterTable.txt -f nxDipCoro1_hifi.fasta -n 1000

# Rename the txt files with the occupancy data with the species abbreviations.
# eg. TAAGGG.occupancy.n1000.all.bar.txt > Dpa.TAAGGG.occupancy.n1000.all.bar.txt
# The occupancy data are ready to be graphed using R ggplot2.

