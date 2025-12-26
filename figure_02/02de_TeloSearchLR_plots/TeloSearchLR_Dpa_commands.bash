# Use TeloSearchLR to find the telomeric repeat motif in the Porechopped Dpachys ONT reads (Dpachys_porechopped.fasta)
TeloSearchLR.py -f Dpachys_porechopped.fasta -k 4 -K 20 -t 1000 -n 1000 -m 1 -M 100 -c 10

# We find a mirrored occupancy for the pattern TAAGGGTAAGGC, the 38th most frequently found
# pattern in the terminal 1000 bp. 
# Upon closer inspection, there appears to be some alternation between the TAAGGG and TAAGGC motifs.
# Fetch occupancy data for TAAGGG as well using the TideHunter table generated in the last step
# (TeloSearchLR.k4.K20.Dpachys_porechopped.TideHunterTable.txt)
TeloSearchLR.py -s TAAGGG -T TeloSearchLR.k4.K20.Dpachys_porechopped.TideHunterTable.txt -f Dpachys_porechopped.fasta -n 1000


# Proceed with graphing the occupancies using R
