##########################################################################################
# Plotting the fraction of restriction enzyme sites along the lengths of Dco chromosomes #
##########################################################################################
#
# An Arima HiC experiment was done on D. coronatus chromosomes to aid the assembly scaffolding.
# Arima HiC restriction enzymes cut at GATC and GANTC. We used the script below to identify and
# calculate the occurrence of these RE sequences to see if the HiC results corresponded with
# RE cut frequencies (they do not).
#
# First, make a file "Arima_HiC_restriction_sites.txt" with every recognition sequence (and reverse
# complement) separated by a return character.
#
# Then, run "fraction_small_multiseq_along_length.py" to count up all the RE sites in intervals
# of 100,000 bps.
#
# Arguments:
#    -f: The genomic scaffolds FASTA file (="Dcor.genomic.fasta")
#    -s: The text file (= "Arima_HiC_restriction_sites.txt") with all the RE sites
#    -p: output file prefix
#    -w: window (interval) size to calculate the RE site occupancy (=100000)
#    -i: increment size (=100000)
#    -x: x-axis name
#    -y: y-axis name
python3 fraction_small_multiseq_along_length.py -f Dcor.genomic.fasta -s Arima_HiC_restriction_sites.txt -p Dcor -w 100000 -i 100000 -x scaffold_position -y fraction_Arima_RE

##########################################################################################
# Plotting the fraction of restriction enzyme sites along the lengths of Dpa chromosomes #
##########################################################################################
#
# D. pachys assembly was scaffolded using PoreC data - which only uses NlaIII.
# run "fraction_small_palindrome_along_length.py" to count up all the RE sites in intervals
# of 100,000 bps.
#
# Arguments:
#    -f: The genomic scaffolds FASTA file (="Dpac.genomic.fasta")
#    -s: RE recognition sequence (for this, NlaIII site GATC)
#    -p: output file prefix
#    -w: window (interval) size to calculate the RE site occupancy (=100000)
#    -i: increment size (=100000)
#    -x: x-axis name
#    -y: y-axis name
python3 fraction_small_palindrome_along_length.py -f Dpac.genomic.fasta -s GATC -p Dpac -w 100000 -i 100000 -x scaffold_position -y fraction_NlaIII