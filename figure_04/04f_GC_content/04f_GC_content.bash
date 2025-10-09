###############################################################################
# GC content along chromosomes
###############################################################################
#
# D. coronatus
# Genome fasta = "Dcor.genomic.fasta"
python3 fraction_GC_along_length.py -f Dcor.genomic.fasta -x scaffold_position -y fraction_GC -w 100000 -i 100000

# D. pachys
# GEnome fasta = "Dpac.genomic.fasta"
python3 fraction_GC_along_length.py -f Dpac.genomic.fasta -x scaffold_position -y fraction_GC -w 100000 -i 100000