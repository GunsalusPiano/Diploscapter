# run the BUSCO search on the proteins from individual chromosomes, rather than the diploid genome, as below.
busco -i Dpa-canu-het00c-YaHS-v202304.chrom_scaffolds.scaffold_1.proteins.faa -l nematoda_odb10 -o DpaA_BUSCO -m protein -c 14
busco -i Dpa-canu-het00c-YaHS-v202304.chrom_scaffolds.scaffold_2.proteins.faa -l nematoda_odb10 -o DpaB_BUSCO -m protein -c 14
busco -i nxDipCoro1_1.curated_primary.SUPER_2.proteins.faa -l nematoda_odb10 -o DcoA_BUSCO -m protein -c 14
busco -i nxDipCoro1_1.curated_primary.SUPER_1.proteins.faa -l nematoda_odb10 -o DcoB_BUSCO -m protein -c 14

# take the full_table.tsv from each run then prefix with species and chrom scaffold name, eg. Dpac.scaffold_1.full_table.tsv
# where
#      scaffold_1, scaffold_2 = DpaA, DpaB
#           SUPER_1, SUPER_2 = DcoB, DcoA

# then look up the Augustus gtf file for the locations of the encoding genes, then record these down in a new file.
# python script BUSCO_protein_full_table_to_BUSCO_genomic_full_table_v2.py
#
#    -b   busco full_table.tsv file path/name
#    -g   gtf file from BRAKER2/Augustus
#    -p   new genomic equivalent BUSCO table's file prefix
python3 BUSCO_protein_full_table_to_BUSCO_genomic_full_table_v2.py -b Dpac.scaffold_1.full_table.tsv -g Dpa-canu-het00c-YaHS-v202304.chrom_scaffolds.gtf -p Dpac.scaffold_1

python3 BUSCO_protein_full_table_to_BUSCO_genomic_full_table_v2.py -b Dpac.scaffold_2.full_table.tsv -g Dpa-canu-het00c-YaHS-v202304.chrom_scaffolds.gtf -p Dpac.scaffold_2
# because these annotations from AUGUSTUS did not contain the scaffold_2 homozygous region, change the 'scaffold_1' into 'scaffold_2' in the tsv
sed -i 's/scaffold_1/scaffold_2/g' Dpac.scaffold_2.genomic_busco_equiv_of_Dpac.scaffold_2.full_table.tsv

python3 BUSCO_protein_full_table_to_BUSCO_genomic_full_table_v2.py -b Dcor.SUPER_1.full_table.tsv -g nxDipCoro1_1.curated_primary.gtf -p Dcor.SUPER_1

python3 BUSCO_protein_full_table_to_BUSCO_genomic_full_table_v2.py -b Dcor.SUPER_2.full_table.tsv -g nxDipCoro1_1.curated_primary.gtf -p Dcor.SUPER_2

# the genomic equivalent BUSCO table is now ready for Nigon painting 