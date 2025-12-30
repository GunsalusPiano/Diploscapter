# Chromosome center and arm characteristics for P. exspectatus, showing typical Rhabditid chromosome arm-center differentiation.

*Pristionchus exspectatus* genome assembly and annotations can be found at [Pristionchus.org/download](http://pristionchus.org/download/) (linked under "*P. exspectatus* from Yoshida *et al. Nature Evol. Ecology.* 2023). I prepared the data for analysis by having a genomic FASTA of only the nuclear chromosomes ("Pexspectatus_nuclear_genome.fasta"). I used the annotations ("exspectatus_canonical.gff3") and the protein FASTA ("exspectatus_nr_proteins.fa") as provided on Pristionchus.org. The *P. exspectatus* Hi-C mapping data were provided by Kohta Yoshida and used with permission in the manuscript.
- Fraction repeats:
  ```{bash}
  python3 fraction_RepeatScout_repeats_along_length_v2.py -f Pexspectatus_nuclear_genome.fasta -w 100000 -i 100000 -p Pexsp
  ```
- Fraction exonic ```{bash}python3 fraction_exonic_DNA_along_length_v2.py -f Pexspectatus_nuclear_genome.fasta -g exspectatus_canonical.gff3 -p Pexsp -w 100000 -i 10000```
- Fraction GC ```{bash}python3 fraction_GC_along_length_v2.py -f Pexspectatus_nuclear_genome.fasta -p Pexsp -w 100000 -i 100000```
- BUSCO search ```{bash}busco -i exspectatus_nr_proteins.fa -l nematoda_odb10 -o Pexspectatus_protein_busco -m prot -c 14```

