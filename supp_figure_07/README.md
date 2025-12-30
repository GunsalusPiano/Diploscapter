# Chromosome center and arm characteristics for P. exspectatus, showing typical Rhabditid chromosome arm-center differentiation.

*Pristionchus exspectatus* genome assembly and annotations can be found at [Pristionchus.org/download](http://pristionchus.org/download/) (linked under "*P. exspectatus* from Yoshida *et al. Nature Evol. Ecology.* 2023). I prepared the data for analysis by having a genomic FASTA of only the nuclear chromosomes ("Pexspectatus_nuclear_genome.fasta"). I used the annotations ("exspectatus_canonical.gff3") and the protein FASTA ("exspectatus_nr_proteins.fa") as provided on Pristionchus.org. The *P. exspectatus* Hi-C mapping data were provided by Kohta Yoshida and used with permission in the manuscript. Most of the scripts used below can be found under **Fig. 4**.
- Fraction repeats:
  ```{bash}
  python3 fraction_RepeatScout_repeats_along_length_v2.py -f Pexspectatus_nuclear_genome.fasta -w 100000 -i 100000 -p Pexsp
  ```
- Fraction exonic
  ```{bash}
  python3 fraction_exonic_DNA_along_length_v2.py -f Pexspectatus_nuclear_genome.fasta -g exspectatus_canonical.gff3 -p Pexsp -w 100000 -i 10000
  ```
- Fraction GC
  ```{bash}
  python3 fraction_GC_along_length_v2.py -f Pexspectatus_nuclear_genome.fasta -p Pexsp -w 100000 -i 100000
  ```
- BUSCO search
  ```{bash}
  busco -i exspectatus_nr_proteins.fa -l nematoda_odb10 -o Pexspectatus_protein_busco -m prot -c 14
  ```
  then use the python script below to turn the protein busco table into the genomic equivalent.
  ```{bash}
  python3 BUSCO_protein_full_table_to_BUSCO_genomic_full_table_Pexsp_Yoshida_version.py
  ```
  
The data were plotted by running the R code in ```supp_figure_2024-12-28.md```, with the resulting plots below.
<img width="1170" height="931" alt="image" src="https://github.com/user-attachments/assets/d7c34ab3-4aa3-4745-96cf-c94df863ce13" />

