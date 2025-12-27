# Calculate the fraction of exonic DNA along the length of chromosomes
Use the following commands and the corresponding GTFs (available at Zenodo) to get the fraction of exonic DNA every 100,000 bps.

```python
python3 fraction_exonic_DNA_along_length.py -f ${GENOMIC_FASTA} -g ${GTF_ANNOTATIONS} -p Dpac -x scaffold_position -y fraction_exonic -w 100000 -i 100000
python3 fraction_exonic_DNA_along_length.py -f ${GENOMIC_FASTA} -g ${GTF_ANNOTATIONS} -p Dcor -x scaffold_position -y fraction_exonic -w 100000 -i 100000
```
