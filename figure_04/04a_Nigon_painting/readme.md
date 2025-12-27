# BUSCO searches and Nigon painting

For BUSCO searches ([Manni & al 2021. _Current Protocols_](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.323)), since high-quality annotations based on RNAseq data already exist (see Zenodo) for both _D. pachys_ and _D. coronatus_, _ab initio_ gene prediction by metaeuk or Augustus was not used. (Metaeuk, the default for BUSCO searches, was found to be unreliable in gene prediction for BUSCO searches.) Instead, BUSCO searches were done in the protein mode, and the protein BUSCO tables were turned into equivalent genomic BUSCO tables using the gene annotations.

```BUSCO_analysis.bash``` contains the commands to run BUSCO searches, and the commands to run ```BUSCO_protein_full_table_to_BUSCO_genomic_full_table_v2.py``` to turn protein BUSCO table to its genomic equivalent (for Nigon painting).
