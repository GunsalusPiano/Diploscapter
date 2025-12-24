**Figure 4b:** HiC and PoreC principal components for *Diploscapter* contact data.

The ```.bash``` file contains the run commands to get the PCs for contact data using hicexplorer (https://hicexplorer.readthedocs.io/en/latest/). For *D. pachys*, this involves an additional step: as the homozygous region of the genome is represented only once in the assembly, the read depth in this region is doubled. I include the ```Dpa-canu-het00c-YaHS-v202304-hethomosep_PoreC_correction.py``` script to re-map reads mapping to the homozygous region to either DpaA or DpaB.
