# Counting the restriction enzyme sites along the *Diploscapter* chromosomes

The purpose of looking for restriction sites used by Arima Hi-C and Pore-C is to determine if the chromosome contacts seen are an artefact of the non-random distribution of restriction enzyme sites. For *D. pachys*, the PoreC protocol used NlaIII (recognition site CATG), while the *D. coronatus* used Arima's proprietary Hi-C cocktail.
- The search for NlaIII sites on DpaA and DpaB was done using ```fraction_small_palindrome_along_length.py```, which finds palindromes in sequence FASTAs.
- The search for Arima's RE sites on DcoA and DcoB was done using ```fraction_small_multiseq_along_length.py```, with the specific restriction sites listed in ```Arima_HiC_restriction_sites.txt```.

The commands to run these scripts are in the ```04c_RE_sites.bash``` file. Plotting of the resulting data was done with the R ggplot2 code in ```genomic_features_visualizations.md``` in the main **Fig. 4** folder.
