# Quantification of traits inside and outside center-like, high frequency contact (HFC) regions

We use the script ```Diploscapter_traits_inside_outisde_HFCs_v4.py``` to extract measurements of
- RE site density
- repeat density
- exonic DNA density
- GC content

from **Fig. 4**, and classify these as either within or outside HFC regions. The HFC coordinates are defined in ```DpaA_HFCs.txt```, ```DpaB_HFCs.txt```, ```DcoA_HFCs.txt``` and ```DcoB_HFCs.txt```. To classify the data as inside or outside the HFCs, run
```{bash}
python3 Diploscapter_traits_inside_outisde_HFCs_v4.py
````

The results are stored in ```Diploscapter.exons.data.txt```, ```Diploscapter.GC_content.data.txt```, ```Diploscapter.RE.data.txt```, and ```Diploscapter.repeats.data.txt```. The R code in ```significance_of_HFC_values.md``` generates box-and-whisker plots from these results, and determines the statistical significance.

<img width="642" height="851" alt="image" src="https://github.com/user-attachments/assets/48cbf9f4-84de-48ff-ac6a-fef0adcaa5a6" />

