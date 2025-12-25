R Notebook
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

``` r
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.3     ✔ tibble    3.2.1
    ## ✔ purrr     1.0.2     ✔ tidyr     1.3.1
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
axis_text_size <- 8


selected_telomeres_df = read_tsv("all_nematode_telomeres.headerless.txt", col_names=c("accession", "species", "strain", "description", "read_name", "telomere_length"),
                            col_types=c("c", "c", "c", "c", "c", "i"), comment = "#")


telomere_dotplot_horizontal <- ggplot(selected_telomeres_df, aes(description, telomere_length)) +
                                #geom_boxplot(width = 0.5,
                                #   linewidth = 0.2) +
                                geom_dotplot(binaxis = "y",
                                   stackdir = "center",
                                   fill = "black",
                                   dotsize = 1,
                                   binwidth = 50,
                                   stackratio = 0.7,
                                   alpha = 0.7,
                                   stroke = 0) +
                                theme_bw() +   
                                xlab("experiment") +
                                ylab("double-stranded telomere length (bps)") +
                                scale_x_discrete(limits = c("C. elegans VC2010/PD1074 re-sequencing (Yoshimura 2019)",
                                                            "C. elegans CB4856 re-sequencing (Kim 2019)",
                                                            "C. briggsae QX1410 sequencing (Stevens 2022)",
                                                            "D. coronatus PacBio HiFi genomic library (this work)",
                                                            "D. pachys ONT R9.4 genomic library (this work)"),
                                                 labels = c("C. elegans VC2010",
                                                            "C. elegans CB4856",
                                                            "C. briggsae QX1410",
                                                            "D. coronatus PDL0010",
                                                            "D. pachys PF1309")) +
                                scale_y_continuous(limits = c(0,8000)) +
                                #stat_n_text() +
                                #stat_median_iqr_text() +
                                #geom_text(data = telomere_stats, aes(description, number_of_reads, label = number_of_reads),
                                #          position = position_dodge(width = 0.8), size = 3, vjust = -0.5) +
                                #geom_dotplot(binaxis='y', stackdir='center', dotsize = .5) +
                                #scale_fill_discrete(breaks=c("D. pachys",
                                #                            "C. elegans",
                                #                            "C. briggsae",
                                #                            "P. pacificus",
                                #                            "O. tipulae")) +
                                #scale_fill_brewer(palette="Pastel1")  +
                                theme(legend.key.size = unit(0.2, 'cm')) +
                                coord_flip() +
                                theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                                theme(axis.title.x = element_text(size = axis_text_size), axis.title.y = element_text(size = axis_text_size)) +
                                theme(legend.text=element_text(size=axis_text_size), legend.title=element_text(size=axis_text_size))
                                
                                #scale_fill_brewer(palette="BuPu")


telomere_dotplot_horizontal
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`stat_bindot()`).

![](telomere_lengths_violin_plot_2025-12-20_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("telomere_lengths_dot_box_plots_horizontal_2025-12.pdf", telomere_dotplot_horizontal, width = 9, height = 4, units = "cm")
```

    ## Warning: Removed 9 rows containing missing values or values outside the scale range
    ## (`stat_bindot()`).

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or
by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output
will be saved alongside it (click the *Preview* button or press
*Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the
editor. Consequently, unlike *Knit*, *Preview* does not run any R code
chunks. Instead, the output of the chunk when it was last run in the
editor is displayed.
