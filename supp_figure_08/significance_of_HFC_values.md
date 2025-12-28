R Notebook
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(ggpubr)
```

    ## 
    ## Attaching package: 'ggpubr'
    ## 
    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

``` r
# RE data
chrom_traits_df = read_tsv("Diploscapter.RE.data.txt")
```

    ## Rows: 3438 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): chromosome, HFC
    ## dbl (2): coordinate, RE_density
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
RE_density_dotplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=RE_density, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_dotplot(binaxis = "y",
                                   stackdir = "center",
                                   fill = "black",
                                   dotsize = 0.5,
                                   binwidth = 0.0005,
                                   stackratio = 0.3,
                                   alpha = 0.7,
                                   stroke = 0) +
                                theme_bw()
                               


RE_density_dotplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# RE boxplot with stats
compare_means(RE_density ~ HFC, data = chrom_traits_df)
```

    ## # A tibble: 1 × 8
    ##   .y.        group1  group2         p    p.adj p.format p.signif method  
    ##   <chr>      <chr>   <chr>      <dbl>    <dbl> <chr>    <chr>    <chr>   
    ## 1 RE_density non-HFC HFC    0.0000946 0.000095 9.5e-05  ****     Wilcoxon

``` r
RE_density_boxplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=RE_density, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_boxplot(fill = "darkorchid4", alpha = 0.5, lwd = 0.2, outlier.size = 0.5, width = 0.5) +
                                theme_bw() +
                                stat_compare_means(label="p.format") +
                                theme(axis.title.x = element_blank(),
                                      axis.text = element_text(size = 8),
                                      axis.text.x = element_blank())

RE_density_boxplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# repeat data
chrom_traits_df = read_tsv("Diploscapter.repeats.data.txt")
```

    ## Rows: 3438 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): chromosome, HFC
    ## dbl (2): coordinate, repeat_density
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
repeat_density_dotplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=repeat_density, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_dotplot(binaxis = "y",
                                   stackdir = "center",
                                   fill = "black",
                                   dotsize = 0.8,
                                   binwidth = 0.01,
                                   stackratio = 0.4,
                                   alpha = 0.7,
                                   stroke = 0) +
                                theme_bw()
                               


repeat_density_dotplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# repeat density boxplot with stats
repeat_density_boxplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=repeat_density, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_boxplot(fill = "darkgreen", alpha = 0.5, lwd = 0.2, outlier.size = 0.5, width = 0.5) +
                                theme_bw() +
                                stat_compare_means(label="p.format") +
                                theme(axis.title.x = element_blank(),
                                      strip.text = element_blank(),
                                      axis.text = element_text(size = 8),
                                      axis.text.x = element_blank())

repeat_density_boxplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# exon density data
chrom_traits_df = read_tsv("Diploscapter.exons.data.txt")
```

    ## Rows: 3438 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): chromosome, HFC
    ## dbl (2): coordinate, exon_density
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
exon_density_dotplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=exon_density, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_dotplot(binaxis = "y",
                                   stackdir = "center",
                                   fill = "black",
                                   dotsize = 0.5,
                                   binwidth = 0.01,
                                   stackratio = 0.4,
                                   alpha = 0.7,
                                   stroke = 0) +
                                theme_bw()
                               


exon_density_dotplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# exon density boxplot with stats
exon_density_boxplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=exon_density, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_boxplot(fill = "darkorange2", alpha = 0.5, lwd = 0.2, outlier.size = 0.5, width = 0.5) +
                                theme_bw() +
                                stat_compare_means(label="p.format") +
                                theme(axis.title.x=element_blank(),
                                      strip.text = element_blank(),
                                      axis.text = element_text(size = 8),
                                      axis.text.x = element_blank())

exon_density_boxplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# GC content data
chrom_traits_df = read_tsv("Diploscapter.GC_content.data.txt")
```

    ## Rows: 3438 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): chromosome, HFC
    ## dbl (2): coordinate, GC_content
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
GC_content_dotplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=GC_content, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_dotplot(binaxis = "y",
                                   stackdir = "center",
                                   fill = "black",
                                   dotsize = 0.5,
                                   binwidth = 0.002,
                                   stackratio = 0.4,
                                   alpha = 0.7,
                                   stroke = 0) +
                                theme_bw()
                               


GC_content_dotplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# GC content boxplot with stats
GC_content_boxplot <- ggplot(chrom_traits_df, aes(x=factor(HFC, levels=c("non-HFC","HFC")), y=GC_content, scales = "free")) +
                                facet_grid(. ~ factor(chromosome, levels=c("DpaA","DpaB","DcoA","DcoB"))) +
                                geom_boxplot(fill = "blue", alpha = 0.5, lwd = 0.2, outlier.size = 0.5, width = 0.5) +
                                theme_bw() +
                                stat_compare_means(label="p.format", size=3, hjust=-0.4, vjust = 6) +
                                theme(axis.title.x = element_blank(),
                                      strip.text = element_blank(),
                                      axis.text = element_text(size = 8))

GC_content_boxplot
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
stats_grid <-plot_grid(RE_density_boxplot,
                      repeat_density_boxplot,
                      exon_density_boxplot,
                      GC_content_boxplot,
                      ncol = 1, nrow = 4, align = "v", axis = "lr", rel_heights = c(1.1,1,1,1.1))

stats_grid
```

![](significance_of_HFC_values_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave('genomic_features_stats.v5.pdf', stats_grid, width = 12, height = 16, units = 'cm')
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or
by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output
will be saved alongside it (click the *Preview* button or press
*Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the
editor. Consequently, unlike *Knit*, *Preview* does not run any R code
chunks. Instead, the output of the chunk when it was last run in the
editor is displayed.
