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

HFC_fill_alpha <- 0.2

DpaA_cen1s <- 10000000
DpaA_cen1e <- 13000000
DpaA_cen2s <- 25750000
DpaA_cen2e <- 28750000
DpaA_cen3s <- 32300000
DpaA_cen3e <- 35300000
DpaA_cen4s <- 37540000
DpaA_cen4e <- 40540000
DpaA_cen5s <- 45750000
DpaA_cen5e <- 48750000
DpaA_cen6s <- 59500000
DpaA_cen6e <- 62500000
DpaA_cen7s <- 72000000
DpaA_cen7e <- 75000000

DpaB_cen2s <- 25300000
DpaB_cen2e <- 28300000
DpaB_cen3s <- 33250000
DpaB_cen3e <- 36250000
DpaB_cen4s <- 39750000
DpaB_cen4e <- 42750000
DpaB_cen5s <- 48500000
DpaB_cen5e <- 51500000
DpaB_cen6s <- 61500000
DpaB_cen6e <- 64500000
DpaB_cen7s <- 75000000
DpaB_cen7e <- 78000000

DcoA_cen1s <- 9000000
DcoA_cen1e <- 12000000
DcoA_cen2s <- 24500000
DcoA_cen2e <- 27500000
DcoA_cen3s <- 30500000
DcoA_cen3e <- 33500000
DcoA_cen4s <- 36500000
DcoA_cen4e <- 39500000
DcoA_cen5s <- 44500000
DcoA_cen5e <- 47500000
DcoA_cen6s <- 56500000
DcoA_cen6e <- 59500000
DcoA_cen7s <- 69000000
DcoA_cen7e <- 72000000
DcoA_cen8s <- 75500000
DcoA_cen8e <- 78500000

DcoB_cen1s <- 9500000
DcoB_cen1e <- 12500000
DcoB_cen2s <- 25500000
DcoB_cen2e <- 28500000
DcoB_cen3s <- 33500000
DcoB_cen3e <- 36500000
DcoB_cen4s <- 39000000
DcoB_cen4e <- 42000000
DcoB_cen5s <- 47750000
DcoB_cen5e <- 50750000
DcoB_cen6s <- 60500000
DcoB_cen6e <- 63500000
DcoB_cen7s <- 74000000
DcoB_cen7e <- 77000000


input_filename <- "all_vs_all.DpaAB_DcoAB.comparison.percentage_id_along_length.window_100000.sample_freq_50000.reject_threshold_0.7.txt"

df_all_vs_all <- read_tsv(input_filename)
```

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 20616 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): reference, query
    ## dbl (4): coordinate, evaluated, matched, percentage_id
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
graph_point_grid_colorScheme3 <- ggplot(df_all_vs_all) +
  geom_point(mapping = aes(x = coordinate, y = percentage_id, color = query), size = 0.2, alpha = 0.3) +
  facet_grid( reference ~ . ) +
  theme_bw() + #, "#377eb8", "#4daf4a", "#ff7f00"
  scale_color_manual(values=c("#4daf4a", "#377eb8", "#ff7f00", "#e41a1c")) +
  scale_x_continuous(labels = scales::comma) +
  coord_cartesian(ylim = c(87, 100), xlim = c(0, 90000000)) +
  theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=8), axis.text.y=element_text(size=8), axis.title.y=element_text(size=8)) + 
  theme(strip.text = element_text(size = 8)) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = 0, xmax = 17839851, ymin = 0, ymax = 88), alpha = 0.6, fill="grey", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = DpaA_cen1s, xmax = DpaA_cen1e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = DpaA_cen2s, xmax = DpaA_cen2e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = DpaA_cen4s, xmax = DpaA_cen4e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = DpaA_cen5s, xmax = DpaA_cen5e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = DpaA_cen6s, xmax = DpaA_cen6e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaA"), aes(xmin = DpaA_cen7s, xmax = DpaA_cen7e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = 0, xmax = 17839851, ymin = 0, ymax = 88), alpha = 0.6, fill="grey", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = DpaA_cen1s, xmax = DpaA_cen1e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = DpaB_cen2s, xmax = DpaB_cen2e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = DpaB_cen4s, xmax = DpaB_cen4e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = DpaB_cen5s, xmax = DpaB_cen5e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = DpaB_cen6s, xmax = DpaB_cen6e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DpaB"), aes(xmin = DpaB_cen7s, xmax = DpaB_cen7e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoA"), aes(xmin = DcoA_cen1s, xmax = DcoA_cen1e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoA"), aes(xmin = DcoA_cen2s, xmax = DcoA_cen2e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoA"), aes(xmin = DcoA_cen4s, xmax = DcoA_cen4e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoA"), aes(xmin = DcoA_cen5s, xmax = DcoA_cen5e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoA"), aes(xmin = DcoA_cen6s, xmax = DcoA_cen6e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoA"), aes(xmin = DcoA_cen7s, xmax = DcoA_cen7e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoB"), aes(xmin = DcoB_cen1s, xmax = DcoB_cen1e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoB"), aes(xmin = DcoB_cen2s, xmax = DcoB_cen2e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoB"), aes(xmin = DcoB_cen4s, xmax = DcoB_cen4e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoB"), aes(xmin = DcoB_cen5s, xmax = DcoB_cen5e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoB"), aes(xmin = DcoB_cen6s, xmax = DcoB_cen6e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(reference = "DcoB"), aes(xmin = DcoB_cen7s, xmax = DcoB_cen7e, ymin = 80, ymax = 110), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  guides(color = guide_legend(override.aes = list(size=1, alpha=0.8)))

graph_point_grid_colorScheme3
```

    ## Warning: Removed 1969 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](nucmer_percentage_id_along_length_visualization_2025-12_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave('all_vs_all.DpaAB_DcoAB.comparison.percentage_id_along_length.window_100000.sample_freq_50000.reject_threshold_0.7.points_only.16x11cm.changed_colors3.pdf', graph_point_grid_colorScheme3, width = 16, height = 11, units = 'cm')
```

    ## Warning: Removed 1969 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or
by pressing *Ctrl+Alt+I*.

When you save the notebook, an HTML file containing the code and output
will be saved alongside it (click the *Preview* button or press
*Ctrl+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the
editor. Consequently, unlike *Knit*, *Preview* does not run any R code
chunks. Instead, the output of the chunk when it was last run in the
editor is displayed.
