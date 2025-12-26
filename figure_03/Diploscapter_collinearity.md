R Notebook
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(readr)
library(magrittr)
library(GenomicRanges)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'GenomicRanges'

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     subtract

``` r
library(knitr)
library(ggplot2)
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

    ## The following object is masked from 'package:magrittr':
    ## 
    ##     extract

``` r
margins <- 0

#############################################################################################################
# DpaB vs DpaA
#############################################################################################################

chromcompare_01_DpaBvsDpaA = read_tsv("DpaB_vs_DpaA.delta-filter.options_rql0000.coords")
```

    ## Rows: 1935 Columns: 8

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): rid, qid, strand
    ## dbl (5): rs, re, qs, qe, error
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
chromcompare_01_DpaBvsDpaA_plot <- ggplot(chromcompare_01_DpaBvsDpaA, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          #coord_cartesian(xlim=c(0,90000000), ylim=c(0,90000000)) +
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #geom_rect(aes(xmin = 1, xmax = 17839851, ymin = 1, ymax = 17839851), alpha = 0.2, fill = "grey", color = NA, inherit.aes = FALSE) +
          annotate("rect", xmin = 1, xmax = 17839851, ymin = 1, ymax = 17839851, alpha = .2,fill = "grey") +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm"))

#############################################################################################################
# DcoA vs DpaA
#############################################################################################################

chromcompare_02_DcoAvsDpaA = read_tsv("DcoA_vs_DpaA.delta-filter.options_rql0000.coords")
```

    ## Rows: 2187 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): rid, qid, strand
    ## dbl (5): rs, re, qs, qe, error
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
chromcompare_02_DcoAvsDpaA_plot <- ggplot(chromcompare_02_DcoAvsDpaA, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm"))

#############################################################################################################
# DcoB vs DpaA
#############################################################################################################

chromcompare_03_DcoBvsDpaA = read_tsv("DcoB_vs_DpaA.delta-filter.options_rql0000.coords")
```

    ## Rows: 2504 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): rid, qid, strand
    ## dbl (5): rs, re, qs, qe, error
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
chromcompare_03_DcoBvsDpaA_plot <- ggplot(chromcompare_03_DcoBvsDpaA, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm"))


#############################################################################################################
# DcoA vs DpaB
#############################################################################################################

chromcompare_04_DcoAvsDpaB = read_tsv("DcoA_vs_DpaB.delta-filter.options_rql0000.coords")
```

    ## Rows: 2412 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): rid, qid, strand
    ## dbl (5): rs, re, qs, qe, error
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
chromcompare_04_DcoAvsDpaB_plot <- ggplot(chromcompare_04_DcoAvsDpaB, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm"))

#############################################################################################################
# DcoB vs DpaB
#############################################################################################################

chromcompare_05_DcoBvsDpaB = read_tsv("DcoB_vs_DpaB.delta-filter.options_rql0000.coords")
```

    ## Rows: 2524 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): rid, qid, strand
    ## dbl (5): rs, re, qs, qe, error
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
chromcompare_05_DcoBvsDpaB_plot <- ggplot(chromcompare_05_DcoBvsDpaB, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm"))


#############################################################################################################
# DcoB vs DcoA
#############################################################################################################

chromcompare_06_DcoBvsDcoA = read_tsv("DcoB_vs_DcoA.delta-filter.options_rql0000.coords")
```

    ## Rows: 2461 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): rid, qid, strand
    ## dbl (5): rs, re, qs, qe, error
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
chromcompare_06_DcoBvsDcoA_plot <- ggplot(chromcompare_06_DcoBvsDcoA, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position=c(.99,.01), legend.justification=c(1,0)) +
          theme(legend.key.height= unit(3, 'mm'), legend.title = element_text(size=8), legend.text = element_text(size=8)) +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm"))
```

    ## Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
    ## 3.5.0.
    ## ℹ Please use the `legend.position.inside` argument of `theme()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
#############################################################################################################
# placeholder plot
#############################################################################################################

chromcompare_dummy <- ggplot(chromcompare_06_DcoBvsDcoA, aes(x=rs, xend=re, y=qs, yend=qe, colour=strand)) +
          scale_x_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          scale_y_continuous(limits = c(0,90000000), labels = scales::comma, breaks = seq(0, 90000000, by = 20000000))+
          geom_segment(lineend = "round", linejoin = "round") +
          #geom_point(alpha=0.5) +
          theme_bw() +
          #theme(strip.text.y=element_text(angle=180, size=5),
          #      legend.position=c(.99,.01), legend.justification=c(1,0),
          #      strip.background=element_blank()) +
          theme(legend.position = "none") +
          scale_color_manual(values = c("red", "blue")) +
          theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          coord_fixed() +
          theme(plot.margin = unit(c(margins, margins, margins, margins), "mm")) +
          geom_label(x=45000000, y=45000000, label="PLACEHOLDER")


#############################################################################################################
# construct matrix of comparisons
#############################################################################################################

library(patchwork)

Diploscapter_chromcompare_matrix <- chromcompare_01_DpaBvsDpaA_plot + chromcompare_02_DcoAvsDpaA_plot + chromcompare_03_DcoBvsDpaA_plot +
                                    chromcompare_dummy + chromcompare_04_DcoAvsDpaB_plot + chromcompare_05_DcoBvsDpaB_plot +
                                    chromcompare_dummy + chromcompare_dummy + chromcompare_06_DcoBvsDcoA_plot +
                                    plot_layout(ncol = 3, nrow = 3)

Diploscapter_chromcompare_matrix
```

![](Diploscapter_collinearity_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave('chromosome_alignment_9cm_by_9cm.pdf', Diploscapter_chromcompare_matrix, width = 9, height = 9, units = 'cm')
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
