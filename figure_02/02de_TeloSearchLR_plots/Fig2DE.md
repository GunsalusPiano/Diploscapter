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
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(stringr)   

axis_text_size = 8
legned_box_size = 0.2
y_max = 1.6e-5

#turn a value from 1,000,000 to 1
million_times <- function(x){
  x*1000000
}

df <- read_tsv("Dpa.TAAGGGTAAGGC.occupancy.n1000.head.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TAAGGGTAAGGC_head_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            xlab("bp from 5' end") +
                            ylab("fraction occupancy") +
                            scale_color_manual(values = c("red", "blue")) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500)) +
                            scale_y_continuous(limits = c(0, y_max), labels = million_times) +
                            theme(legend.position = "none") +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
                            theme(axis.text.x = element_text(size = axis_text_size),axis.text.y = element_text(size = axis_text_size))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
df <- read_tsv("Dpa.TAAGGGTAAGGC.occupancy.n1000.tail.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TAAGGGTAAGGC_tail_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500), transform = "reverse") +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            xlab("bp from 3' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
                            


df <- read_tsv("Dpa.TAAGGGTAAGGC.occupancy.n1000.all.bar.txt", comment = "#")
```

    ## Rows: 6 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, motif
    ## dbl (1): occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TAAGGGTAAGGC_all_plot <- ggplot(df, aes(fill=motif, y=occupancy, x=region)) +
                            geom_bar(position="dodge", stat="identity") +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_fill_manual(values=c('red', 'blue')) +
                            scale_x_discrete(limits = c("first 1000bp",
                                                        "middle",
                                                        "last 1000bp"),
                                             labels = c("first 1kb",
                                                        "mid",
                                                        "last 1kb")) +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) +
                            #theme(legend.position = c(0.8, 0.8)) +
                            theme(legend.key.size = unit(legned_box_size, 'cm')) +
                            theme(legend.text=element_text(size=axis_text_size)) +
                            theme(legend.title=element_blank())
                            



Dpa_TAAGGGTAAGGC_head_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
Dpa_TAAGGGTAAGGC_tail_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
Dpa_TAAGGGTAAGGC_all_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
y_max = 1.6e-5

df <- read_tsv("Dpa.TAAGGG.occupancy.n1000.head.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TAAGGG_head_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500)) +
                            scale_y_continuous(limits = c(0, y_max), labels = million_times) +
                            xlab("bp from 5' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
                            theme(axis.text.x = element_text(size = axis_text_size),axis.text.y = element_text(size = axis_text_size))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
df <- read_tsv("Dpa.TAAGGG.occupancy.n1000.tail.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TAAGGG_tail_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500), transform = "reverse") +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            xlab("bp from 3' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
                            


df <- read_tsv("Dpa.TAAGGG.occupancy.n1000.all.bar.txt", comment = "#")
```

    ## Rows: 6 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, motif
    ## dbl (1): occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TAAGGG_all_plot <- ggplot(df, aes(fill=motif, y=occupancy, x=region)) +
                            geom_bar(position="dodge", stat="identity") +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_fill_manual(values=c('red', 'blue')) +
                            scale_x_discrete(limits = c("first 1000bp",
                                                        "middle",
                                                        "last 1000bp"),
                                             labels = c("first 1kb",
                                                        "mid",
                                                        "last 1kb")) +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) +
                            #theme(legend.position = c(0.8, 0.8)) +
                            theme(legend.key.size = unit(legned_box_size, 'cm')) +
                            theme(legend.text=element_text(size=axis_text_size)) +
                            theme(legend.title=element_blank())
                            



Dpa_TAAGGG_head_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
Dpa_TAAGGG_tail_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
Dpa_TAAGGG_all_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
y_max = 1.6e-5
df <- read_tsv("Dpa.TTAGGC.occupancy.n1000.head.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TTAGGC_head_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500)) +
                            scale_y_continuous(limits = c(0, y_max), labels = million_times) +
                            xlab("bp from 5' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
                            theme(axis.text.x = element_text(size = axis_text_size),axis.text.y = element_text(size = axis_text_size))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
df <- read_tsv("Dpa.TTAGGC.occupancy.n1000.tail.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TTAGGC_tail_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500), transform = "reverse") +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            xlab("bp from 3' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
                            


df <- read_tsv("Dpa.TTAGGC.occupancy.n1000.all.bar.txt", comment = "#")
```

    ## Rows: 6 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, motif
    ## dbl (1): occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dpa_TTAGGC_all_plot <- ggplot(df, aes(fill=motif, y=occupancy, x=region)) +
                            geom_bar(position="dodge", stat="identity") +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_fill_manual(values=c('red', 'blue')) +
                            scale_x_discrete(limits = c("first 1000bp",
                                                        "middle",
                                                        "last 1000bp"),
                                             labels = c("first 1kb",
                                                        "mid",
                                                        "last 1kb")) +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) +
                            #theme(legend.position = c(0.8, 0.8)) +
                            theme(legend.key.size = unit(legned_box_size, 'cm')) +
                            theme(legend.text=element_text(size=axis_text_size)) +
                            theme(legend.title=element_blank())
                            



Dpa_TTAGGC_head_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
Dpa_TTAGGC_tail_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
Dpa_TTAGGC_all_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
y_max = 3e-5
df <- read_tsv("Dco.TAAGGGTAAGGC.occupancy.n1000.head.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TAAGGGTAAGGC_head_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500)) +
                            scale_y_continuous(limits = c(0, y_max), labels = million_times) +
                            xlab("bp from 5' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
                            theme(axis.text.x = element_text(size = axis_text_size),axis.text.y = element_text(size = axis_text_size))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
df <- read_tsv("Dco.TAAGGGTAAGGC.occupancy.n1000.tail.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TAAGGGTAAGGC_tail_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500), transform = "reverse") +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            xlab("bp from 3' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
                            


df <- read_tsv("Dco.TAAGGGTAAGGC.occupancy.n1000.all.bar.txt", comment = "#")
```

    ## Rows: 6 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, motif
    ## dbl (1): occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TAAGGGTAAGGC_all_plot <- ggplot(df, aes(fill=motif, y=occupancy, x=region)) +
                            geom_bar(position="dodge", stat="identity") +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_fill_manual(values=c('red', 'blue')) +
                            scale_x_discrete(limits = c("first 1000bp",
                                                        "middle",
                                                        "last 1000bp"),
                                             labels = c("first 1kb",
                                                        "mid",
                                                        "last 1kb")) +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) +
                            #theme(legend.position = c(0.8, 0.8)) +
                            theme(legend.key.size = unit(legned_box_size, 'cm')) +
                            theme(legend.text=element_text(size=axis_text_size)) +
                            theme(legend.title=element_blank())
                            



Dco_TAAGGGTAAGGC_head_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
Dco_TAAGGGTAAGGC_tail_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
Dco_TAAGGGTAAGGC_all_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
y_max = 3e-5
df <- read_tsv("Dco.TAAGGG.occupancy.n1000.head.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TAAGGG_head_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500)) +
                            scale_y_continuous(limits = c(0, y_max), labels = million_times) +
                            xlab("bp from 5' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
                            theme(axis.text.x = element_text(size = axis_text_size),axis.text.y = element_text(size = axis_text_size))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
df <- read_tsv("Dco.TAAGGG.occupancy.n1000.tail.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TAAGGG_tail_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500), transform = "reverse") +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            xlab("bp from 3' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
                            


df <- read_tsv("Dco.TAAGGG.occupancy.n1000.all.bar.txt", comment = "#")
```

    ## Rows: 6 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, motif
    ## dbl (1): occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TAAGGG_all_plot <- ggplot(df, aes(fill=motif, y=occupancy, x=region)) +
                            geom_bar(position="dodge", stat="identity") +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_fill_manual(values=c('red', 'blue')) +
                            scale_x_discrete(limits = c("first 1000bp",
                                                        "middle",
                                                        "last 1000bp"),
                                             labels = c("first 1kb",
                                                        "mid",
                                                        "last 1kb")) +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) +
                            #theme(legend.position = c(0.8, 0.8)) +
                            theme(legend.key.size = unit(legned_box_size, 'cm')) +
                            theme(legend.text=element_text(size=axis_text_size)) +
                            theme(legend.title=element_blank())
                            



Dco_TAAGGG_head_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
Dco_TAAGGG_tail_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
Dco_TAAGGG_all_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
y_max = 3e-5
df <- read_tsv("Dco.TTAGGC.occupancy.n1000.head.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TTAGGC_head_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500)) +
                            scale_y_continuous(limits = c(0, y_max), labels = million_times) +
                            xlab("bp from 5' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
                            theme(axis.text.x = element_text(size = axis_text_size),axis.text.y = element_text(size = axis_text_size))
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
df <- read_tsv("Dco.TTAGGC.occupancy.n1000.tail.scatter.txt", comment = "#")
```

    ## Rows: 2000 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): motif
    ## dbl (3): coordinate, count, occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TTAGGC_tail_plot <- ggplot(df, aes(x=coordinate, y=occupancy)) +
                            geom_line(aes(color = motif)) +
                            ylim(c(0, y_max)) +
                            scale_x_continuous(breaks = seq(0, 1000, by=500), transform = "reverse") +
                            theme_bw() +
                            scale_color_manual(values = c("red", "blue")) +
                            xlab("bp from 3' end") +
                            ylab("fraction occupancy") +
                            theme(legend.position = "none") +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank())
                            


df <- read_tsv("Dco.TTAGGC.occupancy.n1000.all.bar.txt", comment = "#")
```

    ## Rows: 6 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (2): region, motif
    ## dbl (1): occupancy
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
Dco_TTAGGC_all_plot <- ggplot(df, aes(fill=motif, y=occupancy, x=region)) +
                            geom_bar(position="dodge", stat="identity") +
                            ylim(c(0, y_max)) +
                            theme_bw() +
                            scale_fill_manual(values=c('red', 'blue')) +
                            scale_x_discrete(limits = c("first 1000bp",
                                                        "middle",
                                                        "last 1000bp"),
                                             labels = c("first 1kb",
                                                        "mid",
                                                        "last 1kb")) +
                            theme(axis.text.x = element_text(size = axis_text_size), axis.text.y = element_text(size = axis_text_size)) +
                            theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank()) +
                            #theme(legend.position = c(0.8, 0.8)) +
                            theme(legend.key.size = unit(legned_box_size, 'cm')) +
                            theme(legend.text=element_text(size=axis_text_size)) +
                            theme(legend.title=element_blank())
                            



Dco_TTAGGC_head_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
Dco_TTAGGC_tail_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
Dco_TTAGGC_all_plot
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
library(patchwork)
```

    ## 
    ## Attaching package: 'patchwork'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     align_plots

``` r
telomere_occupancy_grid <- Dco_TAAGGGTAAGGC_head_plot + Dco_TAAGGGTAAGGC_tail_plot + Dco_TAAGGGTAAGGC_all_plot +
                            Dco_TAAGGG_head_plot + Dco_TAAGGG_tail_plot + Dco_TAAGGG_all_plot +
                            Dco_TTAGGC_head_plot + Dco_TTAGGC_tail_plot + Dco_TTAGGC_all_plot +
                            Dpa_TAAGGGTAAGGC_head_plot + Dpa_TAAGGGTAAGGC_tail_plot + Dpa_TAAGGGTAAGGC_all_plot +
                            Dpa_TAAGGG_head_plot + Dpa_TAAGGG_tail_plot + Dpa_TAAGGG_all_plot +
                            Dpa_TTAGGC_head_plot + Dpa_TTAGGC_tail_plot + Dpa_TTAGGC_all_plot +
                            plot_layout(ncol = 3, nrow = 6)


telomere_occupancy_grid
```

![](Fig2DE_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave('diploscapter_telomeric_motifs_flatter.pdf', telomere_occupancy_grid, width = 12, height = 11, units = 'cm')
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
