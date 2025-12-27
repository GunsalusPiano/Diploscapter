R Notebook
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

``` r
#turn a value from 1,000,000 to 1
millionth<-function(x){
  x/1000000
}

# turn a value from 1 to 1000
thousand<-function(x){
  x*1000
}

label_font_size = 8
```

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
library(pracma)
```

    ## 
    ## Attaching package: 'pracma'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     cross

``` r
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(gtools)
```

    ## 
    ## Attaching package: 'gtools'
    ## 
    ## The following object is masked from 'package:pracma':
    ## 
    ##     logit

``` r
library(ggpubr)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     get_legend
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
HFC_fill_alpha <- 0.2

################################################################################
# These are the locations of centre-like regions                               #
################################################################################

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

################################################################################
# Dpa scaffold_1 Nigon painting
################################################################################

tsv_file <- "Dpac.scaffold_1.genomic_busco_equiv_of_Dpac.scaffold_1.full_table.tsv"
nigon_def <- "gene2Nigon_busco20200927.tsv.gz"
windwSize <- 500000
minimumGenesPerSequence <- 15
spName <- "DpaA"

height <- 4
width <- 10

# Load data
nigonDict <- read_tsv(nigon_def,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(tsv_file,
                  col_names = c("Busco_id", "Status", "Sequence",
                                "start", "end", "strand", "Score", "Length",
                                "OrthoDB_url", "Description"),
                  col_types = c("ccciicdicc"),
                  comment = "#"))

# Specify Nigon colors
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Filter data
fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")

consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence, mxGpos > windwSize * 2)



# Plot
DpaA_plNigon <- group_by(consUsco, Sequence) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), windwSize),
                                            labels = seq(windwSize, max(stPos), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  count(ints, nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
  ggplot(aes(fill=nigon, y=n, x=ints-windwSize)) + 
  facet_grid(scaffold_f ~ ., switch = "y") +
  annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaA_cen2s, xmax=DpaA_cen2e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaA_cen4s, xmax=DpaA_cen4e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaA_cen5s, xmax=DpaA_cen5e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaA_cen6s, xmax=DpaA_cen6e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaA_cen7s, xmax=DpaA_cen7e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  geom_bar(position="stack", stat="identity") +
  ggtitle(spName, ) +
  theme_bw() +
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(labels=millionth, position="top") +
  coord_cartesian(ylim = c(0, 40), xlim = c(0, 90000000)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #panel.border = element_blank(),
        #strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(family="Helvetica"),
        plot.title = ggtext::element_markdown(),
        legend.position="none",
        axis.text.y=element_text(size=label_font_size),
        axis.text.x=element_text(size=label_font_size))


################################################################################
# DpaA %GC
################################################################################

fraction_GC_file = "Dpac.scaffold_1.fraction_GC.window100000.increment100000.txt"

fraction_GC_df <- read_tsv(fraction_GC_file)
```

    ## Rows: 852 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_GC
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# plot
# the dark green is the rDNA cluster
DpaA_fraction_GC  <- ggplot(fraction_GC_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen2s, xmax=DpaA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen4s, xmax=DpaA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen5s, xmax=DpaA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen6s, xmax=DpaA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen7s, xmax=DpaA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_GC), color = "blue", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(0.33, 0.47), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 1, alpha = .2, fill = "darkgreen") +
          stat_smooth(aes(x=scaffold_position, y=fraction_GC), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size),
                axis.text.x=element_text(size=label_font_size))


################################################################################
# DpaA %repeat
################################################################################

fraction_RepeatScout_rep_file = "Dpac.08.scaffold_1.repeats.window100000.increment100000.txt"

fraction_RepeatScout_rep_df <- read_tsv(fraction_RepeatScout_rep_file)
```

    ## Rows: 852 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (4): scaffold_position, Simple_repeat, Unspecified, sum
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaA_fraction_repeats  <- ggplot(fraction_RepeatScout_rep_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen2s, xmax=DpaA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen4s, xmax=DpaA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen5s, xmax=DpaA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen6s, xmax=DpaA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen7s, xmax=DpaA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=sum), color = "darkgreen", size = 0.2, alpha = 0.3) +
          theme_bw() +
          #scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = -1, ymax = 1, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(fraction_repeats_all_annotated, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=sum), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))
          

################################################################################
# DpaA %exonic
################################################################################

exonic_DNA_content_file <- "Dpac.scaffold_1.fraction_exonic.window100000.increment100000.txt"

exonic_DNA_content_df <- read_tsv(exonic_DNA_content_file)
```

    ## Rows: 852 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_exonic
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaA_fraction_exonic <- ggplot(exonic_DNA_content_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen2s, xmax=DpaA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen4s, xmax=DpaA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen5s, xmax=DpaA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen6s, xmax=DpaA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen7s, xmax=DpaA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_exonic), color = "darkorange2", size = 0.2, alpha = 0.3) +
          theme_bw() +
          #scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=fraction_exonic), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# Dpa %NlaIII
################################################################################

RE_coverage_file <- "Dpac.scaffold_1.fraction_NlaIII.window100000.increment100000.txt"

RE_coverage_df <- read_tsv(RE_coverage_file)
```

    ## Rows: 852 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_NlaIII
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaA_fraction_RE <- ggplot(RE_coverage_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen2s, xmax=DpaA_cen2e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen4s, xmax=DpaA_cen4e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen5s, xmax=DpaA_cen5e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen6s, xmax=DpaA_cen6e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen7s, xmax=DpaA_cen7e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_NlaIII), color = "darkorchid4", size = 0.2, alpha = 0.3) +
          theme_bw() +
          #scale_x_continuous(labels = scales::comma) +
          scale_y_continuous(labels = thousand) +
          coord_cartesian(ylim = c(0.007, 0.016), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 0.025, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(NlaIII_site_cov, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=fraction_NlaIII), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))
          

################################################################################
# DpaA PoreC PC3
################################################################################

PC_file = "Dpac.scaffold_1.PC3.window100000.increment100000.txt"

PC_df <- read_tsv(PC_file)
```

    ## Rows: 852 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, PC3
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaA_PoreC_PC3  <- ggplot(PC_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen2s, xmax=DpaA_cen2e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen4s, xmax=DpaA_cen4e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen5s, xmax=DpaA_cen5e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen6s, xmax=DpaA_cen6e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaA_cen7s, xmax=DpaA_cen7e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=PC3), color = "red", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(-0.1, 0.1), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=PC3), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))


################################################################################
# Construct DpaA plots grid
################################################################################

DpaA_grid <- plot_grid(DpaA_plNigon,
                      DpaA_fraction_GC,
                      DpaA_fraction_repeats,
                      DpaA_fraction_exonic,
                      DpaA_fraction_RE,
                      DpaA_PoreC_PC3,
                      ncol = 1, nrow = 6, align = "v", axis = "lr", rel_heights = c(5,2,2,2,2,2))
```

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
DpaA_grid
```

![](genomic_features_visualizations_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
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

################################################################################
# Dpa scaffold_2 Nigon painting
################################################################################

tsv_file <- "Dpac.scaffold_2.genomic_busco_equiv_of_Dpac.scaffold_2.full_table.tsv"
nigon_def <- "gene2Nigon_busco20200927.tsv.gz"
windwSize <- 500000
minimumGenesPerSequence <- 15
spName <- "DpaB"

height <- 4
width <- 10

# Load data
nigonDict <- read_tsv(nigon_def,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(tsv_file,
                  col_names = c("Busco_id", "Status", "Sequence",
                                "start", "end", "strand", "Score", "Length",
                                "OrthoDB_url", "Description"),
                  col_types = c("ccciicdicc"),
                  comment = "#"))

# Specify Nigon colors
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Filter data
fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")

consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence, mxGpos > windwSize * 2)



# Plot
DpaB_plNigon <- group_by(consUsco, Sequence) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), windwSize),
                                            labels = seq(windwSize, max(stPos), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  count(ints, nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
  ggplot(aes(fill=nigon, y=n, x=ints-windwSize)) + 
  facet_grid(scaffold_f ~ ., switch = "y") +
  annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaB_cen2s, xmax=DpaB_cen2e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  #annotate("rect", xmin=DpaB_cen3s, xmax=DpaB_cen3e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaB_cen4s, xmax=DpaB_cen4e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaB_cen5s, xmax=DpaB_cen5e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaB_cen6s, xmax=DpaB_cen6e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DpaB_cen7s, xmax=DpaB_cen7e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  geom_bar(position="stack", stat="identity") +
  ggtitle(spName) +
  theme_bw() +
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(labels=millionth, position="top") +
  coord_cartesian(ylim = c(0, 40), xlim = c(0, 90000000)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #panel.border = element_blank(),
        #strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(family="Helvetica"),
        plot.title = ggtext::element_markdown(),
        legend.position="none",
        axis.text.y=element_text(size=label_font_size),
        axis.text.x=element_text(size=label_font_size))


################################################################################
# DpaB %GC
################################################################################

fraction_GC_file = "Dpac.scaffold_2.fraction_GC.window100000.increment100000.txt"

fraction_GC_df <- read_tsv(fraction_GC_file)
```

    ## Rows: 888 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_GC
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# plot
# the dark green is the rDNA cluster
DpaB_fraction_GC  <- ggplot(fraction_GC_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen2s, xmax=DpaB_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaB_cen3s, xmax=DpaB_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen4s, xmax=DpaB_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen5s, xmax=DpaB_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen6s, xmax=DpaB_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen7s, xmax=DpaB_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_GC), color = "blue", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(0.33, 0.47), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 1, alpha = .2, fill = "darkgreen") +
          theme(t) +
          stat_smooth(aes(x=scaffold_position, y=fraction_GC), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size),
                axis.text.x=element_text(size=label_font_size))


################################################################################
# DpaB %repeat
################################################################################

fraction_RepeatScout_rep_file = "Dpac.08.scaffold_2.repeats.window100000.increment100000.txt"

fraction_RepeatScout_rep_df <- read_tsv(fraction_RepeatScout_rep_file)
```

    ## Rows: 888 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (4): scaffold_position, Simple_repeat, Unspecified, sum
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaB_fraction_repeats  <- ggplot(fraction_RepeatScout_rep_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen2s, xmax=DpaB_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaB_cen3s, xmax=DpaB_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen4s, xmax=DpaB_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen5s, xmax=DpaB_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen6s, xmax=DpaB_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen7s, xmax=DpaB_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=sum), color = "darkgreen", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = -1, ymax = 1, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(fraction_repeats_all_annotated, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=sum), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))


################################################################################
# DpaB %exonic
################################################################################

exonic_DNA_content_file <- "Dpac.scaffold_2.fraction_exonic.window100000.increment100000.txt"

exonic_DNA_content_df <- read_tsv(exonic_DNA_content_file)
```

    ## Rows: 888 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_exonic
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaB_fraction_exonic <- ggplot(exonic_DNA_content_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen2s, xmax=DpaB_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaB_cen3s, xmax=DpaB_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen4s, xmax=DpaB_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen5s, xmax=DpaB_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen6s, xmax=DpaB_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen7s, xmax=DpaB_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_exonic), color = "darkorange2", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=fraction_exonic), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# DpaB %NlaIII
################################################################################

RE_coverage_file <- "Dpac.scaffold_2.fraction_NlaIII.window100000.increment100000.txt"

RE_coverage_df <- read_tsv(RE_coverage_file)
```

    ## Rows: 888 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_NlaIII
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaB_fraction_RE <- ggplot(RE_coverage_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen2s, xmax=DpaB_cen2e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaB_cen3s, xmax=DpaB_cen3e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen4s, xmax=DpaB_cen4e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen5s, xmax=DpaB_cen5e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen6s, xmax=DpaB_cen6e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen7s, xmax=DpaB_cen7e, ymin=0, ymax=0.03, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_NlaIII), color = "darkorchid4", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          scale_y_continuous(labels = thousand) +
          coord_cartesian(ylim = c(0.007, 0.016), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 0.025, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(NlaIII_site_cov, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=fraction_NlaIII), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# DpaB PoreC PC4
################################################################################

PC_file = "Dpac.scaffold_2.PC4.window100000.increment100000.txt"

PC_df <- read_tsv(PC_file)
```

    ## Rows: 888 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, PC4
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DpaB_PoreC_PC4  <- ggplot(PC_df) +
          annotate("rect", xmin=DpaA_cen1s, xmax=DpaA_cen1e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen2s, xmax=DpaB_cen2e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DpaB_cen3s, xmax=DpaB_cen3e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen4s, xmax=DpaB_cen4e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen5s, xmax=DpaB_cen5e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen6s, xmax=DpaB_cen6e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DpaB_cen7s, xmax=DpaB_cen7e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=-PC4), color = "red", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(-0.1, 0.1), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=-PC4), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))


################################################################################
# Construct DpaA plots grid
################################################################################

DpaB_grid <- plot_grid(DpaB_plNigon,
                      DpaB_fraction_GC,
                      DpaB_fraction_repeats,
                      DpaB_fraction_exonic,
                      DpaB_fraction_RE,
                      DpaB_PoreC_PC4,
                      ncol = 1, nrow = 6, align = "v", axis = "lr", rel_heights = c(5,2,2,2,2,2))
```

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
DpaB_grid
```

![](genomic_features_visualizations_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
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

################################################################################
# DcoA Nigon painting
################################################################################

tsv_file <- "Dcor.SUPER_2.genomic_busco_equiv_of_Dcor.SUPER_2.full_table.tsv"
nigon_def <- "gene2Nigon_busco20200927.tsv.gz"
windwSize <- 500000
minimumGenesPerSequence <- 15
spName <- "DcoA"

height <- 4
width <- 10

# Load data
nigonDict <- read_tsv(nigon_def,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(tsv_file,
                  col_names = c("Busco_id", "Status", "Sequence",
                                "start", "end", "strand", "Score", "Length",
                                "OrthoDB_url", "Description"),
                  col_types = c("ccciicdicc"),
                  comment = "#"))

# Specify Nigon colors
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Filter data
fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")

consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence, mxGpos > windwSize * 2)



# Plot
DcoA_plNigon <- group_by(consUsco, Sequence) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), windwSize),
                                            labels = seq(windwSize, max(stPos), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  count(ints, nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
  ggplot(aes(fill=nigon, y=n, x=ints-windwSize)) + 
  facet_grid(scaffold_f ~ ., switch = "y") +
  annotate("rect", xmin=DcoA_cen1s, xmax=DcoA_cen1e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoA_cen2s, xmax=DcoA_cen2e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  #annotate("rect", xmin=DcoA_cen3s, xmax=DcoA_cen3e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoA_cen4s, xmax=DcoA_cen4e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoA_cen5s, xmax=DcoA_cen5e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoA_cen6s, xmax=DcoA_cen6e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoA_cen7s, xmax=DcoA_cen7e, ymin=0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  geom_bar(position="stack", stat="identity") +
  ggtitle(spName) +
  theme_bw() +
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(labels=millionth, position="top") +
  coord_cartesian(ylim = c(0, 40), xlim = c(0, 90000000)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #panel.border = element_blank(),
        #strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(family="Helvetica"),
        plot.title = ggtext::element_markdown(),
        legend.position="none",
        axis.text.y=element_text(size=label_font_size),
        axis.text.x=element_text(size=label_font_size))

DcoA_plNigon
```

![](genomic_features_visualizations_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
################################################################################
# DcoA %GC
################################################################################

fraction_GC_file = "Dcor.SUPER_2.fraction_GC.window100000.increment100000.txt"

fraction_GC_df <- read_tsv(fraction_GC_file)
```

    ## Rows: 820 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_GC
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# plot
# the dark green is the rDNA cluster
DcoA_fraction_GC  <- ggplot(fraction_GC_df) +
          annotate("rect", xmin=DcoA_cen1s, xmax=DcoA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen2s, xmax=DcoA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoA_cen3s, xmax=DcoA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen4s, xmax=DcoA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen5s, xmax=DcoA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen6s, xmax=DcoA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen7s, xmax=DcoA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_GC), color = "blue", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(0.33, 0.47), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 1, alpha = .2, fill = "darkgreen") +
          stat_smooth(aes(x=scaffold_position, y=fraction_GC), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size),
                axis.text.x=element_text(size=label_font_size))


################################################################################
# DcoA %repeat
################################################################################

fraction_RepeatScout_rep_file = "Dcor.08.SUPER_2.repeats.window100000.increment100000.txt"

fraction_RepeatScout_rep_df <- read_tsv(fraction_RepeatScout_rep_file)
```

    ## Rows: 820 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (4): scaffold_position, Unspecified, Simple_repeat, sum
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoA_fraction_repeats  <- ggplot(fraction_RepeatScout_rep_df) +
          annotate("rect", xmin=DcoA_cen1s, xmax=DcoA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen2s, xmax=DcoA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoA_cen3s, xmax=DcoA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen4s, xmax=DcoA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen5s, xmax=DcoA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen6s, xmax=DcoA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen7s, xmax=DcoA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=sum), color = "darkgreen", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = -1, ymax = 1, alpha = .2, fill = "darkgreen")+
          theme(text = element_text(family="Helvetica")) +
          #geom_line(aes(x=scaffold_position, y=whittaker(fraction_repeats_all_annotated, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=sum), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# DcoA %exonic
################################################################################

exonic_DNA_content_file <- "Dcor.SUPER_2.fraction_exonic.window100000.increment100000.txt"

exonic_DNA_content_df <- read_tsv(exonic_DNA_content_file)
```

    ## Rows: 820 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_exonic
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoA_fraction_exonic <- ggplot(exonic_DNA_content_df) +
          annotate("rect", xmin=DcoA_cen1s, xmax=DcoA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen2s, xmax=DcoA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoA_cen3s, xmax=DcoA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen4s, xmax=DcoA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen5s, xmax=DcoA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen6s, xmax=DcoA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen7s, xmax=DcoA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_exonic), color = "darkorange2", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=fraction_exonic), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# DcoA RE (Arima Genomics REs)
################################################################################

RE_coverage_file <- "Dcor.SUPER_2.fraction_Arima_RE.window100000.increment100000.txt"

RE_coverage_df <- read_tsv(RE_coverage_file)
```

    ## Rows: 820 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_Arima_RE
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoA_fraction_RE  <- ggplot(RE_coverage_df) +
          annotate("rect", xmin=DcoA_cen1s, xmax=DcoA_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen2s, xmax=DcoA_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoA_cen3s, xmax=DcoA_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen4s, xmax=DcoA_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen5s, xmax=DcoA_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen6s, xmax=DcoA_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen7s, xmax=DcoA_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_Arima_RE), color = "darkorchid4", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          scale_y_continuous(labels = thousand) +
          coord_cartesian(ylim = c(0.024, 0.044), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 0.025, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(NlaIII_site_cov, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=fraction_Arima_RE), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))


################################################################################
# DcoA HiC PC1 
################################################################################

PC_file = "Dcor.SUPER_2.PC1.window100000.increment100000.txt"

PC_df <- read_tsv(PC_file)
```

    ## Rows: 820 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, PC1
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoA_HiC  <- ggplot(PC_df) +
          annotate("rect", xmin=DcoA_cen1s, xmax=DcoA_cen1e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen2s, xmax=DcoA_cen2e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoA_cen3s, xmax=DcoA_cen3e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen4s, xmax=DcoA_cen4e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen5s, xmax=DcoA_cen5e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen6s, xmax=DcoA_cen6e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoA_cen7s, xmax=DcoA_cen7e, ymin=-0.2, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=PC1), color = "red", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(-0.08, 0.14), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=PC1), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))


################################################################################
# Construct DcoA plots grid
################################################################################

library(patchwork)
```

    ## 
    ## Attaching package: 'patchwork'
    ## 
    ## The following object is masked from 'package:cowplot':
    ## 
    ##     align_plots

``` r
lots_of_empty_spaces <- plot_spacer()

DcoA_grid <- plot_grid(DcoA_plNigon,
                      DcoA_fraction_GC,
                      DcoA_fraction_repeats,
                      DcoA_fraction_exonic,
                      DcoA_fraction_RE,
                      DcoA_HiC,
                      ncol = 1, nrow = 6, align = "v", axis = "lr", rel_heights = c(5,2,2,2,2,2))
```

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
DcoA_grid
```

![](genomic_features_visualizations_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
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

################################################################################
# DcoB Nigon painting
################################################################################

tsv_file <- "Dcor.SUPER_1.genomic_busco_equiv_of_Dcor.SUPER_1.full_table.tsv"
nigon_def <- "gene2Nigon_busco20200927.tsv.gz"
windwSize <- 500000
minimumGenesPerSequence <- 15
spName <- "DcoB"

height <- 4
width <- 10

# Load data
nigonDict <- read_tsv(nigon_def,
                      col_types = c(col_character(), col_character()))
busco <- suppressWarnings(read_tsv(tsv_file,
                  col_names = c("Busco_id", "Status", "Sequence",
                                "start", "end", "strand", "Score", "Length",
                                "OrthoDB_url", "Description"),
                  col_types = c("ccciicdicc"),
                  comment = "#"))

# Specify Nigon colors
cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

# Filter data
fbusco <- filter(busco, !Status %in% c("Missing")) %>%
  left_join(nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
  mutate(nigon = ifelse(is.na(nigon), "-", nigon),
         stPos = start) %>%
  filter(nigon != "-")

consUsco <- group_by(fbusco, Sequence) %>%
  mutate(nGenes = n(),
         mxGpos = max(stPos)) %>%
  ungroup() %>%
  filter(nGenes > minimumGenesPerSequence, mxGpos > windwSize * 2)



# Plot
DcoB_plNigon <- group_by(consUsco, Sequence) %>%
  mutate(ints = as.numeric(as.character(cut(stPos,
                                            breaks = seq(0, max(stPos), windwSize),
                                            labels = seq(windwSize, max(stPos), windwSize)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
  count(ints, nigon) %>%
  ungroup() %>%
  mutate(scaffold_f = factor(Sequence,
                             levels = mixedsort(unique(Sequence)))) %>%
  ggplot(aes(fill=nigon, y=n, x=ints-windwSize)) + 
  facet_grid(scaffold_f ~ ., switch = "y") +
  annotate("rect", xmin=DcoB_cen1s, xmax=DcoB_cen1e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoB_cen2s, xmax=DcoB_cen2e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  #annotate("rect", xmin=DcoB_cen3s, xmax=DcoB_cen3e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoB_cen4s, xmax=DcoB_cen4e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoB_cen5s, xmax=DcoB_cen5e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoB_cen6s, xmax=DcoB_cen6e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  annotate("rect", xmin=DcoB_cen7s, xmax=DcoB_cen7e, ymin= 0, ymax=50, alpha=HFC_fill_alpha, fill="red") +
  geom_bar(position="stack", stat="identity") +
  ggtitle(spName) +
  theme_bw() +
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(labels=millionth, position="top") +
  coord_cartesian(ylim = c(0, 40), xlim = c(0, 90000000)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        #panel.border = element_blank(),
        #strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(family="Helvetica"),
        plot.title = ggtext::element_markdown(),
        legend.position="none",
        axis.text.y=element_text(size=label_font_size),
        axis.text.x=element_text(size=label_font_size))


################################################################################
# DcoB %GC
################################################################################

fraction_GC_file = "Dcor.SUPER_1.fraction_GC.window100000.increment100000.txt"

fraction_GC_df <- read_tsv(fraction_GC_file)
```

    ## Rows: 878 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_GC
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# plot
# the dark green is the rDNA cluster
DcoB_fraction_GC  <- ggplot(fraction_GC_df) +
          annotate("rect", xmin=DcoB_cen1s, xmax=DcoB_cen1e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen2s, xmax=DcoB_cen2e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoB_cen3s, xmax=DcoB_cen3e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen4s, xmax=DcoB_cen4e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen5s, xmax=DcoB_cen5e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen6s, xmax=DcoB_cen6e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen7s, xmax=DcoB_cen7e, ymin= 0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_GC), color = "blue", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(0.33, 0.47), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 1, alpha = .2, fill = "darkgreen") +
          stat_smooth(aes(x=scaffold_position, y=fraction_GC), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                #axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size),
                axis.text.x=element_text(size=label_font_size))


################################################################################
# DcoB %repeat
################################################################################

fraction_RepeatScout_rep_file = "Dcor.08.SUPER_1.repeats.window100000.increment100000.txt"

fraction_RepeatScout_rep_df <- read_tsv(fraction_RepeatScout_rep_file)
```

    ## Rows: 878 Columns: 4
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (4): scaffold_position, Unspecified, Simple_repeat, sum
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoB_fraction_repeats  <- ggplot(fraction_RepeatScout_rep_df) +
          annotate("rect", xmin=DcoB_cen1s, xmax=DcoB_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen2s, xmax=DcoB_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoB_cen3s, xmax=DcoB_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen4s, xmax=DcoB_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen5s, xmax=DcoB_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen6s, xmax=DcoB_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen7s, xmax=DcoB_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=sum), color = "darkgreen", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = -1, ymax = 1, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(fraction_repeats_all_annotated, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=sum), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# DcoB %exonic
################################################################################

exonic_DNA_content_file <- "Dcor.SUPER_1.fraction_exonic.window100000.increment100000.txt"

exonic_DNA_content_df <- read_tsv(exonic_DNA_content_file)
```

    ## Rows: 878 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_exonic
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoB_fraction_exonic <- ggplot(exonic_DNA_content_df) +
          annotate("rect", xmin=DcoB_cen1s, xmax=DcoB_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen2s, xmax=DcoB_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoB_cen3s, xmax=DcoB_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen4s, xmax=DcoB_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen5s, xmax=DcoB_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen6s, xmax=DcoB_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen7s, xmax=DcoB_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_exonic), color = "darkorange2", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(0, 0.5), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=fraction_exonic), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# DcoA RE (Arima Genomics REs)
################################################################################

RE_coverage_file <- "Dcor.SUPER_1.fraction_Arima_RE.window100000.increment100000.txt"

RE_coverage_df <- read_tsv(RE_coverage_file)
```

    ## Rows: 878 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, fraction_Arima_RE
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoB_fraction_RE  <- ggplot(RE_coverage_df) +
          annotate("rect", xmin=DcoB_cen1s, xmax=DcoB_cen1e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen2s, xmax=DcoB_cen2e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoB_cen3s, xmax=DcoB_cen3e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen4s, xmax=DcoB_cen4e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen5s, xmax=DcoB_cen5e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen6s, xmax=DcoB_cen6e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen7s, xmax=DcoB_cen7e, ymin=0, ymax=1, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=fraction_Arima_RE), color = "darkorchid4", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          scale_y_continuous(labels = thousand) +
          coord_cartesian(ylim = c(0.024, 0.044), xlim = c(0, 90000000)) +
          #annotate("rect", xmin = 2108219, xmax = 2264010, ymin = 0, ymax = 0.025, alpha = .2, fill = "darkgreen")+
          #geom_line(aes(x=scaffold_position, y=whittaker(NlaIII_site_cov, lambda = 200)), size = 0.2) +
          stat_smooth(aes(x=scaffold_position, y=fraction_Arima_RE), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))


################################################################################
# DcoB HiC PC1 
################################################################################

PC_file = "Dcor.SUPER_1.PC1.window100000.increment100000.txt"

PC_df <- read_tsv(PC_file)
```

    ## Rows: 878 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## dbl (2): scaffold_position, PC1
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
DcoB_HiC  <- ggplot(PC_df) +
          annotate("rect", xmin=DcoB_cen1s, xmax=DcoB_cen1e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen2s, xmax=DcoB_cen2e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          #annotate("rect", xmin=DcoB_cen3s, xmax=DcoB_cen3e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen4s, xmax=DcoB_cen4e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen5s, xmax=DcoB_cen5e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen6s, xmax=DcoB_cen6e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          annotate("rect", xmin=DcoB_cen7s, xmax=DcoB_cen7e, ymin=-0.1, ymax=0.2, alpha=HFC_fill_alpha, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=PC1), color = "red", size = 0.2, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = millionth) +
          coord_cartesian(ylim = c(-0.08, 0.14), xlim = c(0, 90000000)) +
          stat_smooth(aes(x=scaffold_position, y=PC1), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_text(size=label_font_size))

################################################################################
# Construct DcoA plots grid
################################################################################

library(patchwork)

lots_of_empty_spaces <- plot_spacer()

DcoB_grid <- plot_grid(DcoB_plNigon,
                      DcoB_fraction_GC,
                      DcoB_fraction_repeats,
                      DcoB_fraction_exonic,
                      DcoB_fraction_RE,
                      DcoB_HiC,
                      ncol = 1, nrow = 6, align = "v", axis = "lr", rel_heights = c(5,2,2,2,2,2))
```

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
DcoB_grid
```

![](genomic_features_visualizations_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
all_grid <- plot_grid(DpaA_plNigon, DpaB_plNigon, DcoA_plNigon, DcoB_plNigon,
                      DpaA_PoreC_PC3, DpaB_PoreC_PC4, DcoA_HiC, DcoB_HiC,
                      DpaA_fraction_RE, DpaB_fraction_RE, DcoA_fraction_RE, DcoB_fraction_RE,
                      DpaA_fraction_repeats, DpaB_fraction_repeats, DcoA_fraction_repeats, DcoB_fraction_repeats,
                      DpaA_fraction_exonic, DpaB_fraction_exonic, DcoA_fraction_exonic, DcoB_fraction_exonic,
                      DpaA_fraction_GC, DpaB_fraction_GC, DcoA_fraction_GC, DcoB_fraction_GC,
                      ncol = 4, nrow = 6, align = "v", axis = "lr", rel_heights = c(1.75,1,1,1,1,1.2))
```

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
all_grid
```

![](genomic_features_visualizations_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave('fig04_genomic_features.v4.pdf', all_grid, width = 19.8, height = 12, units = 'cm')
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
