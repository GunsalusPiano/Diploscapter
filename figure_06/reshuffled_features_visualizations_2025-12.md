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
label_font_size <- 8

################################################################################
# These are the locations of centre-like regions                               #
################################################################################

cen1s <- 10000000
cen1e <- 13000000
cen2s <- 25750000
cen2e <- 28750000
cen3s <- 32300000
cen3e <- 35300000
cen4s <- 37540000
cen4e <- 40540000
cen5s <- 45750000
cen5e <- 48750000
cen6s <- 59500000
cen6e <- 62500000
cen7s <- 72000000
cen7e <- 75000000

################################################################################
# DpaA Nigon painting
################################################################################

tsv_file <- "Dpac.scaffold_1.genomic_busco_equiv_of_Dpac.scaffold_1.full_table.tsv"
nigon_def <- "gene2Nigon_busco20200927.tsv.gz"
windwSize <- 500000
minimumGenesPerSequence <- 15
spName <- "DpaA"

height <- 4
width <- 10

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
plNigon <- group_by(consUsco, Sequence) %>%
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
  annotate("rect", xmin=cen1s, xmax=cen1e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  annotate("rect", xmin=cen2s, xmax=cen2e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  #annotate("rect", xmin=cen3s, xmax=cen3e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  annotate("rect", xmin=cen4s, xmax=cen4e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  annotate("rect", xmin=cen5s, xmax=cen5e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  annotate("rect", xmin=cen6s, xmax=cen6e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  annotate("rect", xmin=cen7s, xmax=cen7e, ymin=0, ymax=50, alpha=0.1, fill="red") +
  geom_bar(position="stack", stat="identity", just=0) +
  ggtitle(spName) + 
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(labels = scales::comma, position="top") +
  coord_cartesian(ylim = c(0, 40), xlim = c(0, 110000000)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme_void() +
  theme(legend.position="none")

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
contacts  <- ggplot(PC_df) +
          annotate("rect", xmin=cen1s, xmax=cen1e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          annotate("rect", xmin=cen2s, xmax=cen2e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          #annotate("rect", xmin=DpaA_cen3s, xmax=DpaA_cen3e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          annotate("rect", xmin=cen4s, xmax=cen4e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          annotate("rect", xmin=cen5s, xmax=cen5e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          annotate("rect", xmin=cen6s, xmax=cen6e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          annotate("rect", xmin=cen7s, xmax=cen7e, ymin=-0.2, ymax=0.2, alpha=0.1, fill="red") +
          geom_point(mapping = aes(x=scaffold_position, y=PC3), color = "red", size = 0.1, alpha = 0.3) +
          coord_cartesian(ylim = c(-0.1, 0.1), xlim = c(0, 110000000)) +
          theme_void()



################################################################################
# DpaA re-shuffled Nigon painting
################################################################################

tsv_file <- "Dpac.scaffold_1.reshuffled.Nigon.tsv"
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
plNigon_reshuffled <- group_by(consUsco, Sequence) %>%
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
  geom_bar(position="stack", stat="identity", just=0) +
  ggtitle(spName) +
  theme_bw() +
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(labels = scales::comma, position="top") +
  coord_cartesian(ylim = c(0, 40), xlim = c(0, 110000000)) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme_void() +
  theme(legend.position="none")


################################################################################
# DpaA PoreC PC3 reshuffled
################################################################################

PC_file = "Dpac.scaffold_1.reshuffled.PC3.window100000.increment100000.txt"

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
contacts_reshuffled  <- ggplot(PC_df) +
          geom_point(mapping = aes(x=scaffold_position, y=PC3), color = "red", size = 0.1, alpha = 0.3) +
          theme_bw() +
          scale_x_continuous(labels = scales::comma) +
          coord_cartesian(ylim = c(-0.1, 0.1), xlim = c(0, 110000000)) +
          theme(text = element_text(family="Helvetica")) +
          theme(axis.title.y=element_text(angle=0),) +
          #stat_smooth(aes(x=scaffold_position, y=PC3), method="loess",n=1000,fullrange=TRUE, span=0.1, color = "red", linewidth = 0.5) +
          #theme(axis.title.x=element_blank(), axis.text.x=element_blank())
          theme_void()



################################################################################
# Construct DpaA plots grid
################################################################################

DpaA_grid <- plot_grid(plNigon,
                       contacts,
                       plNigon_reshuffled,
                       contacts_reshuffled,
                       ncol = 1, nrow = 4, align = "v", axis = "lr", rel_heights = c(2,1,2,1))

DpaA_grid
```

![](reshuffled_features_visualizations_2025-12_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave('DpaA.pdf', DpaA_grid, width = 10, height = 6, units = 'cm')
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
