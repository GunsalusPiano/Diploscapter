R Notebook
================

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you
execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk
or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*.

``` r
# This code attempts to visualize P. pacificus chromosomes (I, II etc.) data with respect to
#    - Nigon painting
#    - GC content
#    - repeat content (determined by RepeatMasker)
#    - exonic DNA content

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
library(gtools)
library(ggplot2)
library(gggenes)
library(ggh4x)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
label_font_size <- 8
#turn a value from 1,000,000 to 1
millionth<-function(x){
  x/1000000
}



ChrIV_end <- 35816911
ChrX_end <- 31480924
ChrII_end <- 25511976
ChrV_end <- 24657445
ChrIII_end <- 24638849
ChrI_end <- 24218945

cenIs <- 11000000
cenIIs <- 12000000
cenIIIs <- 10000000
cenIVs <- 21000000
cenVs <- 10000000
cenXs1 <- 9000000
cenXs2 <- 23000000

cenIe <- 14000000
cenIIe <- 15000000
cenIIIe <- 13000000
cenIVe <- 24000000
cenVe <- 13000000
cenXe1 <- 12000000
cenXe2 <- 26000000


GC_yaxis_min = 0.35
GC_yaxis_max = 0.50

repeats_yaxis_min = 0
repeats_yaxis_max = 0.5

exonic_yaxis_min = 0
exonic_yaxis_max = 0.6

HFC_fill_alpha <- 0.2


millionths<-function(x){
  x/1000000
}

current_ymin <- 0
current_ymax <- 40

tsv_file <- "Pexsp_genomic_busco_equiv_of_full_table_from_BUSCO.tsv"
nigon_def <- "gene2Nigon_busco20200927.tsv.gz"
windwSize <- 500000
minimumGenesPerSequence <- 15
spName <- "Pexsp"

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
  facet_grid(. ~ scaffold_f, switch = "y") +
  geom_rect(data = data.frame(scaffold_f = "ChrI"), aes(xmin = cenIs, xmax = cenIe, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(scaffold_f = "ChrII"), aes(xmin = cenIIs, xmax = cenIIe, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(scaffold_f = "ChrIII"), aes(xmin = cenIIIs, xmax = cenIIIe, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(scaffold_f = "ChrIV"), aes(xmin = cenIVs, xmax = cenIVe, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(scaffold_f = "ChrV"), aes(xmin = cenVs, xmax = cenVe, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(scaffold_f = "ChrX"), aes(xmin = cenXs1, xmax = cenXe1, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_rect(data = data.frame(scaffold_f = "ChrX"), aes(xmin = cenXs2, xmax = cenXe2, ymin = -5, ymax = 45), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
  geom_bar(position="stack", stat="identity") +
  ggtitle(spName) +
  theme_bw() +coord_cartesian(ylim=c(current_ymin, current_ymax)) +  
  #scale_y_continuous(breaks = scales::pretty_breaks(4), position = "left") +
  scale_x_continuous(position="bottom") +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(ncol = 1, title = "Nigon")) +
  theme(text = element_text(family="Helvetica"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = ggtext::element_markdown(),
        axis.text=element_text(size=label_font_size))
        #legend.key.size = unit(0.3, 'cm')) +

plNigon
```

![](supp_figure_2024-12-28_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# fraction GC plot for CeleI
current_ymin <- 0.36
current_ymax <- 0.48


fraction_GC_file = "Pexsp.fraction_GC.window100000.increment100000.txt"

fraction_GC_df <- read_tsv(fraction_GC_file)
```

    ## Rows: 1667 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): scaffold
    ## dbl (2): coordinate, fraction_GC
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# plot
fraction_GC  <- ggplot(fraction_GC_df) +
          geom_point(mapping = aes(x=coordinate, y=fraction_GC), color = "blue", size = 0.2, alpha = 0.3) +
          facet_grid(. ~ scaffold) +
          theme_bw() +
          geom_rect(data = data.frame(scaffold= "ChrI"), aes(xmin = cenIs, xmax = cenIe, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrII"), aes(xmin = cenIIs, xmax = cenIIe, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIII"), aes(xmin = cenIIIs, xmax = cenIIIe, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIV"), aes(xmin = cenIVs, xmax = cenIVe, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrV"), aes(xmin = cenVs, xmax = cenVe, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrX"), aes(xmin = cenXs1, xmax = cenXe1, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold="ChrX"), aes(xmin = cenXs2, xmax = cenXe2, ymin = 0.34, ymax = 0.5), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          coord_cartesian(ylim=c(current_ymin, current_ymax)) +        
          stat_smooth(aes(x=coordinate, y=fraction_GC), method="loess", n=1000, fullrange=TRUE, span=0.35, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                strip.text=element_blank(),
                axis.text=element_text(size=label_font_size))




fraction_GC
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

![](supp_figure_2024-12-28_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# fraction RepeatScout repeats plot for CeleI

current_ymin <- 0.05
current_ymax <- 0.6

fraction_RepeatScout_rep_file = "Pexsp.08.repeats.window100000.increment100000.txt"

fraction_RepeatScout_rep_df <- read_tsv(fraction_RepeatScout_rep_file)
```

    ## Rows: 1667 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): scaffold
    ## dbl (4): coordinate, Simple_repeat, Unspecified, sum
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
fraction_RepeatScoutRep <- ggplot(fraction_RepeatScout_rep_df) +
          geom_point(mapping = aes(x=coordinate, y=sum), color = "darkgreen", size = 0.2, alpha = 0.3) +
          facet_grid(. ~ scaffold) +
          theme_bw() +
          geom_rect(data = data.frame(scaffold= "ChrI"), aes(xmin = cenIs, xmax = cenIe, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrII"), aes(xmin = cenIIs, xmax = cenIIe, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIII"), aes(xmin = cenIIIs, xmax = cenIIIe, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIV"), aes(xmin = cenIVs, xmax = cenIVe, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrV"), aes(xmin = cenVs, xmax = cenVe, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrX"), aes(xmin = cenXs1, xmax = cenXe1, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold="ChrX"), aes(xmin = cenXs2, xmax = cenXe2, ymin = 0.03, ymax = 0.62), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          coord_cartesian(ylim=c(current_ymin, current_ymax)) +        
          stat_smooth(aes(x=coordinate, y=sum), method="loess",n=1000,fullrange=TRUE, span=0.35, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                strip.text=element_blank(),
                axis.text=element_text(size=label_font_size))

fraction_RepeatScoutRep
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

![](supp_figure_2024-12-28_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Fraction exonic DNA for CeleI

current_ymin <- 0
current_ymax <- 0.5

exonic_DNA_content_file <- "Pexsp.fraction_exonic_DNA.window100000.increment100000.txt"

exonic_DNA_content_df <- read_tsv(exonic_DNA_content_file)
```

    ## Rows: 1667 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): scaffold
    ## dbl (2): coordinate, fraction_exonic_DNA
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
fraction_exonicDNA  <- ggplot(exonic_DNA_content_df) +
          geom_point(mapping = aes(x=coordinate, y=fraction_exonic_DNA), color = "darkorange2", size = 0.2, alpha = 0.3) +
          facet_grid(. ~ scaffold) +
          theme_bw() +
          geom_rect(data = data.frame(scaffold= "ChrI"), aes(xmin = cenIs, xmax = cenIe, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrII"), aes(xmin = cenIIs, xmax = cenIIe, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIII"), aes(xmin = cenIIIs, xmax = cenIIIe, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIV"), aes(xmin = cenIVs, xmax = cenIVe, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrV"), aes(xmin = cenVs, xmax = cenVe, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrX"), aes(xmin = cenXs1, xmax = cenXe1, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold="ChrX"), aes(xmin = cenXs2, xmax = cenXe2, ymin = -0.02, ymax = 0.52), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          coord_cartesian(ylim=c(current_ymin, current_ymax)) +
          stat_smooth(aes(x=coordinate, y=fraction_exonic_DNA), method="loess",n=1000,fullrange=TRUE, span=0.2, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                axis.text.x=element_blank(),
                strip.text=element_blank(),
                axis.text=element_text(size=label_font_size))

fraction_exonicDNA
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

![](supp_figure_2024-12-28_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
################################################################################
# Pexsp HiC PC1 
################################################################################

current_ymin = -0.01
current_ymax = 0.08

PC_file = "Pexsp.PC1.window100000.increment100000.txt"

PC_df <- read_tsv(PC_file)
```

    ## Rows: 1667 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (1): scaffold
    ## dbl (2): coordinate, PC1
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
HiC_PCA  <- ggplot(PC_df) +
          geom_point(mapping = aes(x=coordinate, y=PC1), color = "red", size = 0.2, alpha = 0.3) +
          facet_grid(. ~ scaffold) +
          theme_bw() +
          geom_rect(data = data.frame(scaffold= "ChrI"), aes(xmin = cenIs, xmax = cenIe, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrII"), aes(xmin = cenIIs, xmax = cenIIe, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIII"), aes(xmin = cenIIIs, xmax = cenIIIe, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrIV"), aes(xmin = cenIVs, xmax = cenIVe, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrV"), aes(xmin = cenVs, xmax = cenVe, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold= "ChrX"), aes(xmin = cenXs1, xmax = cenXe1, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          geom_rect(data = data.frame(scaffold="ChrX"), aes(xmin = cenXs2, xmax = cenXe2, ymin = -0.03, ymax = 0.1), alpha = HFC_fill_alpha, fill="red", inherit.aes = FALSE) +
          coord_cartesian(ylim=c(current_ymin, current_ymax)) +
          scale_x_continuous(labels=millionth) +
          stat_smooth(aes(x=coordinate, y=PC1), method="loess",n=1000,fullrange=TRUE, span=0.1, linewidth=0.2, color="black") +
          theme(text=element_text(family="Helvetica"),
                axis.title=element_blank(),
                strip.text=element_blank(),
                axis.text=element_text(size=label_font_size))

HiC_PCA
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

![](supp_figure_2024-12-28_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
library(ggpubr)
```

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

``` r
Pexsp_grid <- plot_grid(plNigon,
                      fraction_RepeatScoutRep,
                      fraction_exonicDNA,
                      fraction_GC,
                      HiC_PCA,
                      ncol = 1, nrow = 5, align = "v", axis = "lr", rel_heights = c(2,1,1,1,1.2))
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1362 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

``` r
Pexsp_grid
```

![](supp_figure_2024-12-28_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("Pexsp_traits.pdf", Pexsp_grid, width = 15, height = 12, units = "cm")
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
