-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
    -   [Creating folders using
        “initProjectFolder()”](#creating-folders-using-initprojectfolder)
-   [Demographic Data](#demographic-data)
    -   [tidying the demographic data](#tidying-the-demographic-data)
-   [TACs and Blood Data](#tacs-and-blood-data)
    -   [Summary Statistics again](#summary-statistics-again)
-   [creating logan\_tstar](#creating-logan_tstar)
    -   [Plotting logan\_tstar](#plotting-logan_tstar)
    -   [All 4 tstar plots on a single
        page](#all-4-tstar-plots-on-a-single-page)
-   [Define functions for fitting the
    models](#define-functions-for-fitting-the-models)
-   [Fit kinetic models](#fit-kinetic-models)
    -   [Plot 2tcm](#plot-2tcm)
-   [new plot 2tcm](#new-plot-2tcm)
    -   [Plot ma1](#plot-ma1)
    -   [Plot Loganplot](#plot-loganplot)
-   [Test-retest](#test-retest)
    -   [trt preparation](#trt-preparation)
    -   [trt results](#trt-results)
-   [Interregional Correlation](#interregional-correlation)
-   [Vt corellation](#vt-corellation)
    -   [R-squared](#r-squared)

Aims
====

The aim of this assignment is to analyze the PK11195 data in kinfitr

Libraries
=========

CRAN libraries
--------------

installing packages

``` r
install.packages("stringr")
install.packages("corrplot")
install.packages("grid")
install.packages("gridExtra")
install.packages("RColorBrewer")
install.packages("psych")
install.packages("readxl")
install.packages("pracma")
install.packages("lme4")
install.packages("rjags")
install.packages("knitr")
install.packages("cowplot")
install.packages("corrplot")
install.packages("kableExtra")
install.packages("tidyverse")
```

First, the libraries for the analysis and plotting are loaded.

``` r
library(stringr)
library(corrplot)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(psych)
library(readxl)
library(pracma)
library(lme4)
library(rjags)
library(knitr)
library(cowplot)
library(corrplot)
library(kableExtra)
library(tidyverse)
```

Non-CRAN libraries
------------------

The libraries above can be installed from CRAN. Those which cannot are
installed as follows:

``` r
install.packages("devtools")  # If you do not already have devtools
devtools::install_github("mathesong/kinfitr")
devtools::install_github("mathesong/granviller")
devtools::install_github("mvuorre/vmisc")
devtools::install_github("mathesong/kipettools")
devtools::install_github("mathesong/relfeas")
```

### Loading Non\_CRAN libraries and setting theme

``` r
library(kinfitr)
library(vmisc)
library(kipettools)
library(granviller)
library(relfeas)
library(pander)

theme_set(theme_light())
```

Creating folders using “initProjectFolder()”
--------------------------------------------

``` r
initProjectFolder()
```

Demographic Data
================

Here, the demographic data is loaded in.

``` r
demog <- read_excel('DerivedData/TrT_chemistry_demograph.xlsx') 
```

tidying the demographic data
----------------------------

``` r
demog <- demog %>%
  select(Subjname=Akronym, Gender=Sex, Age, Genotype, 
         PET_same_day, `MBq PET1`, `MBq PET2`, bodyMass=Weight_kg) %>%  
  gather(key = PETNo, value = injRad, contains('MBq')) %>% 
  mutate(PETNo = as.numeric(str_match(PETNo, '\\d'))) %>% 
  mutate(PET = paste(Subjname, PETNo, sep='_'))
```

TACs and Blood Data
===================

First, we must read in the TAC and blood data. It should be noted that
the blood data is already dispersion corrected, and the plasma data is
already metabolite corrected, thus the plasma fraction is set to 100% in
the input data frame. It should also be noted that we set all negative
values of both blood and plasma radioactivity concentrations to zero.

``` r
tacdata <- read_csv("DerivedData/tacdata.csv") %>%
  group_by(PET) %>%
  nest(.key = 'tacdata')
blooddata <- read_csv("DerivedData/blooddata.csv") %>%
 mutate(Cbl.disp.corr = ifelse(Cbl.disp.corr < 0, 0, Cbl.disp.corr),
 Cpl..metabcorr. = ifelse(Cpl..metabcorr. < 0, 0, Cpl..metabcorr.)) %>%
  group_by(PET) %>%
  nest(.key='blooddata') %>%
  mutate(input = map(blooddata, ~blood_interp(
            t_blood = .x$ABSS.sec/60, blood=.x$Cbl.disp.corr, 
            t_plasma=.x$ABSS.sec/60, plasma=.x$Cpl..metabcorr., 
             t_parentfrac = 1, parentfrac=1 ) ))
tacs <- inner_join(tacdata, blooddata) %>%
separate(PET, c("Subjname", "PETNo"), sep='_', remove = F, convert=T)

tacs <- tacs %>%
  inner_join(demog) %>%
  arrange(PET)

saveRDS(tacs, 'DerivedData/tacs.rds')
```

\#Read tacs.rds and create delay fit and inputshift

Here, the delay and blood bolume fraction are fitted using the whole
brain ROI using 2TCM.

``` r
tacs <- readRDS('DerivedData/tacs.rds')

tacs <- tacs %>%
  ungroup() %>% 
  arrange(PET) %>% 
  group_by(PET, Subjname, PETNo) %>%
  mutate(delayFit = map2(tacdata,input, 
                         ~twotcm(t_tac = .x$Times/60, tac = .x$WB, input = .y,                                             inpshift.upper = 1, 
                                 frameStartEnd = c(1,33)))) %>%
  mutate(inpshift = map_dbl(delayFit, c("par", "inpshift"))) %>% 
  ungroup()
```

New weights with kinfitr::weights\_create

``` r
#Creating start and end times and 
#removing the old weights, so that there is no mix-up later on

# tacdata <- names(tacs$tacdata[[1]])

tacs <- tacs %>% 
  unnest(tacdata, .drop = F) %>% 
  select(-Weights)

# New weights

tacs$Weights <- weights_create(t_start = tacs$StartTime/60, t_end = (tacs$StartTime + tacs$Duration)/60, tac = tacs$WB, radioisotope = "C11")
#Nest

tacs <- tacs %>%  
  nest(-(PET:Times), .key = 'tacdata')
```

### Summary Statistics again

Below are presented some summary statistics of the demographic data.
What’s up with kurtosis?

``` r
tacs %>%
  select(Age, InjectedRadioactivity = injRad) %>%
  describe() %>%
  pandoc.table(digits=3, caption = "Summary Statistics", split.tables=Inf)
```

creating logan\_tstar
=====================

``` r
logantstar <- function(tacdata, input, inpshift) {
  Logan_tstar(t_tac = tacdata$Times/60, input = input, lowroi = tacdata$FC, medroi =   tacdata$CBL, highroi = tacdata$THA, inpshift = inpshift, gridbreaks = 3, frameStartEnd = c(1,33))
}


tacs <-  tacs %>% 
  group_by(Subjname, PETNo) %>% 
  mutate(logan_tstar = pmap(list(tacdata, input, inpshift), logantstar))
```

\#Rearrangement of the Data into Long Format

``` r
tacs_long <- tacs %>%
  unnest(tacdata, .drop = FALSE) %>%
  select(-Weights) 

tacs_long$weights = weights_create(t_start = tacs_long$StartTime/60, t_end = (tacs_long$StartTime + tacs_long$Duration)/60, tac = tacs_long$WB, radioisotope = "C11")

tacs_long <- tacs_long %>% 
gather(Region, TAC, -(PET:Times), -(StartTime:weights)) %>%
group_by(PET, Subjname, PETNo, Region) %>%
nest(-(blooddata:logan_tstar), .key = 'tacdata') 
  
  

tacs_long <- tacs %>%
  select(PET, Subjname, PETNo, blooddata:logan_tstar) %>%
  inner_join(tacs_long, by = c("PET", "Subjname", "PETNo"))
```

Plotting logan\_tstar
---------------------

``` r
set.seed(123)

tstar_fits <- tacs_long %>%  
  ungroup() %>% 
  select(logan_tstar, PET) %>% 
  select(PET, logan_tstar) %>% 
  sample_n(size =  4, replace = F) 


walk2(list(tstar_fits$logan_tstar), tstar_fits$PET, 
    ~print(plot_grid(plotlist = .x, ncol = 1, nrow = 1, labels = paste('PET:',.y), label_x = 0.5, label_y = 1, label_colour = "blue", vjust = 0.97)))
```

![](pbr28_files/figure-markdown_github/plot_tstar-1.png)![](pbr28_files/figure-markdown_github/plot_tstar-2.png)![](pbr28_files/figure-markdown_github/plot_tstar-3.png)![](pbr28_files/figure-markdown_github/plot_tstar-4.png)

All 4 tstar plots on a single page
----------------------------------

Note: I would like to have the “draw\_figure\_label” be drawn a bit
above the first PET, but that does not seem to be happening. Should I
just skip it? I Placed the plots in a single column because it was a
little hard distinguishing where one ands and another begins in a 2 x 2
format. Perhaps one could put borders and then the 2 x 2 would work
better?

``` r
plot_grid(plotlist = tstar_fits$logan_tstar, ncol = 2, nrow = 1, labels = paste('PET:',tstar_fits$PET), label_x = 0.5, label_y = 1, label_colour = "blue", vjust = 0.97) +
  draw_figure_label("t*", position = "top", fontface = "bold", size = 32, colour = "red")
```

Define functions for fitting the models
=======================================

``` r
# MA1 using the fitted delay and vB from delayFit
fitma1 <- function(tacdata, input, delayFit) {
  ma1(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, tstarIncludedFrames = 6,
      inpshift = delayFit$par$inpshift, weights=tacdata$weights, vB= 0, frameStartEnd = c(1,33))
}

# 2TCM using the fitted delay and vB from delayFit
fit2tcm <- function(tacdata, input, delayFit) {
  twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, 
         inpshift = delayFit$par$inpshift, vB= 0.05, weights=tacdata$weights, frameStartEnd = c(1,33), multstart_iter = 100)
}

#Loganplot
fit_Logan <- function(tacdata, input, delayFit) {
  Loganplot(t_tac = tacdata$Times/60, tac = tacdata$TAC,
         input = input, inpshift = delayFit$par$inpshift, tstarIncludedFrames = 10, vB= 0,frameStartEnd = c(1,33))
}
```

Fit kinetic models
==================

``` r
set.seed(123)

tacs_long <- tacs_long %>% 

  # 2TCM using fitted vB and delay
  mutate(fit_2tcm= pmap(list(tacdata, input, delayFit), fit2tcm)) %>%
  mutate(Vt_2tcm = map_dbl(fit_2tcm, c('par', 'Vt'))) %>% 
  
  # MA1
  mutate(fit_ma1 = pmap(list(tacdata, input, delayFit), fitma1))  %>%
  mutate(Vt_ma1 = map_dbl(fit_ma1, c('par', 'Vt'))) %>%
   
  #Loganplot
  mutate(Loganfit = pmap(list(tacdata, input, delayFit), fit_Logan)) %>% 
  mutate(Vt_Logan = map_dbl(Loganfit, c("par", "Vt")))
```

\#Plot kinetic models

Plot 2tcm
---------

``` r
PETs <- unique(tacs_long$PET)

delayFits <- map(tacs_long$delayFit[tacs_long$Region=='WB'], 
                ~plot_inptac_fit(.x) + ggtitle('Delay')) 

delayFits  <- data.frame(PET = unique(tacs_long$PET)) %>% 
              mutate(fit = delayFits) %>% 
              filter(PET == PETs[5:8])
  

plot_2tcm <- tacs_long %>% 
  filter(PET == PETs[5:8]) %>% 
 group_by(PET, Region) %>% 
  mutate(fit = map2(fit_2tcm, Region, 
      ~ plot_kinfit(.x, roiname = .y))) %>% 
  ungroup() %>% 
  filter(Region %in% c('FC', 'WB', 'STR', 'CBL', 'THA')) %>% 
  select(PET, fit) %>%
  group_by(PET) %>% 
  arrange(PET) %>% 
  bind_rows(delayFits)


walk2(list(plot_2tcm$fit), unique(plot_2tcm$PET), 
    ~print(plot_grid( plotlist = .x, ncol = 2, nrow = 3, align = 'hv') +
  draw_plot_label( label = .y , y = 1, x = 0.5, fontface = "bold", size = 20, colour = "red", hjust = 1, vjust = 1.2)))
```

![](pbr28_files/figure-markdown_github/plot_2tcmfit-1.png)![](pbr28_files/figure-markdown_github/plot_2tcmfit-2.png)![](pbr28_files/figure-markdown_github/plot_2tcmfit-3.png)![](pbr28_files/figure-markdown_github/plot_2tcmfit-4.png)

\#\#\#subplot version

``` r
plot.TAC <- plot_kinfit(tacs_long$fit_2tcm[[4]], roiname = tacs_long$Region[4])

plot.AIF <- plot_kinfit(tacs_long$fit_2tcm[[4]], roiname = tacs_long$Region[4]) +
coord_cartesian(ylim=c(0,3000))

p <- ggdraw() +
  draw_plot(plot.TAC + theme(legend.position = "none"), 0, 0, 1, 1) +
  draw_plot(plot.AIF  + 
              theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()), 0.5, 0.6, 0.5, 0.4) 

legend <- get_legend(plot.TAC)

q <- plot_grid( legend, p , rel_widths = c(.35, 2))

print(q)
```

new plot 2tcm
=============

``` r
k <- tacs_long %>% 
  group_by(PET, Region) %>% 
  mutate(twotcm = map(fit_2tcm, c("tacs"))) %>% 
   select(PET, Region, twotcm) %>% 
  filter(Region %in% c('FC', 'WB', 'STR', 'CBL')) %>% 
  unnest()
         

ggplot(k, aes(x=Time, y=Target, color = Region)) +
  geom_point() + geom_line(aes(y=Target_fitted, color = Region)) + 
  facet_wrap(~ PET , ncol=2)
```

Plot ma1
--------

``` r
plot_MA1 <- tacs_long %>%
  filter(PET == PETs[5:8]) %>%
 group_by(PET, Region) %>% 
  mutate(fit = map2(fit_ma1, Region, 
      ~ plot_kinfit(.x, roiname = .y))) %>% 
  ungroup() %>% 
  filter(Region %in% c('FC', 'WB', 'STR', 'CBL', 'THA')) %>% 
  select(PET, fit) %>%
  group_by(PET) %>% 
  arrange(PET) %>% 
  bind_rows(delayFits)


walk2(list(plot_MA1$fit), unique(plot_MA1$PET), 
    ~print(plot_grid( plotlist = .x, ncol = 2, nrow = 3, align = 'hv') +
  draw_plot_label( label = .y , y = 1, x = 0.5, fontface = "bold", size = 20, colour = "red", hjust = 1, vjust = 1.2)))
```

![](pbr28_files/figure-markdown_github/plot_MA1-1.png)![](pbr28_files/figure-markdown_github/plot_MA1-2.png)![](pbr28_files/figure-markdown_github/plot_MA1-3.png)![](pbr28_files/figure-markdown_github/plot_MA1-4.png)

Plot Loganplot
--------------

``` r
plot_Logan <- tacs_long %>% 
  filter(PET == PETs[5:8]) %>%
 group_by(PET, Region) %>% 
  mutate(fit = map2(Loganfit, Region, 
      ~ plot(.x, roiname = .y))) %>% 
  ungroup() %>% 
  filter(Region %in% c('FC', 'WB', 'CBL', 'THA')) %>% 
  select(PET, fit) %>%
  group_by(PET) %>% 
  arrange(PET) 
  


walk2(list(plot_Logan$fit), unique(plot_Logan$PET), 
    ~print(plot_grid( plotlist = .x, ncol = 2, nrow = 2, align = 'hv') +
  draw_plot_label( label = .y , y = 1, x = 0.5, fontface = "bold", size = 20, colour = "red", hjust = 1, vjust = 1.2)))
```

![](pbr28_files/figure-markdown_github/Plot_Logan-1.png)![](pbr28_files/figure-markdown_github/Plot_Logan-2.png)![](pbr28_files/figure-markdown_github/Plot_Logan-3.png)

Test-retest
===========

trt preparation
---------------

``` r
trt_check <- tacs_long %>%
  select(Subjname, PETNo, Region, Vt_ma1, Vt_2tcm, Vt_Logan) %>% 
  gather(Measure, Value, -Subjname, -PETNo, -Region) %>%
  group_by(Region, Measure) %>% 
  nest(.key = "data")

saveRDS(tacs_long,'DerivedData/raw_kinfit_pbr28.rds')
```

trt results
-----------

``` r
trt_check <- trt_check %>% 
  group_by(Region, Measure) %>% 
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "Subjname")),
         trt_tidy = map(trt, c("tidy")))

trt_table <- select(trt_check, trt_tidy) %>% 
  unnest() %>% 
  ungroup() %>% 
  select(Region:wscv)

measure <- unique(trt_table$Measure) 

trt_table %>% 
  select(-Measure) %>% 
kable(digits=3) 
```

<table>
<thead>
<tr>
<th style="text-align:left;">
Region
</th>
<th style="text-align:right;">
mean
</th>
<th style="text-align:right;">
sd
</th>
<th style="text-align:right;">
cov
</th>
<th style="text-align:right;">
skew
</th>
<th style="text-align:right;">
kurtosis
</th>
<th style="text-align:right;">
icc
</th>
<th style="text-align:right;">
icc\_l
</th>
<th style="text-align:right;">
icc\_u
</th>
<th style="text-align:right;">
wscv
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
FC
</td>
<td style="text-align:right;">
3.131
</td>
<td style="text-align:right;">
1.499
</td>
<td style="text-align:right;">
0.479
</td>
<td style="text-align:right;">
1.243
</td>
<td style="text-align:right;">
0.977
</td>
<td style="text-align:right;">
0.913
</td>
<td style="text-align:right;">
0.736
</td>
<td style="text-align:right;">
0.974
</td>
<td style="text-align:right;">
0.145
</td>
</tr>
<tr>
<td style="text-align:left;">
TC
</td>
<td style="text-align:right;">
3.157
</td>
<td style="text-align:right;">
1.544
</td>
<td style="text-align:right;">
0.489
</td>
<td style="text-align:right;">
1.311
</td>
<td style="text-align:right;">
1.342
</td>
<td style="text-align:right;">
0.927
</td>
<td style="text-align:right;">
0.777
</td>
<td style="text-align:right;">
0.978
</td>
<td style="text-align:right;">
0.135
</td>
</tr>
<tr>
<td style="text-align:left;">
STR
</td>
<td style="text-align:right;">
3.026
</td>
<td style="text-align:right;">
1.465
</td>
<td style="text-align:right;">
0.484
</td>
<td style="text-align:right;">
1.220
</td>
<td style="text-align:right;">
0.982
</td>
<td style="text-align:right;">
0.906
</td>
<td style="text-align:right;">
0.720
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
0.151
</td>
</tr>
<tr>
<td style="text-align:left;">
THA
</td>
<td style="text-align:right;">
4.050
</td>
<td style="text-align:right;">
2.239
</td>
<td style="text-align:right;">
0.553
</td>
<td style="text-align:right;">
1.204
</td>
<td style="text-align:right;">
0.831
</td>
<td style="text-align:right;">
0.905
</td>
<td style="text-align:right;">
0.715
</td>
<td style="text-align:right;">
0.971
</td>
<td style="text-align:right;">
0.174
</td>
</tr>
<tr>
<td style="text-align:left;">
WB
</td>
<td style="text-align:right;">
3.050
</td>
<td style="text-align:right;">
1.469
</td>
<td style="text-align:right;">
0.481
</td>
<td style="text-align:right;">
1.301
</td>
<td style="text-align:right;">
1.192
</td>
<td style="text-align:right;">
0.905
</td>
<td style="text-align:right;">
0.715
</td>
<td style="text-align:right;">
0.971
</td>
<td style="text-align:right;">
0.152
</td>
</tr>
<tr>
<td style="text-align:left;">
CBL
</td>
<td style="text-align:right;">
3.355
</td>
<td style="text-align:right;">
1.680
</td>
<td style="text-align:right;">
0.501
</td>
<td style="text-align:right;">
1.352
</td>
<td style="text-align:right;">
1.367
</td>
<td style="text-align:right;">
0.922
</td>
<td style="text-align:right;">
0.759
</td>
<td style="text-align:right;">
0.977
</td>
<td style="text-align:right;">
0.143
</td>
</tr>
<tr>
<td style="text-align:left;">
FC
</td>
<td style="text-align:right;">
3.083
</td>
<td style="text-align:right;">
1.496
</td>
<td style="text-align:right;">
0.485
</td>
<td style="text-align:right;">
0.922
</td>
<td style="text-align:right;">
0.214
</td>
<td style="text-align:right;">
0.899
</td>
<td style="text-align:right;">
0.699
</td>
<td style="text-align:right;">
0.969
</td>
<td style="text-align:right;">
0.158
</td>
</tr>
<tr>
<td style="text-align:left;">
TC
</td>
<td style="text-align:right;">
3.028
</td>
<td style="text-align:right;">
1.517
</td>
<td style="text-align:right;">
0.501
</td>
<td style="text-align:right;">
1.053
</td>
<td style="text-align:right;">
0.659
</td>
<td style="text-align:right;">
0.917
</td>
<td style="text-align:right;">
0.750
</td>
<td style="text-align:right;">
0.975
</td>
<td style="text-align:right;">
0.147
</td>
</tr>
<tr>
<td style="text-align:left;">
STR
</td>
<td style="text-align:right;">
8.798
</td>
<td style="text-align:right;">
27.297
</td>
<td style="text-align:right;">
3.102
</td>
<td style="text-align:right;">
4.266
</td>
<td style="text-align:right;">
17.026
</td>
<td style="text-align:right;">
0.010
</td>
<td style="text-align:right;">
-0.545
</td>
<td style="text-align:right;">
0.560
</td>
<td style="text-align:right;">
3.087
</td>
</tr>
<tr>
<td style="text-align:left;">
THA
</td>
<td style="text-align:right;">
3.907
</td>
<td style="text-align:right;">
2.196
</td>
<td style="text-align:right;">
0.562
</td>
<td style="text-align:right;">
1.067
</td>
<td style="text-align:right;">
0.401
</td>
<td style="text-align:right;">
0.910
</td>
<td style="text-align:right;">
0.728
</td>
<td style="text-align:right;">
0.973
</td>
<td style="text-align:right;">
0.173
</td>
</tr>
<tr>
<td style="text-align:left;">
WB
</td>
<td style="text-align:right;">
2.919
</td>
<td style="text-align:right;">
1.437
</td>
<td style="text-align:right;">
0.492
</td>
<td style="text-align:right;">
1.031
</td>
<td style="text-align:right;">
0.455
</td>
<td style="text-align:right;">
0.902
</td>
<td style="text-align:right;">
0.707
</td>
<td style="text-align:right;">
0.970
</td>
<td style="text-align:right;">
0.158
</td>
</tr>
<tr>
<td style="text-align:left;">
CBL
</td>
<td style="text-align:right;">
3.054
</td>
<td style="text-align:right;">
1.631
</td>
<td style="text-align:right;">
0.534
</td>
<td style="text-align:right;">
1.227
</td>
<td style="text-align:right;">
0.926
</td>
<td style="text-align:right;">
0.909
</td>
<td style="text-align:right;">
0.716
</td>
<td style="text-align:right;">
0.973
</td>
<td style="text-align:right;">
0.165
</td>
</tr>
<tr>
<td style="text-align:left;">
FC
</td>
<td style="text-align:right;">
2.940
</td>
<td style="text-align:right;">
1.471
</td>
<td style="text-align:right;">
0.501
</td>
<td style="text-align:right;">
1.173
</td>
<td style="text-align:right;">
0.827
</td>
<td style="text-align:right;">
0.917
</td>
<td style="text-align:right;">
0.741
</td>
<td style="text-align:right;">
0.975
</td>
<td style="text-align:right;">
0.147
</td>
</tr>
<tr>
<td style="text-align:left;">
TC
</td>
<td style="text-align:right;">
2.933
</td>
<td style="text-align:right;">
1.480
</td>
<td style="text-align:right;">
0.505
</td>
<td style="text-align:right;">
1.259
</td>
<td style="text-align:right;">
1.198
</td>
<td style="text-align:right;">
0.922
</td>
<td style="text-align:right;">
0.756
</td>
<td style="text-align:right;">
0.977
</td>
<td style="text-align:right;">
0.144
</td>
</tr>
<tr>
<td style="text-align:left;">
STR
</td>
<td style="text-align:right;">
2.816
</td>
<td style="text-align:right;">
1.426
</td>
<td style="text-align:right;">
0.506
</td>
<td style="text-align:right;">
1.156
</td>
<td style="text-align:right;">
0.830
</td>
<td style="text-align:right;">
0.906
</td>
<td style="text-align:right;">
0.719
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
0.158
</td>
</tr>
<tr>
<td style="text-align:left;">
THA
</td>
<td style="text-align:right;">
3.794
</td>
<td style="text-align:right;">
2.062
</td>
<td style="text-align:right;">
0.544
</td>
<td style="text-align:right;">
1.182
</td>
<td style="text-align:right;">
0.781
</td>
<td style="text-align:right;">
0.907
</td>
<td style="text-align:right;">
0.722
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
0.169
</td>
</tr>
<tr>
<td style="text-align:left;">
WB
</td>
<td style="text-align:right;">
2.822
</td>
<td style="text-align:right;">
1.401
</td>
<td style="text-align:right;">
0.497
</td>
<td style="text-align:right;">
1.254
</td>
<td style="text-align:right;">
1.098
</td>
<td style="text-align:right;">
0.907
</td>
<td style="text-align:right;">
0.718
</td>
<td style="text-align:right;">
0.972
</td>
<td style="text-align:right;">
0.155
</td>
</tr>
<tr>
<td style="text-align:left;">
CBL
</td>
<td style="text-align:right;">
3.103
</td>
<td style="text-align:right;">
1.603
</td>
<td style="text-align:right;">
0.517
</td>
<td style="text-align:right;">
1.321
</td>
<td style="text-align:right;">
1.290
</td>
<td style="text-align:right;">
0.915
</td>
<td style="text-align:right;">
0.736
</td>
<td style="text-align:right;">
0.975
</td>
<td style="text-align:right;">
0.154
</td>
</tr>
</tbody>
</table>

Interregional Correlation
=========================

Here the interregional correlations for V<sub>T</sub> are assessed

``` r
Vt_2TCM <- tacs_long %>%
  ungroup() %>% 
  select(PET, Region, Vt_2tcm) %>%
  spread(Region, Vt_2tcm)

Vt_MA1 <- tacs_long %>%
  ungroup() %>%
  select(PET, Region, Vt_ma1) %>%
  spread(Region, Vt_ma1)

Vt_logan <- tacs_long %>%
  ungroup() %>%
  select(PET, Region, Vt_Logan) %>%
  spread(Region, Vt_Logan)

col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))

par(mfrow=c(2,2))

Vt_2TCM %>%
  select(FC:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Vt_2TCM ~ Correlations),
                 mar=c(0,0,1,0))

Vt_logan %>%
  select(FC:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Vt_logan ~ Correlations),
                 mar=c(0,0,1,0))

Vt_MA1 %>%
  select(FC:WB) %>%
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Vt_MA1 ~ Correlations),
                 mar=c(0,0,1,0))
```

![](pbr28_files/figure-markdown_github/corrplots-1.png)

\#Corrplot between measures for a single region

``` r
compare <- tacs_long %>%
  ungroup() %>%
  select(PET, Region, Vt_2tcm, Vt_Logan ,Vt_ma1 ) %>% 
  filter(Region %in% c('FC', 'WB', 'STR', 'THA'))

par(mfrow=c(2,2))

compare %>%
  filter(Region == "FC") %>%
  select(Vt_2tcm:Vt_ma1) %>% 
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Model ~ Correlations ~ Region: FC),
                 mar=c(0,0,1,0)) 

compare %>%
  filter(Region == "THA") %>%
  select(Vt_2tcm:Vt_ma1) %>% 
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Model ~ Correlations ~ Region: OC), mar=c(0,0,1,0))  

compare %>%
  filter(Region == "STR") %>%
  select(Vt_2tcm:Vt_ma1) %>% 
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Model ~ Correlations ~ Region: STR),
                 mar=c(0,0,1,0))  

compare %>%
  filter(Region == "WB") %>%
  select(Vt_2tcm:Vt_ma1) %>% 
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(Model ~ Correlations ~ Region: WB),
                 mar=c(0,0,1,0)) 
```

![](pbr28_files/figure-markdown_github/correlation_plots-1.png)

Vt corellation
==============

R-squared
---------

``` r
trtdata <- tacs_long %>%
  select(PET, Subjname, PETNo, Region, Vt_2tcm, Vt_Logan ,Vt_ma1) %>%
  gather(Measure, Value, -(PET:Region)) %>%
  spread(Region, Value)

trtdata <- trtdata %>%
  gather(Region, Value, -(PET:Measure)) %>%
  unite(Outcome, Measure, Region) %>%
  spread(Outcome, Value)

corout <- trtdata %>%
  gather(Measure, Binding, -(PET:PETNo), -Vt_2tcm_WB) %>%
  group_by(Measure) %>%
  summarise('R^2^'=cor(Binding, Vt_2tcm_WB)^2) %>%
  arrange(Measure) %>%
  ungroup() %>%
  mutate(Measure = str_replace(string=Measure, pattern='_', replacement='~')) %>%
  mutate(Measure = str_replace(string=Measure, pattern='FC', replacement='FC~')) %>%
   mutate(Measure = str_replace(string=Measure, pattern='CBL', replacement='CBL~')) %>%
   mutate(Measure = str_replace(string=Measure, pattern='ACC', replacement='ACC~')) %>%
   mutate(Measure = str_replace(string=Measure, pattern='INS', replacement='INS~')) %>%
   mutate(Measure = str_replace(string=Measure, pattern='THA', replacement='THA~')) %>%
  mutate(Measure = str_replace(string=Measure, pattern='WB', replacement='WB~')) %>%
  mutate(Measure = str_replace(string=Measure, pattern='OC', replacement='OC~')) %>% 
  mutate(Measure = str_replace(string=Measure, pattern='WM', replacement='WM~')) %>% 
  mutate(Measure = str_replace(string=Measure, pattern='GM', replacement='GM~')) %>% 
  mutate(Measure = str_replace(string=Measure, pattern='STR', replacement='STR~')) 

kable(corout, digits=2, caption="Correlations with BP_srtm~WB~")
```

<table>
<caption>
Correlations with BP\_srtm<sub>WB</sub>
</caption>
<thead>
<tr>
<th style="text-align:left;">
Measure
</th>
<th style="text-align:right;">
R<sup>2</sup>
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Vt<sub>2tcm\_CBL</sub>
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>2tcm\_FC</sub>
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>2tcm\_STR</sub>
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt\~2tcm\_TC
</td>
<td style="text-align:right;">
0.99
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>2tcm\_THA</sub>
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>Logan\_CBL</sub>
</td>
<td style="text-align:right;">
0.96
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>Logan\_FC</sub>
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>Logan\_STR</sub>
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt\~Logan\_TC
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>Logan\_THA</sub>
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>Logan\_WB</sub>
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>ma1\_CBL</sub>
</td>
<td style="text-align:right;">
0.95
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>ma1\_FC</sub>
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>ma1\_STR</sub>
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt\~ma1\_TC
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>ma1\_THA</sub>
</td>
<td style="text-align:right;">
0.97
</td>
</tr>
<tr>
<td style="text-align:left;">
Vt<sub>ma1\_WB</sub>
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
</tbody>
</table>

\#Plot of the change between PETNo = 1 and PETNo = 2.

``` r
trtdata <- trtdata %>% 
  gather(Region, Value, -(PET:PETNo)) %>% 
  separate(col = "Region", into = c("outcome","Measure", "Region"), sep = '_')%>% unite(outcome, outcome, Measure, sep = '_') %>% 
  ungroup()

trt_2tcm <- trtdata %>% 
  filter(outcome == 'Vt_2tcm')

ggplot(trt_2tcm, aes(x = PETNo, y = Value, 
                         group = Region, colour=Region)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  facet_wrap( ~ Subjname, ncol = 4)  
```

![](pbr28_files/figure-markdown_github/PET_change-1.png)
