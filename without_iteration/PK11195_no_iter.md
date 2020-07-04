-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
    -   [Creating folders using
        “initProjectFolder()”](#creating-folders-using-initprojectfolder)
-   [Extracting roistats and Tidying
    data.](#extracting-roistats-and-tidying-data.)
-   [Fitting of the Delay and Blood Volume
    Fraction](#fitting-of-the-delay-and-blood-volume-fraction)
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

Aims
====

The aim of this assignment is to analyze the PK11195 data in kinfitr
without the function used to iterate across multiple starting parameters
when fitting 2TCM

Libraries
=========

CRAN libraries
--------------

First, the libraries for the analysis and plotting are loaded.

``` r
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.6.3

    ## -- Attaching packages ------------------------------------ tidyverse 1.3.0 --

    ## v ggplot2 3.3.0     v purrr   0.3.3
    ## v tibble  3.0.0     v dplyr   0.8.5
    ## v tidyr   1.0.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## Warning: package 'ggplot2' was built under R version 3.6.3

    ## Warning: package 'tibble' was built under R version 3.6.3

    ## Warning: package 'purrr' was built under R version 3.6.3

    ## Warning: package 'dplyr' was built under R version 3.6.3

    ## Warning: package 'forcats' was built under R version 3.6.3

    ## -- Conflicts --------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(stringr)
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
library(grid)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(RColorBrewer)
library(psych)
```

    ## Warning: package 'psych' was built under R version 3.6.3

    ## 
    ## Attaching package: 'psych'

    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, alpha

``` r
library(readxl)
library(pracma)
```

    ## Warning: package 'pracma' was built under R version 3.6.3

    ## 
    ## Attaching package: 'pracma'

    ## The following objects are masked from 'package:psych':
    ## 
    ##     logit, polar

    ## The following object is masked from 'package:purrr':
    ## 
    ##     cross

``` r
library(lme4)
```

    ## Warning: package 'lme4' was built under R version 3.6.3

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:pracma':
    ## 
    ##     expm, lu, tril, triu

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(rjags)
```

    ## Warning: package 'rjags' was built under R version 3.6.3

    ## Loading required package: coda

    ## Warning: package 'coda' was built under R version 3.6.3

    ## Linked to JAGS 4.3.0

    ## Loaded modules: basemod,bugs

``` r
library(knitr)
```

    ## Warning: package 'knitr' was built under R version 3.6.3

``` r
library(cowplot)
```

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

``` r
library(corrplot)
library(ggplotify)
```

    ## Warning: package 'ggplotify' was built under R version 3.6.3

``` r
library(epitrix)
```

    ## Warning: package 'epitrix' was built under R version 3.6.3

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
# library(vmisc)
library(kipettools)
library(granviller)
```

    ## 
    ## Attaching package: 'granviller'

    ## The following object is masked from 'package:pracma':
    ## 
    ##     cd

``` r
library(relfeas)
```

    ## 
    ## Attaching package: 'relfeas'

    ## The following object is masked from 'package:granviller':
    ## 
    ##     trt

``` r
theme_set(theme_light())
```

Creating folders using “initProjectFolder()”
--------------------------------------------

``` r
initProjectFolder()
```

Extracting roistats and Tidying data.
=====================================

``` r
#Extracting tac data
tacs <- tibble(Filename = list.files(path = "../RawData/", 
                                          pattern = "roistats.mat")) %>% 
   group_by(Filename) %>% 
   mutate(tacdata = map(Filename, ~kipettools::roistats_getData(
     paste0("../RawData/", .x))))
tacs <- tacs %>% 
  ungroup() %>% 
  mutate(Subjname = map_chr(tacdata, "Subjname"),
         PETNo = map_dbl(tacdata, "PETNo"),
         tacdata = map(tacdata, "tacdata")) %>% 
  select(-Filename) %>% 
  mutate(PET = paste(Subjname, PETNo, sep='_')) 
#loading weights and blood data
oldwd <- getwd()
setwd("../RawData/")
blood <- list()
bloodfiles <- list.files(pattern='blood_processed_pfhill')
for(i in 1:length(bloodfiles)) {
blood[[bloodfiles[i]]] <- read_tsv(bloodfiles[i])
   print(paste0('Progress: ', i, ' / ', length(bloodfiles)))
 }
blood <- tibble(blood)
Weights_list <- list()
weightfiles <- list.files(pattern='_weights2009')
for(i in 1:length(weightfiles)) {
Weights_list[[weightfiles[i]]] <- read_csv(weightfiles[i], col_names = FALSE)
   print(paste0('Progress: ', i, ' / ', length(weightfiles)))
 }
setwd(oldwd)
Weights_list <- tibble(Weights_list)
Weights_list$PET <- weightfiles %>%  
str_replace("_weights2009.txt", " ") 
Weights_list <- Weights_list %>% 
group_by(PET) %>% 
unnest() %>% 
rename(weights = X1) %>% 
ungroup()

#Note: new weights added

tacs <- tacs %>% 
unnest() %>% 
rename(Times = times) %>% 
select(Subjname:PET, Times, durations, WM, GM, WB, FC, OC, THA, STR, TC, ACC, CBL = CER, INS) %>% 
add_column(weights = Weights_list$weights) %>% 
select(-weights) %>% 
mutate("start" = Times - (durations/2)) %>% 
mutate("end" = Times + (durations/2)) 

tacs$weights <- weights_create(t_start = tacs$start, t_end = tacs$end, 
                               tac = tacs$WB, radioisotope = "C11")

tacs <- tacs %>% 
group_by(PET, Subjname, PETNo) %>% 
nest(.key = tacdata) %>% 
ungroup()

blood$PET. <- bloodfiles %>%  
str_replace("_blood_processed_pfhill.txt", " ") 
blood  <- blood %>% 
  unnest() %>% 
 rename(Cbl.disp.corr = "Cbl disp corr", Cpl = "Cpl (nCi/cc)", ABSS.sec = "ABSS sec" ) %>% 
  mutate(Cbl.disp.corr = ifelse(Cbl.disp.corr < 0, 0, Cbl.disp.corr)) %>%
  group_by(PET.) %>%
  nest(.key='blooddata') %>%
  ungroup() %>% 
   mutate(input = map(blooddata, ~blood_interp(
             t_blood = .x$ABSS.sec/60, blood =.x$Cbl.disp.corr, 
             t_plasma=.x$ABSS.sec/60, plasma =.x$Cpl, 
             t_parentfrac = .x$ABSS.sec/60 , parentfrac= .x$parent_fract ) )) 

tacs <- tacs %>% 
  arrange(PET, Subjname, PETNo) %>% 
  bind_cols(blood)%>% 
  select(-PET.)

tacs <- tacs %>% 
  mutate(Subjname = hash_names(Subjname, full = F)) %>% 
  mutate(PET = paste(Subjname, PETNo, sep='_'))

saveRDS(tacs, '../DerivedData/tacs.rds')
```

Fitting of the Delay and Blood Volume Fraction
==============================================

``` r
tacs <- readRDS('../DerivedData/tacs.rds')
tacs <- tacs %>%
  group_by(PET, Subjname, PETNo) %>%
  mutate(delayFit = map2(tacdata,input, 
                         ~twotcm(t_tac = .x$Times, tac = .x$WB, input = .y,                                             inpshift.upper = 1))) %>%
  mutate(inpshift = map_dbl(delayFit, c("par", "inpshift"))) %>% 
  ungroup()
```

    ## Registered S3 methods overwritten by 'car':
    ##   method                          from
    ##   influence.merMod                lme4
    ##   cooks.distance.influence.merMod lme4
    ##   dfbeta.influence.merMod         lme4
    ##   dfbetas.influence.merMod        lme4

    ## Warning in twotcm(t_tac = .x$Times, tac = .x$WB, input = .y, inpshift.upper = 1): Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

creating logan\_tstar
=====================

``` r
logantstar <- function(tacdata, input, inpshift) {
  Logan_tstar(t_tac = tacdata$Times, input = input, lowroi = tacdata$FC, medroi =   tacdata$CBL, highroi = tacdata$THA, inpshift = inpshift, gridbreaks = 3)
}
tacs <-  tacs %>% 
  group_by(Subjname, PETNo) %>% 
  mutate(logan_tstar = pmap(list(tacdata, input, inpshift), logantstar))
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

\#Rearrangement of the Data into Long Format

``` r
tacs_long <- tacs %>%
  select(PET, Subjname, PETNo, tacdata, logan_tstar, inpshift) %>%
  unnest(tacdata, .drop = FALSE) %>%
  gather(Region, TAC, WM:INS) %>%
  filter(Region %in% c("FC", "THA")) %>% 
  group_by(PET, Subjname, PETNo, Region) %>%
  nest(.key = 'tacdata') 
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## All list-columns are now preserved.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

    ## Warning: `.key` is deprecated

``` r
tacs_long <- tacs %>%
  select(PET, Subjname, PETNo, input, delayFit) %>%
  inner_join(tacs_long, by = c("PET", "Subjname", "PETNo"))
```

Plotting logan\_tstar
---------------------

``` r
set.seed(123)
tstar_fits <- tacs_long %>%  
  ungroup() %>% 
  select(tacdata, PET) %>% 
  unnest() %>% 
  select(PET, logan_tstar) %>% 
  sample_n(size =  4, replace = F) 
walk2(list(tstar_fits$logan_tstar), tstar_fits$PET, 
    ~print(plot_grid(plotlist = .x, ncol = 1, nrow = 1, labels = paste('PET:',.y), label_x = 0.5, label_y = 1, label_colour = "blue", vjust = 0.97)))
```

All 4 tstar plots on a single page
----------------------------------

``` r
plot_grid(plotlist = tstar_fits$logan_tstar, ncol = 2, nrow = 1, labels = paste('PET:',tstar_fits$PET), label_x = 0.5, label_y = 1, label_colour = "blue", vjust = 0.97) +
  draw_figure_label("t*", position = "top", fontface = "bold", size = 32, colour = "red")
```

Define functions for fitting the models
=======================================

Comment: when the argumant “multstart\_iter” is equal to 1 for 2TCM, the
estimated microparameters of a number of study participants have a
tendency of approaching or even exceeding their default limits.

``` r
# MA1 using the fitted delay and vB from delayFit
fitma1 <- function(tacdata, input, delayFit) {
  ma1(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, tstarIncludedFrames = 6,
      inpshift = delayFit$par$inpshift, weights=tacdata$weights, vB= 0,
      dur = tacdata$durations)
}
# 2TCM using the fitted delay and vB from delayFit
fit2tcm <- function(tacdata, input, delayFit) {
  twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, 
         inpshift = delayFit$par$inpshift, vB= 0.05, weights=tacdata$weights)
}
#Loganplot
fit_Logan <- function(tacdata, input, delayFit) {
  Loganplot(t_tac = tacdata$Times, tac = tacdata$TAC,
         input = input, 
         inpshift = delayFit$par$inpshift, tstarIncludedFrames = 10, vB= 0,
         dur = tacdata$durations)
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
  mutate(Loganfit = pmap(list(tacdata, input, delayFit), 
                         fit_Logan)) %>% 
  mutate(Vt_Logan = map_dbl(Loganfit, c("par", "Vt")))
```

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

\#Plot kinetic models

Plot 2tcm
---------

``` r
delayFits <- map(tacs_long$delayFit[tacs_long$Region=='WB'], 
                ~plot_inptac_fit(.x) + ggtitle('Delay')) 

delayFits  <- data.frame(PET = unique(tacs_long$PET)) %>% 
              mutate(fit = delayFits)

plot_2tcm <- tacs_long %>% 
 group_by(PET, Region) %>% 
  mutate(fit = map2(fit_2tcm, Region, 
      ~ plot_kinfit(.x, roiname = .y))) %>% 
  ungroup() %>% 
  filter(Region %in% c('FC', 'WB', 'ACC', 'CBL', 'THA')) %>% 
  select(PET, fit) %>%
  group_by(PET) %>% 
  arrange(PET) %>% 
  bind_rows(delayFits)


walk2(list(plot_2tcm$fit), unique(plot_2tcm$PET), 
    ~print(plot_grid( plotlist = .x, ncol = 2, nrow = 3, align = 'hv') +
  draw_plot_label( label = .y , y = 1, x = 0.5, fontface = "bold", size = 20, colour = "red", hjust = 1, vjust = 1.2)))
```

\#\#\#subplot version

``` r
plot.TAC <- plot_kinfit(tacs_long$fit_2tcm[[4]], roiname = tacs_long$Region[4])

plot.AIF <- plot_kinfit(tacs_long$fit_2tcm[[4]], roiname = tacs_long$Region[4]) +
coord_cartesian(ylim=c(0,3000))

p <- ggdraw() +
  draw_plot(plot.TAC + theme(legend.position = "none"), 0, 0, 1, 1) +
  draw_plot(plot.AIF  + 
              theme(legend.position = "none"), 0.5, 0.52, 0.5, 0.4) +  
  draw_plot_label(c("A", "B"), c(0, 0.5), c(1, 0.92), size = 15)

legend <- get_legend(plot.TAC)

q <- plot_grid( legend, p , rel_widths = c(.35, 2))

ggsave(filename = "subplot.twotcm.2.pdf", plot = q, path = '../DerivedData/')
```

new plot 2tcm
=============

``` r
k <- tacs_long %>% 
  group_by(PET, Region) %>% 
  mutate(twotcm = map(fit_2tcm, c("tacs"))) %>% 
   select(PET, Region, twotcm) %>% 
  filter(Region %in% c('FC', 'WB', 'ACC', 'CBL')) %>% 
  unnest()
         

ggplot(k, aes(x=Time, y=Target, color = Region)) +
  geom_point() + geom_line(aes(y=Target_fitted, color = Region)) + 
  facet_wrap(~ PET , ncol=2)
```

Plot ma1
--------

``` r
plot_MA1 <- tacs_long %>% 
 group_by(PET, Region) %>% 
  mutate(fit = map2(fit_ma1, Region, 
      ~ plot_kinfit(.x, roiname = .y))) %>% 
  ungroup() %>% 
  filter(Region %in% c('FC', 'WB', 'ACC', 'CBL', 'THA')) %>% 
  select(PET, fit) %>%
  group_by(PET) %>% 
  arrange(PET) %>% 
  bind_rows(delayFits)


walk2(list(plot_MA1$fit), unique(plot_MA1$PET), 
    ~print(plot_grid( plotlist = .x, ncol = 2, nrow = 3, align = 'hv') +
  draw_plot_label( label = .y , y = 1, x = 0.5, fontface = "bold", size = 20, colour = "red", hjust = 1, vjust = 1.2)))
```

Plot Loganplot
--------------

``` r
plot_Logan <- tacs_long %>% 
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
```

    ## Warning: `.key` is deprecated

``` r
# saveRDS(tacs_long,'../DerivedData/raw_kinfit_pk11195_no_iter.rds')
```
