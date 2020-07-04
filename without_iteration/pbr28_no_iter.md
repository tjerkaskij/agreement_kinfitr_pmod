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
install.packages("epitrix")
```

First, the libraries for the analysis and plotting are loaded.

``` r
library(stringr)
library(corrplot)
```

    ## corrplot 0.84 loaded

``` r
library(grid)
library(gridExtra)
library(RColorBrewer)
library(psych)
```

    ## Warning: package 'psych' was built under R version 3.6.3

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
library(kableExtra)
library(tidyverse)
```

    ## Warning: package 'tidyverse' was built under R version 3.6.3

    ## -- Attaching packages ------------------------------------ tidyverse 1.3.0 --

    ## v ggplot2 3.3.0     v purrr   0.3.3
    ## v tibble  3.0.0     v dplyr   0.8.5
    ## v tidyr   1.0.2     v forcats 0.5.0
    ## v readr   1.3.1

    ## Warning: package 'ggplot2' was built under R version 3.6.3

    ## Warning: package 'tibble' was built under R version 3.6.3

    ## Warning: package 'purrr' was built under R version 3.6.3

    ## Warning: package 'dplyr' was built under R version 3.6.3

    ## Warning: package 'forcats' was built under R version 3.6.3

    ## -- Conflicts --------------------------------------- tidyverse_conflicts() --
    ## x ggplot2::%+%()      masks psych::%+%()
    ## x ggplot2::alpha()    masks psych::alpha()
    ## x dplyr::combine()    masks gridExtra::combine()
    ## x purrr::cross()      masks pracma::cross()
    ## x tidyr::expand()     masks Matrix::expand()
    ## x dplyr::filter()     masks stats::filter()
    ## x dplyr::group_rows() masks kableExtra::group_rows()
    ## x dplyr::lag()        masks stats::lag()
    ## x tidyr::pack()       masks Matrix::pack()
    ## x tidyr::unpack()     masks Matrix::unpack()

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
  arrange(PET) %>% 
  mutate(Subjname = hash_names(Subjname, full = F)) %>% #change subject names
  mutate(PET = paste(Subjname, PETNo, sep='_'))

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

    ## Registered S3 methods overwritten by 'car':
    ##   method                          from
    ##   influence.merMod                lme4
    ##   cooks.distance.influence.merMod lme4
    ##   dfbeta.influence.merMod         lme4
    ##   dfbetas.influence.merMod        lme4

    ## Warning in twotcm(t_tac = .x$Times/60, tac = .x$WB, input = .y, inpshift.upper = 1, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

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

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 5 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 4 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

\#Rearrangement of the Data into Long Format

``` r
tacs_long <- tacs %>%
  unnest(tacdata, .drop = FALSE) %>%
  select(-Weights) 
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## All list-columns are now preserved.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

``` r
tacs_long$weights = weights_create(t_start = tacs_long$StartTime/60, t_end = (tacs_long$StartTime + tacs_long$Duration)/60, tac = tacs_long$WB, radioisotope = "C11")

tacs_long <- tacs_long %>% 
gather(Region, TAC, FC:CBL) %>%
group_by(PET, Subjname, PETNo, Region) %>%
  filter(Region %in% c("FC","THA")) %>% 
nest(c(Times, StartTime, Duration, weights, TAC),
     .key = 'tacdata') 
```

    ## Warning: All elements of `...` must be named.
    ## Did you want `tacdata = c(Times, StartTime, Duration, weights, TAC)`?

``` r
# tacs_long <- tacs %>%
#   select(-tacdata) %>%
#   inner_join(tacs_long, by = c("PET", "Subjname", "PETNo"))
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

![](pbr28_no-iter_files/figure-markdown_github/plot_tstar-1.png)![](pbr28_no-iter_files/figure-markdown_github/plot_tstar-2.png)![](pbr28_no-iter_files/figure-markdown_github/plot_tstar-3.png)![](pbr28_no-iter_files/figure-markdown_github/plot_tstar-4.png)

All 4 tstar plots on a single page
----------------------------------

``` r
plot_grid(plotlist = tstar_fits$logan_tstar, ncol = 2, nrow = 1, labels = paste('PET:',tstar_fits$PET), label_x = 0.5, label_y = 1, label_colour = "blue", vjust = 0.97) +
  draw_figure_label("t*", position = "top", fontface = "bold", size = 32, colour = "red")
```

Define functions for fitting the models
=======================================

``` r
# MA1 using the fitted delay and vB from delayFit
fitma1 <- function(tacdata, input, delayFit) {
  ma1(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, 
      tstarIncludedFrames = 6,
      inpshift = delayFit$par$inpshift, weights=tacdata$weights, vB= 0, 
      frameStartEnd = c(1,33), dur = tacdata$Duration/60)
}

# 2TCM using the fitted delay and vB from delayFit
fit2tcm <- function(tacdata, input, delayFit) {
  twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, 
         inpshift = delayFit$par$inpshift, vB= 0.05, weights=tacdata$weights, frameStartEnd = c(1,33))
}

#Loganplot
fit_Logan <- function(tacdata, input, delayFit) {
  Loganplot(t_tac = tacdata$Times/60, tac = tacdata$TAC,
         input = input, inpshift = delayFit$par$inpshift, tstarIncludedFrames = 10, vB= 0,frameStartEnd = c(1,33), dur = tacdata$Duration/60)
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

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

    ## Warning in twotcm(t_tac = tacdata$Times/60, tac = tacdata$TAC, input = input, : Fitted parameters are hitting upper or lower limit bounds. Consider 
    ## either modifying the upper and lower limit boundaries, or else using 
    ## multstart when fitting the model (see the function documentation).

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

    ## Adding missing grouping variables: `PET`

    ## Warning: `.key` is deprecated

``` r
saveRDS(tacs_long,'DerivedData/raw_kinfit_pbr28_no_iter.rds')
```
