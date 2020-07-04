-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
    -   [Creating folders using
        “initProjectFolder()”](#creating-folders-using-initprojectfolder)
-   [Extracting kinfitresults](#extracting-kinfitresults)
-   [Tidying data.](#tidying-data.)
-   [New weights with
    kinfitr::weights\_create](#new-weights-with-kinfitrweights_create)
-   [Fitting and plotting MRTM1 multiple times for regions FC and
    WB](#fitting-and-plotting-mrtm1-multiple-times-for-regions-fc-and-wb)
-   [t\* values that were fitted by PMOD, in
    minutes](#t-values-that-were-fitted-by-pmod-in-minutes)
-   [Fitting MRTM2 to each region of each
    individual](#fitting-mrtm2-to-each-region-of-each-individual)
-   [save](#save)

Aims
====

The aim of this assignment is to analyze the AZ10419369 data in kinfitr
using the t\* values fitted by PMOD and constant weights.

Libraries
=========

CRAN libraries
--------------

First, the libraries for the analysis and plotting are loaded.

``` r
library(tidyverse)
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
library(viridis)
library(kableExtra)
library(epitrix)
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
# library(vmisc)
library(kipettools)
library(granviller)
library(relfeas)

theme_set(theme_light())
```

Creating folders using “initProjectFolder()”
--------------------------------------------

``` r
initProjectFolder()
```

Extracting kinfitresults
========================

``` r
 tactibble <- tibble(Filename = list.files(path = "../RawData/", 
                                           pattern = "kinfitresults.mat")) %>% 
  group_by(Filename) %>% 
  mutate(tacdata = map(Filename, ~kipettools::kfresults_getData(
    paste0("../RawData/", .x))))

saveRDS(tactibble, '../DerivedData/tactibble.rds')
```

Tidying data.
=============

``` r
#loading data

tacs <- readRDS('../DerivedData/tactibble.rds')

tacs <- tacs %>% 
  ungroup() %>% 
  mutate(Subjname = map_chr(tacdata, "Subjname"),
         PETNo = map_dbl(tacdata, "PETNo"),
         tacdata = map(tacdata, "tacdata")) %>% 
  select(-Filename) 
```

New weights with kinfitr::weights\_create
=========================================

``` r
#Creating start and end times and 
#removing the old weights, so that there is no mix-up later on

tacs <- tacs %>% 
  unnest() %>% 
  mutate("start" = times - (durations/2)) %>% 
  mutate("end" = times + (durations/2)) %>% 
  select(-weights)
```

    ## Warning: `cols` is now required.
    ## Please use `cols = c(tacdata)`

``` r
# New weights

tacs$weights <- weights_create(t_start = tacs$start, t_end = tacs$end, 
                               tac = tacs$WB, radioisotope = "C11")
#Nest

tacs <- tacs %>%
  nest(-c(Subjname, PETNo), .key = 'tacdata')
```

    ## Warning: All elements of `...` must be named.
    ## Did you want `tacdata = c(L_WM, R_WM, WM, GM, WB, gmRefCBL, gmL_fslLST, gmR_fslLST, gmfslLST, 
    ##     gmL_fslAST, gmR_fslAST, gmfslAST, gmL_fslSMST, gmR_fslSMST, 
    ##     gmfslSMST, gmL_fslSTR, gmR_fslSTR, gmfslSTR, gmfslFrontal_Lobe, 
    ##     gmfslInsula, gmL_fslAccumbens, gmL_fslAmygdala, gmL_fslCaudate, 
    ##     gmL_fslHippocampus, gmL_fslPallidum, gmL_fslPutamen, gmL_fslThalamus, 
    ##     gmfslOccipital_Lobe, gmfslParietal_Lobe, gmR_fslAccumbens, 
    ##     gmfslAccumbens, gmR_fslAmygdala, gmfslAmygdala, gmR_fslCaudate, 
    ##     gmfslCaudate, gmR_fslHippocampus, gmfslHippocampus, gmR_fslPallidum, 
    ##     gmfslPallidum, gmR_fslPutamen, gmfslPutamen, gmR_fslThalamus, 
    ##     gmfslThalamus, gmfslTemporal_Lobe, L_fslRetina, R_fslRetina, 
    ##     fslRetina, L_fslRetina_1, R_fslRetina_1, fslRetina_1, L_fslRetina_2, 
    ##     R_fslRetina_2, fslRetina_2, times, durations, start, end, 
    ##     weights)`?

Fitting and plotting MRTM1 multiple times for regions FC and WB
===============================================================

``` r
#K2 prime for striatum to be used to MRTM2 later on

 tacs <- tacs %>% 
  group_by(Subjname, PETNo) %>% 
  mutate(MRTM1fit = map(tacdata, ~mrtm1(t_tac = .x$times, reftac = .x$gmRefCBL,
                                   roitac = .x$gmfslSTR, 
                                      weights = .x$weights,
                                   dur = .x$durations))) %>% 
  ungroup() %>% 
  mutate(k2prime_MRTM1 = map_dbl(MRTM1fit, c("par", 
                                             "k2prime"))) %>%  mutate(MRTM1fit_weights = map(tacdata, ~mrtm1(t_tac = .x$times, reftac = .x$gmRefCBL,
                                   roitac = .x$gmfslSTR, 
                                   dur = .x$durations))) %>% 
  ungroup() %>% 
  mutate(k2prime_MRTM1_weights = map_dbl(MRTM1fit_weights, c("par", 
                                             "k2prime"))) %>%
  group_by(Subjname, PETNo) %>% 
  mutate(logan_tstar = map2(tacdata, k2prime_MRTM1,
                 ~refLogan_tstar(t_tac = .x$times, 
                    reftac = .x$gmRefCBL, lowroi = .x$gmL_fslThalamus, 
                    medroi = .x$gmfslFrontal_Lobe, highroi = .x$gmfslOccipital_Lobe, 
                   k2prime = .y)))
```

t\* values that were fitted by PMOD, in minutes
===============================================

``` r
pmod_t_stars <- readRDS("../DerivedData/pmod_macroparameters_t_star.rds") %>% 
  select(subjname:id,t_star, tracer, model, Region = region) %>% 
  filter(tracer == "az") %>% 
  select(-c(tracer,id)) %>% 
  rename(Subjname = subjname) %>% 
  mutate(PETNo = as.numeric(PETNo)) %>% 
  filter(model != "srtm") %>% 
  spread(key = model, value = t_star) %>%   rename_at(vars(ref_logan, mrtm2), ~ paste(., "_t_star", sep = ""))
```

Fitting MRTM2 to each region of each individual
===============================================

First, let’s select some specific regions. Note: I duplicated the CBL
region column into “CBL” and “Ref”. One of them is used to make the
reference tissue models when nested in tacdata whereas the other is used
for plotting in the new “all regions per PET +facet wrap by PET” - plot.
However, the srtm fitting gave me an error when I did this

``` r
# regions <- c("STR" = "gmfslSTR", "FC" = "gmfslFrontal_Lobe", "WB", "WM", "GM", "OC" = "gmfslOccipital_Lobe", "insula" = "gmfslInsula", "Putamen" = "gmfslPutamen", "Caudate" = "gmfslCaudate", "CBL" = "gmRefCBL", "THA" = "gmfslThalamus", "TC" =
#                "gmfslTemporal_Lobe")

regions <- c("OC" = "gmfslOccipital_Lobe", 
             "CBL" = "gmRefCBL", 
             "FC" = "gmfslFrontal_Lobe")
 
tacs <- tacs %>% 
  select(tacdata, Subjname, PETNo, k2prime_MRTM1,k2prime_MRTM1_weights, logan_tstar) %>% 
  mutate(tacdata = map(tacdata, ~select(.x, all_of(regions), times, 
                                        start, weights, durations))) 


#Long data. By gathering the regions into a single region collumn we can group_by region and then iterate over every region. note, do not select the reference region, as that should be included for every single region and not group_by seperately from the rest!

tacs_long <- tacs %>% 
  unnest(tacdata, .drop = FALSE) %>% 
  gather(key = Region, value = TAC, OC, FC) %>% 
  group_by(Subjname, PETNo, Region) %>% 
  nest(c(CBL, times, weights, logan_tstar, TAC, start,durations, contains("k2prime")),
       .key = "tacdata") %>% 
  mutate(PET = paste(Subjname, PETNo, sep='_')) %>% 
  group_by(Subjname, PETNo, PET, 
                    Region) %>%
  left_join(pmod_t_stars) %>% 
  unnest(tacdata) %>% 
  group_by(Subjname, PETNo, PET, 
                    Region) %>% 
  mutate(frame_number = which.min(abs(mrtm2_t_star-start))) %>% 
  mutate(n = row_number(), max_frame = max(n)) %>% 
  mutate(mrtm2_t_star = (max_frame - frame_number)+1) %>% 
  mutate(frame_number = which.min(abs(ref_logan_t_star-start))) %>% 
  mutate(n = row_number(), max_frame = max(n)) %>% 
  mutate(ref_logan_t_star = (max_frame - frame_number)+1) %>% 
  select(-c(n, frame_number, max_frame)) %>% 
  nest(c(CBL, times, weights, logan_tstar, TAC, start, durations),
       .key = "tacdata") %>% 
   ungroup()
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## All list-columns are now preserved.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

    ## Warning: All elements of `...` must be named.
    ## Did you want `tacdata = c(CBL, times, weights, logan_tstar, TAC, start, durations, k2prime_MRTM1, 
    ##     k2prime_MRTM1_weights)`?

    ## Joining, by = c("Subjname", "PETNo", "Region", "PET")

    ## Warning: All elements of `...` must be named.
    ## Did you want `tacdata = c(CBL, times, weights, logan_tstar, TAC, start, durations)`?

``` r
fitmrtm2_tstar <- function(tacdata, k2prime,
                           t_star) {
  mrtm2(t_tac = tacdata$times, 
      roitac = tacdata$TAC, 
      reftac = tacdata$CBL, 
      k2prime = k2prime,
      tstarIncludedFrames = t_star,
       weights=tacdata$weights, dur = tacdata$durations)
}

fitmrtm2_tstar_weights <- function(tacdata, k2prime,
                           t_star) {
  mrtm2(t_tac = tacdata$times, 
      roitac = tacdata$TAC, 
      reftac = tacdata$CBL, 
      k2prime = k2prime,
      tstarIncludedFrames = t_star, dur = tacdata$durations)
}

fitreflogan_tstar <- function(tacdata, k2prime,
                           t_star) {
  refLogan(t_tac = tacdata$times, 
      roitac = tacdata$TAC, 
      reftac = tacdata$CBL, 
      k2prime = k2prime,
      tstarIncludedFrames = t_star, dur = tacdata$durations)
}
```

``` r
tacs_long <- tacs_long %>% 
    group_by(Subjname, PETNo, Region) %>%
  mutate(fitreflogan_tstar= pmap(list(tacdata, k2prime_MRTM1,
                           ref_logan_t_star),
                           fitreflogan_tstar)) %>%
  mutate(bp_reflogan_tstar = map_dbl(fitreflogan_tstar,
                                  c('par', 'bp'))) %>%
  mutate(MRTM2fit_weights = map2(tacdata, k2prime_MRTM1_weights,
                         ~mrtm2(t_tac=.x$times, 
                                reftac = .x$CBL, 
                               roitac = .x$TAC, 
                               k2prime = .y, dur = .x$durations)),
         bp_MRTM2_weights = map_dbl(MRTM2fit_weights, 
                                    c("par", "bp"))) %>%
  mutate(srtmfit = map(tacdata, ~srtm(t_tac = .x$times, reftac = .x$CBL,      
                                      roitac = .x$TAC))) %>% 
  mutate(bp_srtm_weights = map_dbl(srtmfit, c('par', 'bp')))  %>% 
  mutate(fitmrtm2_tstar= pmap(list(tacdata, k2prime_MRTM1,
                           mrtm2_t_star),
                           fitmrtm2_tstar)) %>%
  mutate(bp_mrtm2_tstar = map_dbl(fitmrtm2_tstar,
                                  c('par', 'bp'))) %>% 
  mutate(fitmrtm2_tstar_weights= pmap(list(tacdata, k2prime_MRTM1_weights,
                           mrtm2_t_star),
                           fitmrtm2_tstar_weights)) %>%
  mutate(bp_mrtm2_tstar_weights = map_dbl(fitmrtm2_tstar_weights,
                                  c('par', 'bp'))) %>% 
  select(-contains("fit"))
```

    ## Registered S3 methods overwritten by 'car':
    ##   method                          from
    ##   influence.merMod                lme4
    ##   cooks.distance.influence.merMod lme4
    ##   dfbeta.influence.merMod         lme4
    ##   dfbetas.influence.merMod        lme4

``` r
tacs_long <- tacs_long %>% 
  select(PET, Subjname ,PETNo, Region, contains("bp"), contains("t_star"))
```

save
====

``` r
saveRDS(tacs_long,'../DerivedData/raw_kinfit_AZ10419369_tstar.rds')

file.copy(from = '../DerivedData/raw_kinfit_AZ10419369_tstar.rds',
          to = '../../PMOD_results/RawData_PMOD_RDS/raw_kinfit_AZ10419369_tstar.rds',
          overwrite = T)
```
