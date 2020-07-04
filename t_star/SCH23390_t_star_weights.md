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

The aim of this assignment is to analyze the SCH23390 data in kinfitr
using the t\* values fitted by PMOD and constant weights

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
  nest(-c(Subjname,PETNo), .key = 'tacdata')
```

    ## Warning: All elements of `...` must be named.
    ## Did you want `tacdata = c(L_WM, R_WM, WM, GM, WB, L_MB, R_MB, MB, gmL_AST_gm, gmR_AST_gm, 
    ##     gmAST_gm, gmL_STR_gm, gmR_STR_gm, gmSTR_gm, gmVST_m_gm, gmpostPU_m_gm, 
    ##     L_AST, R_AST, AST, L_STR, R_STR, STR, VST_m, postPU_m, gmL_CBL, 
    ##     gmR_CBL, gmCBL, gmL_VST, gmR_VST, gmVST, gmL_postCA, gmR_postCA, 
    ##     gmpostCA, gmL_postPU, gmR_postPU, gmpostPU, gmL_preCA, gmR_preCA, 
    ##     gmpreCA, gmL_prePU, gmR_prePU, gmprePU, L_abc, R_abc, abc, 
    ##     help2, help3, L_FSLSFC, R_FSLSFC, FSLSFC, L_FSLSOC, R_FSLSOC, 
    ##     FSLSOC, L_FSLSPC, R_FSLSPC, FSLSPC, L_FSLSTC, R_FSLSTC, FSLSTC, 
    ##     L_FSLSLL, R_FSLSLL, FSLSLL, L_FSLSSTR, R_FSLSSTR, FSLSSTR, 
    ##     L_FSLSACC, R_FSLSACC, FSLSACC, L_FSLSAMG, R_FSLSAMG, FSLSAMG, 
    ##     FSLSBS, L_FSLSCAU, R_FSLSCAU, FSLSCAU, L_FSLSCER, R_FSLSCER, 
    ##     FSLSCER, L_FSLSCERCX, R_FSLSCERCX, FSLSCERCX, L_FSLSHIP, 
    ##     R_FSLSHIP, FSLSHIP, L_FSLSINS, R_FSLSINS, FSLSINS, L_FSLSLFC, 
    ##     R_FSLSLFC, FSLSLFC, L_FSLSLOC, R_FSLSLOC, FSLSLOC, L_FSLSLPC, 
    ##     R_FSLSLPC, FSLSLPC, L_FSLSLTC, R_FSLSLTC, FSLSLTC, L_FSLSMFC, 
    ##     R_FSLSMFC, FSLSMFC, L_FSLSMOC, R_FSLSMOC, FSLSMOC, L_FSLSMPC, 
    ##     R_FSLSMPC, FSLSMPC, L_FSLSMTC, R_FSLSMTC, FSLSMTC, L_FSLSOFC, 
    ##     R_FSLSOFC, FSLSOFC, L_FSLSPAL, R_FSLSPAL, FSLSPAL, L_FSLSPCC, 
    ##     R_FSLSPCC, FSLSPCC, L_FSLSPHIP, R_FSLSPHIP, FSLSPHIP, L_FSLSPUT, 
    ##     R_FSLSPUT, FSLSPUT, L_FSLSSMC, R_FSLSSMC, FSLSSMC, L_FSLSSN, 
    ##     R_FSLSSN, FSLSSN, L_FSLSTHA, R_FSLSTHA, FSLSTHA, L_FSLSVDC, 
    ##     R_FSLSVDC, FSLSVDC, `L_fsldctx-parstriangularis`, `R_fsldctx-parstriangularis`, 
    ##     `fsldctx-parstriangularis`, `L_fsldctx-rostralmiddlefrontal`, 
    ##     `R_fsldctx-rostralmiddlefrontal`, `fsldctx-rostralmiddlefrontal`, 
    ##     `L_fsldctx-caudalmiddlefrontal`, `R_fsldctx-caudalmiddlefrontal`, 
    ##     `fsldctx-caudalmiddlefrontal`, `L_fsldctx-frontalpole`, `R_fsldctx-frontalpole`, 
    ##     `fsldctx-frontalpole`, `L_fsldctx-medialorbitofrontal`, `R_fsldctx-medialorbitofrontal`, 
    ##     `fsldctx-medialorbitofrontal`, `L_fsldctx-superiorfrontal`, 
    ##     `R_fsldctx-superiorfrontal`, `fsldctx-superiorfrontal`, `L_fsldctx-lateralorbitofrontal`, 
    ##     `R_fsldctx-lateralorbitofrontal`, `fsldctx-lateralorbitofrontal`, 
    ##     `L_fsldctx-parsorbitalis`, `R_fsldctx-parsorbitalis`, `fsldctx-parsorbitalis`, 
    ##     `L_fsldctx-lateraloccipital`, `R_fsldctx-lateraloccipital`, 
    ##     `fsldctx-lateraloccipital`, `L_fsldctx-cuneus`, `R_fsldctx-cuneus`, 
    ##     `fsldctx-cuneus`, `L_fsldctx-lingual`, `R_fsldctx-lingual`, 
    ##     `fsldctx-lingual`, `L_fsldctx-pericalcarine`, `R_fsldctx-pericalcarine`, 
    ##     `fsldctx-pericalcarine`, `L_fsldctx-inferiorparietal`, `R_fsldctx-inferiorparietal`, 
    ##     `fsldctx-inferiorparietal`, `L_fsldctx-superiorparietal`, 
    ##     `R_fsldctx-superiorparietal`, `fsldctx-superiorparietal`, 
    ##     `L_fsldctx-supramarginal`, `R_fsldctx-supramarginal`, `fsldctx-supramarginal`, 
    ##     `L_fsldctx-precuneus`, `R_fsldctx-precuneus`, `fsldctx-precuneus`, 
    ##     `L_fsldctx-inferiortemporal`, `R_fsldctx-inferiortemporal`, 
    ##     `fsldctx-inferiortemporal`, `L_fsldctx-middletemporal`, `R_fsldctx-middletemporal`, 
    ##     `fsldctx-middletemporal`, `L_fsldctx-superiortemporal`, `R_fsldctx-superiortemporal`, 
    ##     `fsldctx-superiortemporal`, `L_fsldctx-transversetemporal`, 
    ##     `R_fsldctx-transversetemporal`, `fsldctx-transversetemporal`, 
    ##     `L_fsldctx-entorhinal`, `R_fsldctx-entorhinal`, `fsldctx-entorhinal`, 
    ##     `L_fsldctx-fusiform`, `R_fsldctx-fusiform`, `fsldctx-fusiform`, 
    ##     `L_fsldctx-temporalpole`, `R_fsldctx-temporalpole`, `fsldctx-temporalpole`, 
    ##     L_fsldAmygdala, R_fsldAmygdala, fsldAmygdala, `L_fsldAmygdala-Anterior`, 
    ##     `R_fsldAmygdala-Anterior`, `fsldAmygdala-Anterior`, L_fsldHippocampus, 
    ##     R_fsldHippocampus, fsldHippocampus, `L_fsldctx-parahippocampal`, 
    ##     `R_fsldctx-parahippocampal`, `fsldctx-parahippocampal`, `L_fsldctx-caudalanteriorcingulate`, 
    ##     `R_fsldctx-caudalanteriorcingulate`, `fsldctx-caudalanteriorcingulate`, 
    ##     `L_fsldctx-rostralanteriorcingulate`, `R_fsldctx-rostralanteriorcingulate`, 
    ##     `fsldctx-rostralanteriorcingulate`, `L_fsldctx-isthmuscingulate`, 
    ##     `R_fsldctx-isthmuscingulate`, `fsldctx-isthmuscingulate`, 
    ##     `L_fsldctx-posteriorcingulate`, `R_fsldctx-posteriorcingulate`, 
    ##     `fsldctx-posteriorcingulate`, L_fsldCaudate, R_fsldCaudate, 
    ##     fsldCaudate, `L_fsldAccumbens-area`, `R_fsldAccumbens-area`, 
    ##     `fsldAccumbens-area`, L_fsldPutamen, R_fsldPutamen, fsldPutamen, 
    ##     `L_fsldCerebellum-Cortex`, `R_fsldCerebellum-Cortex`, `fsldCerebellum-Cortex`, 
    ##     L_fsldInsula, R_fsldInsula, fsldInsula, `L_fsldctx-insula`, 
    ##     `R_fsldctx-insula`, `fsldctx-insula`, L_fsldPallidum, R_fsldPallidum, 
    ##     fsldPallidum, L_fsldOperculum, R_fsldOperculum, fsldOperculum, 
    ##     `L_fsldctx-paracentral`, `R_fsldctx-paracentral`, `fsldctx-paracentral`, 
    ##     `L_fsldctx-parsopercularis`, `R_fsldctx-parsopercularis`, 
    ##     `fsldctx-parsopercularis`, `L_fsldctx-postcentral`, `R_fsldctx-postcentral`, 
    ##     `fsldctx-postcentral`, `L_fsldctx-precentral`, `R_fsldctx-precentral`, 
    ##     `fsldctx-precentral`, L_fsldThalamus, R_fsldThalamus, fsldThalamus, 
    ##     `L_fsldThalamus-Proper`, `R_fsldThalamus-Proper`, `fsldThalamus-Proper`, 
    ##     L_fsldVentralDC, R_fsldVentralDC, fsldVentralDC, L_FSPVMPFC, 
    ##     R_FSPVMPFC, FSPVMPFC, L_FSPLPFC, R_FSPLPFC, FSPLPFC, L_FSPDLPFC, 
    ##     R_FSPDLPFC, FSPDLPFC, L_FSPMFC, R_FSPMFC, FSPMFC, L_FSPMOFC, 
    ##     R_FSPMOFC, FSPMOFC, L_FSPOFC, R_FSPOFC, FSPOFC, L_FSLSLFC2, 
    ##     R_FSLSLFC2, FSLSLFC2, L_FSG1ACORT, R_FSG1ACORT, FSG1ACORT, 
    ##     L_FSG1ASC, R_FSG1ASC, FSG1ASC, L_FSGTSC, R_FSGTSC, FSGTSC, 
    ##     L_FSGTCORT, R_FSGTCORT, FSGTCORT, L_FSGD1EXSTR, R_FSGD1EXSTR, 
    ##     FSGD1EXSTR, L_FSGD1CORT, R_FSGD1CORT, FSGD1CORT, times, durations, 
    ##     L_abc1, gmR_abc, help1, start, end, weights)`?

Fitting and plotting MRTM1 multiple times for regions FC and WB
===============================================================

``` r
#K2 prime for striatum to be used to MRTM2 later on

 tacs <- tacs %>% 
  group_by(Subjname, PETNo) %>% 
  mutate(MRTM1fit = map(tacdata, ~mrtm1(t_tac = .x$times, reftac = .x$gmCBL,
                                   roitac = .x$FSLSSTR, 
                                      weights = .x$weights, 
                                   dur = .x$durations))) %>% 
  mutate(MRTM1fit_weights = map(tacdata, ~mrtm1(t_tac = .x$times, reftac = .x$gmCBL,
                                   roitac = .x$FSLSSTR, dur = .x$durations))) %>% 
  ungroup() %>% 
  mutate(k2prime_MRTM1 = map_dbl(MRTM1fit, c("par", "k2prime"))) %>% 
   mutate(k2prime_MRTM1_weights = map_dbl(MRTM1fit_weights, c("par", "k2prime"))) %>% 
  group_by(Subjname, PETNo) %>% 
  mutate(logan_tstar = map2(tacdata, k2prime_MRTM1,
                 ~refLogan_tstar(t_tac = .x$times, 
                    reftac = .x$gmCBL, lowroi = .x$FSLSFC, 
                    medroi = .x$FSLSINS, highroi = .x$FSLSSTR, 
                   k2prime = .y)))
```

t\* values that were fitted by PMOD, in minutes
===============================================

``` r
pmod_t_stars <- readRDS("../DerivedData/pmod_macroparameters_t_star.rds") %>% 
  select(subjname:id,t_star, tracer, model, Region = region) %>% 
  filter(tracer == "sch") %>% 
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
# regions <- c("STR" = "FSLSSTR", "FC" = "FSLSFC", "WB", "WM", "GM", "OC" = "FSLSOC", "insula" = "FSLSINS", "Putamen" = "FSLSPUT", "ACC" = "FSLSACC", "CBL" = "gmCBL",
#              "THA" = "FSLSTHA", "TC" = "FSLSTC")

regions <- c("STR" = "FSLSSTR", "FC" = "FSLSFC", "CBL" = "gmCBL")     
 
tacs <- tacs %>% 
  select(tacdata, Subjname, PETNo, k2prime_MRTM1, logan_tstar, k2prime_MRTM1_weights) %>% 
  mutate(tacdata = map(tacdata, ~select(.x, all_of(regions), times, weights, start, durations))) 


#Long data. By gathering the regions into a single region collumn we can group_by region and then iterate over every region. note, do not select the reference region, as that should be included for every single region and not group_by seperately from the rest!


tacs_long <- tacs %>% 
  unnest(tacdata, .drop = FALSE) %>% 
  gather(key = Region, value = TAC, STR, FC) %>% 
  group_by(Subjname, PETNo, Region) %>% 
  nest(c(CBL, times, weights, logan_tstar, TAC, start, durations, contains("k2prime")),
       .key = "tacdata") %>% 
  mutate(PET = paste(Subjname, PETNo, sep='_')) %>% 
  group_by(PET, Subjname, PETNo,  
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

functions

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
saveRDS(tacs_long,'../DerivedData/raw_kinfit_SCH23390_tstar.rds')

file.copy(from = '../DerivedData/raw_kinfit_SCH23390_tstar.rds',
          to = '../../PMOD_results/RawData_PMOD_RDS/raw_kinfit_SCH23390_tstar.rds',
          overwrite = T)
```

    ## [1] TRUE
