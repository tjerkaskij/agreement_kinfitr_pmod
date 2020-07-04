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
-   [Define functions for fitting the
    models](#define-functions-for-fitting-the-models)
-   [Fit kinetic models](#fit-kinetic-models)
-   [filtering unnecessary variables](#filtering-unnecessary-variables)
-   [save](#save)

Aims
====

The aim of this assignment is to analyze the PK11195 data in kinfitr
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
library(ggplotify)
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
nest(.key = "tacdata") %>% 
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
  # mutate(Subjname = hash_names(Subjname, full = F)) %>% 
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

``` r
pmod_t_stars <- readRDS("../DerivedData/pmod_macroparameters_t_star.rds") %>% 
  select(subjname:id,t_star, tracer, model, Region = region) %>% 
  filter(tracer == "pk") %>% 
  select(-c(tracer,id)) %>% 
  rename(Subjname = subjname) %>% 
  mutate(PETNo = as.numeric(PETNo)) %>% 
  filter(model != "two_tcm") %>% 
  spread(key = model, value = t_star) %>% rename_at(vars(logan, ma1), ~ paste(., "_t_star", sep = ""))
```

\#Rearrangement of the Data into Long Format

``` r
tacs_long <-  tacs %>%
  group_by(PET, Subjname, PETNo) %>% 
  unnest(tacdata, .drop = FALSE) %>% 
gather(Region, TAC, FC:CBL) %>%
   filter(Region %in% c( "THA", "FC")) %>% 
group_by(PET, Subjname, PETNo, 
         Region) %>%
nest(c(Times, start, durations, end, weights, TAC),
     .key = 'tacdata') %>% 
  left_join(pmod_t_stars) %>% 
  unnest(tacdata) %>% 
  mutate(frame_number = which.min(abs(ma1_t_star-start))) %>% 
  mutate(n = row_number(), max_frame = max(n)) %>% 
  mutate(ma1_t_star = (max_frame - frame_number)+1) %>% 
  mutate(frame_number = which.min(abs(logan_t_star-start))) %>% 
  mutate(n = row_number(), max_frame = max(n)) %>% 
  mutate(logan_t_star = (max_frame - frame_number)+1) %>% 
  select(PET, Subjname, PETNo, Region, logan_t_star, ma1_t_star, Times, start, durations, end, weights, TAC) %>% nest(c(Times, start, durations, end, weights, TAC),
     .key = 'tacdata') 

tacs_long <- tacs %>%
  select(-tacdata) %>%
  inner_join(tacs_long, by = c("PET", "Subjname", "PETNo"))
```

Define functions for fitting the models
=======================================

``` r
#2TCM using the fitted delay and vB from delayFit

# fit2tcm <- function(tacdata, input, delayFit) {
#   twotcm(t_tac = tacdata$Times,
#          tac = tacdata$TAC,
#          input = input,
#          inpshift = delayFit$par$inpshift, vB= 0.05, multstart_iter = 100)
# }

fit2tcm <- function(tacdata, input, delayFit) {
  twotcm(t_tac = tacdata$Times,
         tac = tacdata$TAC,
         input = input,
         inpshift = delayFit$par$inpshift, vB= 0.05)
}


# MA1 using the fitted delay and vB from delayFit
fitma1_tstar <- function(tacdata, input, delayFit, ma1_t_star) {
  ma1(t_tac = tacdata$Times, 
      tac = tacdata$TAC, input = input, tstarIncludedFrames = ma1_t_star,
      inpshift = delayFit$par$inpshift, weights=tacdata$weights, vB= 0,
      dur = tacdata$durations)
}

#Loganplot
fit_Logan_tstar <- function(tacdata, input, delayFit, logan_t_star) {
  Loganplot(t_tac = tacdata$Times,
            tac = tacdata$TAC,
         input = input, 
         inpshift = delayFit$par$inpshift, tstarIncludedFrames = logan_t_star, 
         vB= 0, dur = tacdata$durations)
}

fitma1_weights <- function(tacdata, input, delayFit) {
  ma1(t_tac = tacdata$Times, 
      tac = tacdata$TAC, input = input, tstarIncludedFrames = 6,
      inpshift = delayFit$par$inpshift, vB= 0, dur = tacdata$durations)
}

fitma1_tstar_weights <- function(tacdata, input, delayFit, ma1_t_star) {
  ma1(t_tac = tacdata$Times,
      tac = tacdata$TAC, input = input, tstarIncludedFrames = ma1_t_star,
      inpshift = delayFit$par$inpshift, vB= 0, dur = tacdata$durations)
}
```

Fit kinetic models
==================

``` r
set.seed(123)

tacs_long <- tacs_long %>% 

  # 2TCM using fitted vB and delay
  mutate(fit_2tcm_weight= pmap(list(tacdata, input, delayFit), fit2tcm)) %>%
  mutate(Vt_2tcm_weight = map_dbl(fit_2tcm_weight,
      c('par', 'Vt'))) %>%
  
  # MA1
  mutate(fit_ma1_tstar = pmap(list(tacdata, input, delayFit, ma1_t_star), fitma1_tstar))  %>%
  mutate(Vt_ma1_tstar = map_dbl(fit_ma1_tstar, c('par', 'Vt'))) %>%
  
# MA1
  mutate(fit_ma1_weight = pmap(list(tacdata, input, delayFit), fitma1_weights))  %>%
  mutate(Vt_ma1_weight = map_dbl(fit_ma1_weight, c('par', 'Vt'))) %>%
  
# MA1
  mutate(fit_ma1_tstar_weight = pmap(list(tacdata, input, delayFit, ma1_t_star), fitma1_tstar_weights))  %>%
  mutate(Vt_ma1_tstar_weight = map_dbl(fit_ma1_tstar_weight, c('par', 'Vt'))) %>%
   
  #Loganplot
  mutate(Loganfit_tstar = pmap(list(tacdata, input, delayFit, logan_t_star), fit_Logan_tstar)) %>% 
  mutate(Vt_Logan_tstar = map_dbl(Loganfit_tstar, 
          c("par", "Vt")))
```

filtering unnecessary variables
===============================

``` r
tacs_long <- tacs_long %>% 
  select(PET, Subjname, PETNo, Region, contains("Vt"), contains("t_star"))
```

save
====

``` r
saveRDS(tacs_long,'../DerivedData/raw_kinfit_pk11195_tstar.rds')

file.copy(from = '../DerivedData/raw_kinfit_pk11195_tstar.rds',
          to = '../../PMOD_results/RawData_PMOD_RDS/raw_kinfit_pk11195_tstar.rds',
          overwrite = T)
```

    ## [1] TRUE
