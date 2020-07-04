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
                         ~twotcm(t_tac = .x$Times/60, tac = .x$WB, input = .y,        inpshift.upper = 1, 
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
  rename(weights = Weights) %>% 
  mutate(StartTime = StartTime/60) %>% 
  mutate(Duration = Duration/60) %>%
  nest(c(Times, StartTime, Duration, weights, FC:CBL), .key = 'tacdata')
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

``` r
pmod_t_stars <- readRDS("DerivedData/pmod_macroparameters_t_star.rds") %>% 
  select(subjname:id,t_star, tracer, model, Region = region) %>% 
  filter(tracer == "pbr28") %>% 
  select(-c(tracer,id)) %>% 
  rename(Subjname = subjname) %>% 
  mutate(PETNo = as.numeric(PETNo)) %>% 
  filter(model != "two_tcm") %>% 
  spread(key = model, value = t_star) %>%   rename_at(vars(logan, ma1), ~ paste(., "_t_star", sep = ""))
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
nest(c(Times, StartTime, Duration, weights, TAC),
     .key = 'tacdata') %>% 
  left_join(pmod_t_stars) %>% 
  unnest(tacdata) %>% 
  mutate(frame_number = which.min(abs(ma1_t_star-StartTime))) %>% 
  mutate(n = row_number(), max_frame = max(n)) %>% 
  mutate(ma1_t_star = (max_frame - frame_number)+1) %>% 
  mutate(frame_number = which.min(abs(logan_t_star-StartTime))) %>% 
  mutate(n = row_number(), max_frame = max(n)) %>% 
  mutate(logan_t_star = (max_frame - frame_number)+1) %>% 
  select(PET, Subjname, PETNo, Region, logan_t_star, ma1_t_star, Times, StartTime, Duration, weights, TAC) %>% nest(c(Times, StartTime, Duration, weights, TAC),
     .key = 'tacdata') %>% 
  group_by(PET, Subjname, PETNo) 

tacs_long <- tacs %>%
  select(-tacdata) %>%
  group_by(PET, Subjname, PETNo) %>% 
  inner_join(tacs_long) %>% 
  ungroup()
```

Define functions for fitting the models
=======================================

``` r
# 2TCM using the fitted delay and vB from delayFit
# 
fit2tcm <- function(tacdata, input, delayFit) {
  twotcm(t_tac = tacdata$Times/60, 
         tac = tacdata$TAC, 
         input = input, 
         inpshift = delayFit$par$inpshift, vB= 0.05, frameStartEnd = c(1,33),
         multstart_iter = 100)
}


# MA1 using the fitted delay and vB from delayFit
fitma1_tstar <- function(tacdata, input, delayFit, ma1_t_star) {
  ma1(t_tac = tacdata$Times/60, 
      tac = tacdata$TAC, input = input, tstarIncludedFrames = ma1_t_star,
      inpshift = delayFit$par$inpshift, weights=tacdata$weights, 
      vB= 0, frameStartEnd = c(1,33), dur = tacdata$Duration)
}

#Loganplot
fit_Logan_tstar <- function(tacdata, input, delayFit, logan_t_star) {
  Loganplot(t_tac = tacdata$Times/60,
            tac = tacdata$TAC,
         input = input, 
         inpshift = delayFit$par$inpshift, tstarIncludedFrames = logan_t_star, 
         vB= 0, frameStartEnd = c(1,33),
         dur = tacdata$Duration)
}

fitma1_weights <- function(tacdata, input, delayFit) {
  ma1(t_tac = tacdata$Times/60, 
      tac = tacdata$TAC, input = input, tstarIncludedFrames = 6,
      inpshift = delayFit$par$inpshift, vB= 0, frameStartEnd = c(1,33),
      dur = tacdata$Duration)
}

fitma1_tstar_weights <- function(tacdata, input, delayFit, ma1_t_star) {
  ma1(t_tac = tacdata$Times/60,
      tac = tacdata$TAC, input = input, tstarIncludedFrames = ma1_t_star,
      inpshift = delayFit$par$inpshift, vB= 0, frameStartEnd = c(1,33), 
      dur = tacdata$Duration)
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
saveRDS(tacs_long,'DerivedData/raw_kinfit_pbr28_tstar.rds')

file.copy(from = '../DerivedData/raw_kinfit_pbr28_tstar.rds',
          to = '../PMOD_results/RawData_PMOD_RDS/raw_kinfit_pbr28_tstar.rds',
          overwrite = T)
```
