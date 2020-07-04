-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)

Aims

The aim of this code is to convert the data from the kinfitr analysis to
the same format as the data from the PMOD analysis is in, such that all
of this data may be combined into a single table.

Libraries
=========

CRAN libraries
--------------

package installation

``` r
install.packages("tidyverse")
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
install.packages("viridis")
install.packages("janitor")
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
library(viridis)
library(janitor)
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
# library(vmisc)
library(kipettools)
library(granviller)
library(relfeas)

theme_set(theme_light())
```

\#load pbr28 data

``` r
pbr28 <- readRDS('RawData_PMOD_RDS/raw_kinfit_pbr28_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region = Region, contains("Vt"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save pbr28 data just with the macroparameter Vt in a similar format to
the PMOD

NOTE: COULD NOT FIND STANDARD ERROR

``` r
pbr28 <-  pbr28 %>%
  mutate(tracer = replicate(length(PET), "pbr28")) %>% 
  mutate(Measure = replicate(length(PET), "Vt")) %>% 
  gather(model, Value, contains("Vt")) %>% 
  mutate(model = str_replace(string= model, pattern="^.*?_", replacement = ""))
```

\#load pk11195 data

``` r
pk11195 <- readRDS('RawData_PMOD_RDS/raw_kinfit_pk11195_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, contains("Vt"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save pk11195 data just with the macroparameter Vt in a similar format to
the PMOD

``` r
pk11195 <-  pk11195 %>%
  mutate(tracer = replicate(length(PET), "pk")) %>% 
  mutate(Measure = replicate(length(PET), "Vt")) %>% 
  gather(model, Value, contains("Vt")) %>% 
  mutate(model = str_replace(string= model, pattern="^.*?_", replacement = ""))
```

\#load SCH data

``` r
sch <- readRDS('RawData_PMOD_RDS/raw_kinfit_SCH23390_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, contains("bp"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save SCH23390 data just with the macroparameter Vt in a similar format
to the PMOD

``` r
sch <-  sch %>%
  mutate(tracer = replicate(length(PET), "sch")) %>% 
  mutate(Measure = replicate(length(PET), "bp")) %>% 
  gather(model, Value, contains("bp")) %>% 
  mutate(model = str_replace(string= model, pattern="^.*?_", replacement = ""))
```

\#load AZ data

``` r
az <- readRDS('RawData_PMOD_RDS/raw_kinfit_AZ10419369_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, contains("bp"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save AZ10419369 data just with the macroparameter Vt in a similar format
to the PMOD

``` r
az <-  az %>%
  mutate(tracer = replicate(length(PET), "az")) %>% 
  mutate(Measure = replicate(length(PET), "bp")) %>% 
  gather(model, Value, contains("bp")) %>% 
  mutate(model = str_replace(string= model, pattern="^.*?_", replacement = ""))
```

\#bind all four tracers into a single tibble

``` r
kinfitr_macro <- bind_rows(az, sch, pbr28, pk11195) %>%
  mutate(software = replicate(length(PET), "kinfitR")) %>% 
  rename(measure = Measure)
```

\#\#Prepare kinfitR macroparameter dataset

As different regions are required for different tracers, I will create 4
different datasets, one for each tracer, and then filter by the regions
I want and then combine the datasets once again. OC and FC for “9369”
(az), THA and FC for pbr28 and pk11195 as both are TSPO-ligands. STR and
FC for SCH.

``` r
kinfitr_az <- kinfitr_macro %>% 
  filter(tracer == "az") %>% 
  filter(region %in% c( "OC", "FC"))

kinfitr_sch <- kinfitr_macro %>% 
  filter(tracer == "sch") %>% 
  filter(region %in% c( "STR", "FC"))

kinfitr_pk11195 <- kinfitr_macro %>% 
  filter(tracer == "pk") %>% 
  filter(region %in% c( "THA", "FC"))

kinfitr_pbr28 <- kinfitr_macro %>% 
  filter(tracer == "pbr28") %>% 
  filter(region %in% c( "THA", "FC"))

kinfitr_macro <- bind_rows(kinfitr_az, kinfitr_sch, kinfitr_pbr28, kinfitr_pk11195) %>% 
  ungroup()

# saveRDS(kinfitr_macro, 'Data_PMOD_RDS/kinfitr_macroparameters_tstar.rds')
```

\#Microparameters

\#\#Preparing microparameters for ligands requiring blood sampling

``` r
pbr28_micro <- readRDS('RawData_PMOD_RDS/raw_kinfit_pbr28.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region = Region, fit_2tcm) %>% 
  mutate(K1 = map_dbl(fit_2tcm, c('par', 'K1'))) %>% 
  mutate(k2 = map_dbl(fit_2tcm, c('par', 'k2'))) %>%
  mutate(k3 = map_dbl(fit_2tcm, c('par', 'k3'))) %>%
  mutate(k4 = map_dbl(fit_2tcm, c('par', 'k4'))) %>% 
  select(-fit_2tcm) %>% 
  mutate(Ligand = replicate(length(PET), "PBR28")) %>% 
  filter(region == "FC") %>% 
  group_by(subjname) %>% 
  nest() %>% 
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()

pk11195_micro <- readRDS('RawData_PMOD_RDS/raw_kinfit_pk11195.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region = Region,  fit_2tcm) %>% 
  mutate(K1 = map_dbl(fit_2tcm, c('par', 'K1'))) %>% 
  mutate(k2 = map_dbl(fit_2tcm, c('par', 'k2'))) %>%
  mutate(k3 = map_dbl(fit_2tcm, c('par', 'k3'))) %>%
  mutate(k4 = map_dbl(fit_2tcm, c('par', 'k4'))) %>% 
  select(-fit_2tcm) %>% 
  mutate(Ligand = replicate(length(PET), "PK11195")) %>% 
  filter(region == "FC") %>% 
  group_by(subjname) %>% 
  nest() %>% 
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()

micro_aif_kinfitr <- bind_rows(pbr28_micro, pk11195_micro) %>% 
  mutate(software = replicate(length(PET), "kinfitR")) %>% 
    mutate(model = replicate(length(PET), "two_tcm")) 

# saveRDS(micro_aif_kinfitr, 'Data_PMOD_RDS/micro_aif_kinfitr.rds' )
```

\#\#loading the raw data for ligands that use reference tissue models

``` r
sch_micro <- readRDS('RawData_PMOD_RDS/raw_kinfit_SCH23390.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region = Region, srtmfit) %>% 
  mutate(par = map(srtmfit, 'par')) %>% 
  unnest(par, .drop = TRUE) %>% 
  select(- bp) %>% 
  mutate(Ligand = replicate(length(PET), "SCH23390")) %>% 
  filter(region == "STR")%>% 
  group_by(subjname) %>% 
  nest() %>% 
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()

az_micro <- readRDS('RawData_PMOD_RDS/raw_kinfit_AZ10419369.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region = Region, srtmfit) %>% 
  mutate(par = map(srtmfit, 'par')) %>% 
  unnest(par, .drop = TRUE) %>% 
  select(- bp) %>%
  mutate(Ligand = replicate(length(PET), "AZ10419369")) %>% 
  filter(region == "OC")%>% 
  group_by(subjname) %>% 
  nest() %>% 
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()

micro_ref_kinfitr <- bind_rows(az_micro, sch_micro) %>% 
  mutate(software = replicate(length(PET), "kinfitR")) %>% 
    mutate(model = replicate(length(PET), "SRTM")) 

# saveRDS(micro_ref_kinfitr, 'Data_PMOD_RDS/micro_ref_kinfitr.rds' )
```
