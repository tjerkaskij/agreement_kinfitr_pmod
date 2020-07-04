-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)

Aims

The aim of this code is to convert the data from the kinfitr analysis
which was run without iteration for the radioligands PK11195 and PBR28
to the same format as the data from the PMOD analysis is in, such that
all of this data may be combined into a single table.

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
pbr28 <- readRDS('RawData_PMOD_RDS/raw_kinfit_pbr28_no_iter.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region = Region, Vt_two_tcm = Vt_2tcm, Vt_ma1, Vt_logan = Vt_Logan)
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save pbr28 data just with the macroparameter Vt in a similar format to
the PMOD

NOTE: COULD NOT FIND STANDARD ERROR

``` r
pbr28 <-  pbr28 %>%
  mutate(tracer = replicate(length(PET), "pbr28")) %>% 
  mutate(Measure = replicate(length(PET), "Vt")) %>% 
  gather(model, Value, Vt_two_tcm:Vt_logan) %>% 
  mutate(model = str_replace(string= model, pattern="^.*?_", replacement = ""))

# saveRDS(pbr28, 'Data_PMOD_RDS/macro_kinfit_pbr28_no_iter.rds' )
```

\#load pk11195 data

``` r
pk11195 <- readRDS('RawData_PMOD_RDS/raw_kinfit_pk11195_no_iter.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, Vt_two_tcm = Vt_2tcm, Vt_ma1, Vt_logan = Vt_Logan)
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save pk11195 data just with the macroparameter Vt in a similar format to
the PMOD

``` r
pk11195 <-  pk11195 %>%
  mutate(tracer = replicate(length(PET), "pk")) %>% 
  mutate(Measure = replicate(length(PET), "Vt")) %>% 
  gather(model, Value, Vt_two_tcm:Vt_logan) %>% 
  mutate(model = str_replace(string= model, pattern="^.*?_", replacement = ""))

# saveRDS(pk11195, 'Data_PMOD_RDS/macro_kinfit_pk11195_no_iter.rds' )
```

\#bind all four tracers into a single tibble

``` r
kinfitr_macro <- bind_rows(pbr28, pk11195) %>%
  mutate(software = replicate(length(PET), "no_iter")) %>% 
  rename(measure = Measure)

macro <- kinfitr_macro %>% 
    filter(region %in% c("FC", "TC", "THA", "WB", "STR", "OC"))

# saveRDS(macro, 'Data_PMOD_RDS/kinfitr_macroparameters_more.rds')
```

\#\#Prepare kinfitR macroparameter dataset

As different regions are required for different tracers, I will create 4
different datasets, one for each tracer, and then filter by the regions
I want and then combine the datasets once again. OC and FC for “9369”
(az), THA and FC for pbr28 and pk11195 as both are TSPO-ligands. STR and
FC for SCH.

``` r
kinfitr_pk11195 <- kinfitr_macro %>% 
  filter(tracer == "pk") %>% 
  filter(region %in% c( "THA", "FC"))

kinfitr_pbr28 <- kinfitr_macro %>% 
  filter(tracer == "pbr28") %>% 
  filter(region %in% c( "THA", "FC"))

kinfitr_macro <- bind_rows(kinfitr_pbr28, kinfitr_pk11195) %>% 
  ungroup()

 # saveRDS(kinfitr_macro, 'Data_PMOD_RDS/kinfitr_macroparameters_no_iter.rds')
```
