-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
-   [invasive models and reference tissue
    models](#invasive-models-and-reference-tissue-models)
-   [median t\* for the invasive
    models](#median-t-for-the-invasive-models)

Aims

The purpose of this code is to calculate the median time for t\* fitted
by PMOD.

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
library(kableExtra)
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
  select(PET, subjname = Subjname, PETNo, region = Region, contains("t_star"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save pbr28 data just with the macroparameter Vt in a similar format to
the PMOD

NOTE: COULD NOT FIND STANDARD ERROR

``` r
pbr28 <-  pbr28 %>%
  mutate(tracer = replicate(length(PET), "pbr28")) %>% 
  mutate(Measure = replicate(length(PET), "Vt")) %>% 
  gather(model, Value, contains("t_star")) %>% 
  mutate(model = str_replace(string= model, pattern="_.*?$", replacement = ""))
```

\#load pk11195 data

``` r
pk11195 <- readRDS('RawData_PMOD_RDS/raw_kinfit_pk11195_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, contains("t_star"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save pk11195 data just with the macroparameter Vt in a similar format to
the PMOD

``` r
pk11195 <-  pk11195 %>%
  mutate(tracer = replicate(length(PET), "pk")) %>% 
  mutate(Measure = replicate(length(PET), "Vt")) %>% 
  gather(model, Value, contains("t_star")) %>% 
  mutate(model = str_replace(string= model, pattern="_.*?$", replacement = ""))
```

\#load SCH data

``` r
sch <- readRDS('RawData_PMOD_RDS/raw_kinfit_SCH23390_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, contains("t_star"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save SCH23390 data just with the macroparameter Vt in a similar format
to the PMOD

``` r
sch <-  sch %>%
  mutate(tracer = replicate(length(PET), "sch")) %>% 
  mutate(Measure = replicate(length(PET), "bp")) %>% 
  gather(model, Value, contains("t_star")) %>% 
  mutate(model = str_replace(string= model, pattern="_.*?$", replacement = ""))%>% 
  mutate(model = str_replace(string = model, pattern = "ref", replacement = "ref_logan"))
```

\#load AZ data

``` r
az <- readRDS('RawData_PMOD_RDS/raw_kinfit_AZ10419369_tstar.rds') %>% 
  ungroup() %>% 
  select(PET, subjname = Subjname, PETNo, region =  Region, contains("t_star"))
```

\#\#adding tracer and model columns. Gather Vt’s into a single column

Save AZ10419369 data just with the macroparameter Vt in a similar format
to the PMOD

``` r
az <-  az %>%
  mutate(tracer = replicate(length(PET), "az")) %>% 
  mutate(Measure = replicate(length(PET), "bp")) %>% 
 gather(model, Value, contains("t_star")) %>% 
  mutate(model = str_replace(string= model, pattern="_.*?$", replacement = "")) %>% 
  mutate(model = str_replace(string = model, pattern = "ref", replacement = "ref_logan"))
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
```

invasive models and reference tissue models
===========================================

``` r
macro_aif <- kinfitr_macro %>% 
  filter(measure == "Vt")

macro_ref <- kinfitr_macro %>% 
  filter(measure == "bp")
```

median t\* for the invasive models
==================================

``` r
macro_aif <- macro_aif %>% 
  spread(key = model, value = Value) %>% 
  select( region, tracer, logan , ma1) %>% 
  rename(Region = region, Ligand = tracer, Logan = "logan", MA1 = "ma1") %>% 
  group_by(Region, Ligand) %>% 
  summarize(median_logan = median(Logan), iqr_logan = IQR(Logan), min_logan = min(Logan), max_logan = max(Logan), median_ma1 = median(MA1), iqr_ma1 = IQR(MA1),
            min_ma1 = min(MA1),max_ma1 = max(MA1)) %>% 
  select(Ligand, Region, everything())

macro_aif <- kable(macro_aif, format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "t* values fitted by PMOD for the invasive models",
      col.names = c("Ligand", "Region", "Median", "IQR", "Min",
                    "Max", "Median", "IQR", "Min",
                    "Max"), align= 'c') %>%
kable_styling(full_width = F) %>%
add_header_above(c(" " = 2, "Logan" = 4, "MA1" = 4))

# save_kable(macro_aif, "t_star_medians_aif.html")
```

``` r
macro_ref <- macro_ref %>% 
  spread(key = model, value = Value) %>% 
  select(Region = region, Ligand = tracer, ref_logan , mrtm2) %>% 
  group_by(Region, Ligand) %>% 
  summarize(median_logan = median(ref_logan), iqr_logan = IQR(ref_logan), min_logan = min(ref_logan), max_logan = max(ref_logan), median_ma1 = median(mrtm2), iqr_ma1 = IQR(mrtm2),
            min_ma1 = min(mrtm2),max_ma1 = max(mrtm2)) %>% 
  select(Ligand, Region, everything())

macro_aif <- kable(macro_ref, format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "t* values fitted by PMOD for the reference tissue models",
      col.names = c("Ligand", "Region", "Median", "IQR", "Min",
                    "Max", "Median", "IQR", "Min",
                    "Max"), align= 'c') %>%
kable_styling(full_width = F) %>%
add_header_above(c(" " = 2, "ref Logan" = 4, "MRTM2" = 4))

# save_kable(macro_aif, "t_star_medians_ref.html")
```
