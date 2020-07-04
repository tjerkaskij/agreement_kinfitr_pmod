-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)

Aims
====

The aim of this assignment is to convert the data from the PMOD analysis
into a format which can be combined with the data from the kinfitr
analysis in the form of a table.

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
install.packages("remotes")
remotes::install_github("mvuorre/vmisc")
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

\#PK11195

``` r
#collect all information from the filename, which is not included in the PMOD data

pk_strings <- as.tibble(stringr::str_match(list.files(path = "tab_delim_PMOD_pk11195/", pattern = ".kinPar"), pattern = "(\\w*)_(\\d)_(\\w*).kinPar")) %>% 
  rename(filename = V1, subjname = V2, PETNo = V3, model = V4) %>% 
  mutate(tracer = replicate(length(filename), "pk"))

#Read all the files into the tibble created above

pk11195 <- pk_strings %>% 
   group_by(filename) %>% 
   mutate(data = map(filename, ~ read_tsv(file.path("tab_delim_PMOD_pk11195/", .), skip = 5))) %>% 
  ungroup()

#use janitor package to fix names? nah, not necessary, I think.
#First, create a subjname+PETno column. Second, group_by the new column
#Next, spread the data column so that each model gets its own column.

pk11195 <- pk11195 %>% 
  mutate(PET = paste(subjname, PETNo, sep='_')) %>% 
  select(-filename) %>% 
  group_by(PET) %>% 
  spread(model, data) %>% 
  ungroup() %>% 
  rename("two_tcm" =`2tcm`)

# dir.create("RawData_PMOD_RDS")

# saveRDS(pk11195, 'RawData_PMOD_RDS/pmod_pk11195_raw.rds')


#selecting some variables from the model summaries using map ~select

pk11195 <-  pk11195 %>% 
   mutate(two_tcm = map(two_tcm, ~select(.x, "region" = "VoiName(Region) [string]", Vt = "Vt [ml/ccm]", se_Vt = "SE_Vt [%]", "k1" = "K1 [ml/ccm/min]","se_k1" = "SE_K1 [%]", k2 = "k2 [1/min]",
                                         se_k2 = "SE_k2 [%]", k3 = "k3 [1/min]",se_k3  = "SE_k3 [%]", k4 = "k4 [1/min]", se_k4 = "SE_k4 [%]"))) %>% 
  mutate(logan = map(logan, ~select(.x, "region" = "VoiName(Region) [string]", Vt = "Vt [ml/ccm]", se_Vt = "SE_Vt [%]", t_star = "t* [min]"))) %>% 
  mutate(ma1 = map(ma1, ~select(.x, "region" = "VoiName(Region) [string]", Vt = "Vt [ml/ccm]", se_Vt = "SE_Vt [%]", "k1" = "K1 [ml/ccm/min]","se_k1" = "SE_K1 [%]",  t_star = "t* [min]")))

# dir.create("Data_PMOD_RDS")
  
# saveRDS(pk11195, 'Data_PMOD_RDS/pmod_pk11195.rds')
```

\#load pk11195 data

``` r
pk11195 <- readRDS('Data_PMOD_RDS/pmod_pk11195.rds') %>% 
  group_by(subjname) %>% 
  nest() %>% 
  ungroup() %>%
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()
```

\#Pbr28

``` r
#collect all information from the filename, which is not included in the PMOD data

pbr_strings <- as.tibble(stringr::str_match(list.files(path = "tab_delim_PMOD_pbr28/", pattern = ".kinPar"), pattern = "(\\w*)_(\\d)_(\\w*).kinPar")) %>% 
  rename(filename = V1, subjname = V2, PETNo = V3, model = V4) %>% 
  mutate(tracer = replicate(length(filename), "pbr28"))

#Read all the files into the tibble created above

pbr28 <- pbr_strings %>% 
   group_by(filename) %>% 
   mutate(data = map(filename, ~ read_tsv(file.path("tab_delim_PMOD_pbr28/", .), skip = 5))) %>% 
  ungroup

#use janitor package to fix names? nah, not necessary, I think.
#First, create a subjname+PETno column. Second, group_by the new column
#Next, spread the data column so that each model gets its own column.

pbr28 <- pbr28 %>% 
  mutate(PET = paste(subjname, PETNo, sep='_')) %>% 
  select(-filename) %>% 
  group_by(PET) %>% 
  spread(model, data) %>% 
  ungroup() %>% 
  rename("two_tcm" =`2tcm`)

# saveRDS(pbr28, 'RawData_PMOD_RDS/pmod_pbr28_raw.rds')


#selecting some variables from the model summaries using map ~select

pbr28 <- pbr28 %>% 
   mutate(two_tcm = map(two_tcm, ~select(.x, "region" = "VoiName(Region) [string]", Vt = "Vt [ml/ccm]", se_Vt = "SE_Vt [%]", "k1" = "K1 [ml/ccm/min]","se_k1" = "SE_K1 [%]", k2 = "k2 [1/min]",
                                         se_k2 = "SE_k2 [%]", k3 = "k3 [1/min]",se_k3  = "SE_k3 [%]", k4 = "k4 [1/min]", se_k4 = "SE_k4 [%]"))) %>% 
  mutate(logan = map(logan, ~select(.x, "region" = "VoiName(Region) [string]", Vt = "Vt [ml/ccm]", se_Vt = "SE_Vt [%]", t_star = "t* [min]"))) %>% 
  mutate(ma1 = map(ma1, ~select(.x, "region" = "VoiName(Region) [string]", Vt = "Vt [ml/ccm]", se_Vt = "SE_Vt [%]", "k1" = "K1 [ml/ccm/min]","se_k1" = "SE_K1 [%]",  t_star = "t* [min]")))
  
# saveRDS(pbr28, 'Data_PMOD_RDS/pmod_pbr28.rds')
```

\#load pbr28 data

``` r
pbr28 <- readRDS('Data_PMOD_RDS/pmod_pbr28.rds')%>% 
  group_by(subjname) %>% 
  nest() %>% 
  ungroup() %>%
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()
```

\#SCH

``` r
#collect all information from the filename, which is not included in the PMOD data

sch_strings <- as.tibble(stringr::str_match(list.files(path = "tab_delim_PMOD_SCH23390/", pattern = ".kinPar"), pattern = "(\\w*)_(\\d)_(\\w*).kinPar")) %>% 
  rename(filename = V1, subjname = V2, PETNo = V3, model = V4) %>% 
  mutate(tracer = replicate(length(filename), "sch"))

#Read all the files into the tibble created above

sch <- sch_strings %>% 
   group_by(filename) %>% 
   mutate(data = map(filename, ~ read_tsv(file.path("tab_delim_PMOD_SCH23390/", .), skip = 5))) %>% 
  ungroup

#use janitor package to fix names? nah, not necessary, I think.
#First, create a subjname+PETno column. Second, group_by the new column
#Next, spread the data column so that each model gets its own column.

sch <- sch %>% 
  mutate(PET = paste(subjname, PETNo, sep='_')) %>% 
  select(-filename) %>% 
  group_by(PET) %>% 
  spread(model, data) %>% 
  ungroup() 

# saveRDS(sch, 'RawData_PMOD_RDS/pmod_SCH23390_raw.rds')


#selecting some variables from the model summaries using map ~select

sch <- sch %>% 
   mutate(srtm = map(srtm, ~select(.x, "region" = "VoiName(Region) [string]", bp = "BPnd [1/1]", se_bp = "SE_BPnd [%]", k2 = "k2 [1/min]", se_k2 =  "SE_k2 [%]", R1 = "R1 [1/1]", se_R1 = "SE_R1 [%]"))) %>% 
  mutate(mrtm2 = map(mrtm2, ~select(.x, "region" = "VoiName(Region) [string]", bp = "BPnd [1/1]", se_bp = "SE_BPnd [%]", k2_prime = "k2' [1/min]" ,t_star = "t* [min]"))) %>% 
  mutate(ref_logan = map(ref_logan, ~select(.x, "region" = "VoiName(Region) [string]", bp = "BPnd [1/1]", se_bp = "SE_BPnd [%]", k2_prime = "k2' [1/min]" ,t_star = "t* [min]", DVR = "DVR [1/1]", se_DVR = "SE_DVR [%]")))

#filtering out CBL

sch <- sch %>% 
 mutate(srtm = map(srtm, ~filter(.x, region != "CBL"))) %>% 
 mutate(mrtm2 = map(mrtm2, ~filter(.x, region != "CBL"))) %>% 
 mutate(ref_logan = map(ref_logan, ~filter(.x, region != "CBL")))

# saveRDS(sch, 'Data_PMOD_RDS/pmod_SCH23390.rds')
```

\#load sch data

``` r
sch <- readRDS('Data_PMOD_RDS/pmod_SCH23390.rds') %>% 
  group_by(subjname) %>% 
  nest() %>% 
  ungroup() %>%
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup()
```

\#AZ9369

``` r
#collect all information from the filename, which is not included in the PMOD data

az_strings <- as.tibble(stringr::str_match(list.files(path = "tab_delim_PMOD_AZ10419369/", pattern = ".kinPar"), pattern = "(\\w*)_(\\d)_(\\w*).kinPar")) %>% 
  rename(filename = V1, subjname = V2, PETNo = V3, model = V4) %>% 
  mutate(tracer = replicate(length(filename), "az"))

#Read all the files into the tibble created above

az <- az_strings %>% 
   group_by(filename) %>% 
   mutate(data = map(filename, ~ read_tsv(file.path("tab_delim_PMOD_AZ10419369/", .), skip = 5))) %>% 
  ungroup

#use janitor package to fix names? nah, not necessary, I think.
#First, create a subjname+PETno column. Second, group_by the new column
#Next, spread the data column so that each model gets its own column.

az <- az %>% 
  mutate(PET = paste(subjname, PETNo, sep='_')) %>% 
  select(-filename) %>% 
  group_by(PET) %>% 
  spread(model, data) %>% 
  ungroup() 

saveRDS(az, 'RawData_PMOD_RDS/pmod_AZ10419369_raw.rds')


#selecting some variables from the model summaries using map ~select

az <- az %>% 
   mutate(srtm = map(srtm, ~select(.x, "region" = "VoiName(Region) [string]", bp = "BPnd [1/1]", se_bp = "SE_BPnd [%]", k2 = "k2 [1/min]", se_k2 =  "SE_k2 [%]", R1 = "R1 [1/1]", se_R1 = "SE_R1 [%]"))) %>% 
  mutate(mrtm2 = map(mrtm2, ~select(.x, "region" = "VoiName(Region) [string]", bp = "BPnd [1/1]", se_bp = "SE_BPnd [%]", k2_prime = "k2' [1/min]" ,t_star = "t* [min]"))) %>% 
  mutate(ref_logan = map(ref_logan, ~select(.x, "region" = "VoiName(Region) [string]",bp = "BPnd [1/1]", se_bp = "SE_BPnd [%]", k2_prime = "k2' [1/min]" ,t_star = "t* [min]", DVR = "DVR [1/1]", se_DVR = "SE_DVR [%]")))

#filtering out CBL

az <- az %>% 
 mutate(srtm = map(srtm, ~filter(.x, region != "CBL"))) %>% 
 mutate(mrtm2 = map(mrtm2, ~filter(.x, region != "CBL"))) %>% 
 mutate(ref_logan = map(ref_logan, ~filter(.x, region != "CBL")))

saveRDS(az, 'Data_PMOD_RDS/pmod_AZ10419369.rds')
```

\#load az9369 data

``` r
az <- readRDS('Data_PMOD_RDS/pmod_AZ10419369.rds') %>% 
  group_by(subjname) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(id = row_number()) %>% 
  unnest() %>% 
  ungroup() 
```

\#microparameters

``` r
pbr28_micro <- pbr28 %>% 
  select(-logan, -ma1) %>% 
  mutate(two_tcm = map(two_tcm, ~select(.x, region, k1, k2,k3,k4))) %>% 
  unnest() %>% 
  rename(Ligand = tracer, K1 = k1)%>% 
  filter(region == "FC") %>% 
  mutate(Ligand = str_replace( string = Ligand, pattern = 'pbr28', replacement = 'PBR28'))
pk11195_micro <- pk11195 %>% 
  select(-logan, -ma1) %>% 
  mutate(two_tcm = map(two_tcm, ~select(.x, region, k1, k2,k3,k4))) %>% 
  unnest() %>% 
  rename(Ligand = tracer, K1 = k1)%>% 
  filter(region == "FC") %>% 
  mutate(Ligand = str_replace( string = Ligand, pattern = 'pk', replacement = 'PK11195'))
az_micro <- az %>% 
  select(-mrtm2, -ref_logan) %>% 
  mutate(srtm = map(srtm, ~select(.x, region, R1, k2))) %>% 
  unnest() %>% 
  rename(Ligand = tracer)%>% 
  filter(region == "OCC") %>% 
  mutate(region = str_replace( string = region, pattern = 'OCC', replacement = 'OC')) %>% 
  mutate(Ligand = str_replace( string = Ligand, pattern = 'az', replacement = "AZ10419369"))
sch_micro <- sch %>% 
  select(-mrtm2, -ref_logan) %>% 
  mutate(srtm = map(srtm, ~select(.x, region, R1, k2))) %>% 
  unnest() %>% 
  rename(Ligand = tracer)%>% 
  filter(region == "STR") %>% 
  mutate(Ligand = str_replace( string = Ligand, pattern = 'sch', replacement = "SCH23390"))
micro_aif_pmod <- bind_rows(pbr28_micro, pk11195_micro) %>% 
  mutate(software = replicate(length(PET), "PMOD")) %>% 
    mutate(model = replicate(length(PET), "two_tcm")) 
micro_ref_pmod <- bind_rows(az_micro, sch_micro) %>% 
  mutate(software = replicate(length(PET), "PMOD")) %>% 
    mutate(model = replicate(length(PET), "SRTM")) 
# saveRDS(micro_aif_pmod, 'Data_PMOD_RDS/micro_aif_pmod.rds' )
#  
# saveRDS(micro_ref_pmod, 'Data_PMOD_RDS/micro_ref_pmod.rds' )
```

\#Combine all macroparameters for all the tracers

This should work for all figures except the microparameter one

``` r
az <- az %>% 
  mutate(srtm = map(srtm, ~select(.x, region, bp, se_bp))) %>% 
  mutate(mrtm2 = map(mrtm2, ~select(.x, region, bp, se_bp, t_star))) %>% 
  mutate(ref_logan = map(ref_logan, ~select(.x, region, bp, se_bp, 
                t_star))) %>%
  group_by(PET, subjname, PETNo) %>% 
  gather(model, data, mrtm2:srtm) %>% 
  unnest(data) %>% 
  rename(result = bp, SE = se_bp) %>% 
  ungroup() %>% 
  mutate(metric = replicate(length(PET), "bp"))

sch <- sch %>% 
  mutate(srtm = map(srtm, ~select(.x, region, bp, se_bp))) %>% 
  mutate(mrtm2 = map(mrtm2, ~select(.x, region, bp, se_bp, t_star))) %>% 
  mutate(ref_logan = map(ref_logan, ~select(.x, region, bp, se_bp,
                    t_star))) %>%
  group_by(PET, subjname, PETNo) %>% 
  gather(model, data, mrtm2:srtm) %>% 
  unnest(data) %>% 
  rename(result = bp, SE = se_bp) %>% 
  ungroup() %>% 
  mutate(metric = replicate(length(PET), "bp"))

pk11195 <- pk11195 %>% 
  mutate(two_tcm = map(two_tcm, ~select(.x, region, Vt, se_Vt))) %>% 
  mutate(ma1 = map(ma1, ~select(.x, region, Vt, se_Vt, t_star))) %>% 
  mutate(logan = map(logan, ~select(.x, region, Vt, se_Vt, t_star))) %>%
  group_by(PET, subjname, PETNo) %>% 
  gather(model, data, two_tcm:ma1) %>% 
  unnest(data) %>% 
  rename(result = Vt, SE = se_Vt) %>% 
  ungroup() %>% 
  mutate(metric = replicate(length(PET), "Vt"))
  
pbr28 <- pbr28 %>% 
  mutate(two_tcm = map(two_tcm, ~select(.x, region, Vt, se_Vt))) %>% 
  mutate(ma1 = map(ma1, ~select(.x, region, Vt, se_Vt, t_star))) %>% 
  mutate(logan = map(logan, ~select(.x, region, Vt, se_Vt, t_star))) %>%
  group_by(PET, subjname, PETNo) %>% 
  gather(model, data, two_tcm:ma1) %>% 
  unnest(data) %>% 
  rename(result = Vt, SE = se_Vt) %>% 
  ungroup() %>%  
  mutate(metric = replicate(length(PET), "Vt"))

pmod_macro <- bind_rows(az, sch, pbr28, pk11195) %>%
  mutate(software = replicate(length(PET), "pmod")) %>% 
  rename(Value = result, measure = metric)
  
macro <- pmod_macro %>%
  filter(region %in% c("FC", "TC", "THA", "WB", "STR", "OCC"))%>% 
  mutate(region = str_replace( string = region, pattern = 'OCC', replacement = 'OC'))

# saveRDS(macro, 'Data_PMOD_RDS/pmod_macroparameters_more.rds')
```

\#\#Prepare PMOD macroparameter dataset

As different regions are required for different tracers, I will create 4
different datasets, one for each tracer, and then filter by the regions
I want and then combine the datasets once again. OC and FC for “9369”
(az), THA and FC for pbr28 and pk11195 as both are TSPO-ligands. STR and
FC for SCH.

``` r
pmod_macro <- pmod_macro %>%
  ungroup() %>% 
  mutate(region = str_replace( string = region, pattern = 'OCC', replacement = 'OC')) %>% 
  select(-SE) 

pmod_az <- pmod_macro %>% 
  filter(tracer == "az") %>% 
  filter(region %in% c( "OC", "FC"))

pmod_sch <- pmod_macro %>% 
  filter(tracer == "sch") %>% 
  filter(region %in% c( "STR", "FC"))

pmod_pk11195 <- pmod_macro %>% 
  filter(tracer == "pk") %>% 
  filter(region %in% c( "THA", "FC"))

pmod_pbr28 <- pmod_macro %>% 
  filter(tracer == "pbr28") %>% 
  filter(region %in% c( "THA", "FC"))

pmod_macro <- bind_rows(pmod_az, pmod_sch, pmod_pbr28, pmod_pk11195) %>% 
  ungroup()

# saveRDS(pmod_macro, 'Data_PMOD_RDS/pmod_macroparameters_t_star.rds')
```
