-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
-   [Macroparameters - Tables 1 and
    2](#macroparameters---tables-1-and-2)
    -   [Load both datasets](#load-both-datasets)
    -   [Comparison between kinfitr and PMOD - models requiring blood
        sampling](#comparison-between-kinfitr-and-pmod---models-requiring-blood-sampling)

Aims
====

The aim of this code is to evaluate the performance of kinfitr and PMOD
without the use of iteration when fitting 2TCM in the kinfitr analysis.

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
install.packages("magick")
install.packages("webshot")
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
library(corrr)
library(magick)
library(webshot)
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

Macroparameters - Tables 1 and 2
================================

Load both datasets
------------------

``` r
kinfitr_macro_no_iter <- readRDS('Data_PMOD_RDS/kinfitr_macroparameters_no_iter.rds') %>% filter(model == "two_tcm") %>% 
  mutate(software = str_replace(string = software, pattern = "no_iter", replacement = "pmod")) 
  
aif_weights <- readRDS('Data_PMOD_RDS/kinfitr_macroparameters.rds') %>%
  filter(measure == "Vt") %>% 
  filter(model == "two_tcm") %>% 
  bind_rows(kinfitr_macro_no_iter) %>% 
  rename(Ligand = tracer, Model = model)
```

Comparison between kinfitr and PMOD - models requiring blood sampling
---------------------------------------------------------------------

``` r
Corr_aif_weights <- aif_weights %>%
  select(Ligand, region, Model, PET, software, Value) %>%
  group_by(Ligand, Model, region) %>%
  spread(software, Value) %>% 
 mutate(Cor = cor(kinfitR, pmod)) %>% 
  gather(key = software, value = Value, kinfitR:pmod) %>% 
  group_by(Ligand, Model, region) %>% 
  summarise(cormean = mean(Cor)) %>% 
  rename(Cor = cormean) %>% 
  mutate(region = ifelse(region == "THA", "Cor_1", "Cor_2"))%>% 
  spread(region, Cor) %>% 
  ungroup() 

Bias_aif_weights <- aif_weights %>%
  select(Ligand, region, Model, PET, software, Value) %>%
  group_by(Ligand, Model, region) %>%
  spread(software, Value) %>% 
  mutate(Bias = (kinfitR-pmod)/pmod) %>% 
  mutate(Bias = Bias * 100) %>% 
  gather(key = software, value = Value, kinfitR:pmod) %>% 
  group_by(Ligand, Model, region) %>% 
  summarise(Bias = mean(Bias)) %>% 
  mutate(region = ifelse(region == "THA", "Bias_1", "Bias_2"))%>% 
  spread(region, Bias) %>% 
  ungroup()

aif_weights <- aif_weights %>% 
  select(Ligand, region, Model, PET, software, Value ) %>% 
  group_by(Ligand, Model, region) %>% 
  nest(.key = "data")

#Using trt from the relfeas package to obtain ICC values
#Note: the values were not arranged by the "subjectname+PETNo" or "PET" variable
#Therefore, I had to use the "rater" argument of trt to specify "software", as 
#the default behaviour is to take the value under as the second rater.

table_2_weights <- aif_weights %>% 
  group_by(Ligand, Model, region) %>% 
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "PET",
                                       rater = "software")),
         trt_tidy = map(trt, c("tidy"))) %>% 
  select(trt_tidy) %>%  
  unnest() 

#Renaming THA = Thalamus and FC = Frontal Cortex

table_2_weights <- table_2_weights %>%
  select(Ligand, Model, region, icc) %>%
  ungroup() %>% 
  mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "ICC (Region 1)")) %>% 
  mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "ICC (Region 2)"))
  

#Rename the Ligands and then use Spread() on the regions and the ICC

table_2_weights <- table_2_weights %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>% 
  spread(region, icc) %>% 
  mutate("Cor (Region 1)" = Corr_aif_weights$Cor_1) %>% 
  mutate("Cor (Region 2)" = Corr_aif_weights$Cor_2) %>% 
  mutate("Bias (Region 1)" = Bias_aif_weights$Bias_1) %>% 
  mutate("Bias (Region 2)" = Bias_aif_weights$Bias_2) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[$^{11}$C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[$^{11}$C]PBR28")) %>% 
  select(Ligand, Model, contains("Cor"), everything())

table_2 <- kable(table_2_weights, format = "html", booktabs = T,  digits = 2, 
                 escape = F,
             caption = "Agreement between 2TCM using KinfitR with and without iteration",
      col.names = c("Ligand", "Model", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(full_width = F) %>%
collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>% 
  add_header_above(c(" " = 2, "Pearson's r" = 2,  "ICC" = 2, "Bias (%)" = 2))

# save_kable(table_2, 'figs/table_1_no_iter.html')
```

means

``` r
median_vt <- aif_weights %>% 
  unnest() %>% 
  group_by(Ligand, region, software) %>% 
  summarize(median = median(Value), IQR = IQR(Value), Min = min(Value), 
            Max = max(Value)) %>% 
  ungroup() %>% 
  mutate(software = ifelse(software == "kinfitR", "Yes", "No")) %>% 
  rename(Iteration = software, Region = region, "Median V[t]" = median)

# median_vt <- kable(median_vt, format = "html", booktabs = T,  digits = 3,
#                  escape = F,
#              caption = "Median Vt values with and without iteration for 2TCM using kinfitR") %>%
#   save_kable(., 'figs/vt_medians_no_iter.html')
```
