-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
-   [Macroparameters - Tables 1 and
    2](#macroparameters---tables-1-and-2)
    -   [Load both datasets](#load-both-datasets)
    -   [Combine datasets, matching
        weights](#combine-datasets-matching-weights)
    -   [Split into reference tissue models and invasive
        models](#split-into-reference-tissue-models-and-invasive-models)
    -   [Comparison between kinfitr and PMOD - models requiring blood
        sampling](#comparison-between-kinfitr-and-pmod---models-requiring-blood-sampling)
        -   [New combined tables 1 and 2](#new-combined-tables-1-and-2)
-   [Session info](#session-info)

Aims
====

The aim of this assignment is to evaluate the performance of kinfitr and
PMOD while using the t\* values computed by PMOD in the kinfitr analysis
and while using constant weights in the kinfitr analysis.

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
kinfitr_macro <- readRDS('Data_PMOD_RDS/kinfitr_macroparameters_tstar.rds') 

pmod_macro <- readRDS('Data_PMOD_RDS/pmod_macroparameters.rds')

pmod_macro <- pmod_macro %>% 
  mutate(PETNo = as.numeric(PETNo)) %>% 
  select(-id)
```

Combine datasets, matching weights
----------------------------------

``` r
kinfitr_macro <- kinfitr_macro %>% 
  filter(model != "mrtm2_tstar" & model != "MRTM2_weights") %>% 
  filter(model != "ma1_tstar" & model != "ma1_weight")

macro <- bind_rows(kinfitr_macro, pmod_macro) %>% 
  rename(Ligand = tracer, Model = model) %>% 
  mutate(Model = str_replace(string = Model, pattern = "mrtm2_tstar_weights", replacement = "MRTM2")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "srtm_weights", replacement = "SRTM")) %>% mutate(Model = str_replace(string = Model, pattern = "ma1_tstar_weight", replacement = "MA1")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "2tcm_weight", replacement = "2TCM")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "mrtm2", replacement = "MRTM2")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "srtm", replacement = "SRTM")) %>% mutate(Model = str_replace(string = Model, pattern = "ma1", replacement = "MA1")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "two_tcm", replacement = "2TCM")) %>%
  mutate(Model = str_replace(string = Model, pattern = "reflogan_tstar", replacement = "ref Logan")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "Logan_tstar", replacement = "Logan")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "ref_logan", replacement = "ref Logan")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "logan", replacement = "Logan")) %>%
  filter(subjname != "rwrd")
```

Split into reference tissue models and invasive models
------------------------------------------------------

``` r
ref <- macro %>% 
  filter(measure == "bp") 

aif <- macro %>% 
  filter(measure == "Vt")
```

\#\#Comparison between kinfitr and PMOD - reference tissue models

``` r
#preparing the data for trt. 

Corr_ref <- ref %>%
  select(Ligand, region, Model, PET, software, Value) %>%
  group_by(Ligand, Model, region) %>%
  spread(software, Value) %>% 
  mutate(Cor = cor(kinfitR, pmod)) %>% 
  gather(key = software, value = Value, kinfitR:pmod) %>% 
  group_by(Ligand, Model, region) %>% 
  summarise(cormean = mean(Cor)) %>% 
  rename(Cor = cormean) %>% 
  mutate(region = case_when(
    region == "OC" ~ "Cor_1",
    region == "STR" ~ "Cor_1",
    TRUE ~ "Cor_2"
  ))  %>% 
  spread(region, Cor) %>% 
  ungroup() 
 
Bias_ref <- ref %>%
  select(Ligand, region, Model, PET, software, Value) %>%
  group_by(Ligand, Model, region) %>%
  spread(software, Value) %>% 
  mutate(Bias = (kinfitR-pmod)/pmod) %>% 
  mutate(Bias = Bias * 100) %>% 
  gather(key = software, value = Value, kinfitR:pmod) %>% 
  group_by(Ligand, Model, region) %>% 
  summarise(Bias = mean(Bias)) %>% 
  mutate(region = case_when(
    region == "OC" ~ "Bias_1",
    region == "STR" ~ "Bias_1",
    TRUE ~ "Bias_2"
  ))  %>% 
  spread(region, Bias) %>% 
  ungroup()  

ref <- ref %>%
  select(Ligand, region, Model, PET, software, Value) %>%
  group_by(Ligand, Model, region) %>%
  nest(.key = "data")

#Using trt from the relfeas package to obtain ICC values
#Note: the values were not arranged by the "subjectname+PETNo" or "PET" variable
#Therefore, I had to use the "rater" argument of trt to specify "software", as 
#the default behaviour is to take the value under as the second rater.

table_1 <- ref %>% 
  group_by(Ligand, Model, region) %>% 
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "PET",
                                       rater = "software")),
         trt_tidy = map(trt, c("tidy"))) %>% 
  select(trt_tidy) %>%  
  unnest() 

#Renaming OC and STR "Region 1" and making FC "Region 2"

table_1 <- table_1 %>%
  select(Ligand, Model, region, icc) %>%
  ungroup() %>% 
  mutate(region = case_when(
    region == "OC" ~ "ICC (Region 1)",
    region == "STR" ~ "ICC (Region 1)",
    TRUE ~ "ICC (Region 2)"
  ))  
 

#Rename the Ligands and then use Spread() on the regions and the ICC

table_1 <- table_1 %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "SCH23390")) %>% 
  spread(region, icc) %>% 
  mutate("Cor (Region 1)" = Corr_ref$Cor_1) %>% 
  mutate("Cor (Region 2)" = Corr_ref$Cor_2) %>% 
  mutate("Bias (Region 1)" = Bias_ref$Bias_1) %>% 
  mutate("Bias (Region 2)" = Bias_ref$Bias_2) %>% 
  arrange(Ligand, desc(Model))

# saveRDS(table_1, "tstars_images/ref_tstar_weights_results.rds")

# kable(table_1, format = "latex", booktabs = T,  digits = 3, escape = F) %>%
#   kable_styling(latex_options = c("striped"),
#                 full_width = F) %>%
#   footnote(symbol = c("Region 1 corresponds to the Occipital cortex in the case of AZ10419369 And the Striatum for SCH23390", "Region 2 is the Frontal cortex for both radioligands"),threeparttable = T)

# save_kable(table_1, 'figs/table_1.pdf')
```

Comparison between kinfitr and PMOD - models requiring blood sampling
---------------------------------------------------------------------

``` r
Corr_aif <- aif %>%
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

Bias_aif <- aif %>%
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

aif <- aif %>% 
  select(Ligand, region, Model, PET, software, Value ) %>% 
  group_by(Ligand, Model, region) %>% 
  nest(.key = "data")

#Using trt from the relfeas package to obtain ICC values
#Note: the values were not arranged by the "subjectname+PETNo" or "PET" variable
#Therefore, I had to use the "rater" argument of trt to specify "software", as 
#the default behaviour is to take the value under as the second rater.

table_2 <- aif %>% 
  group_by(Ligand, Model, region) %>% 
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "PET",
                                       rater = "software")),
         trt_tidy = map(trt, c("tidy"))) %>% 
  select(trt_tidy) %>%  
  unnest() 

#Renaming THA = Thalamus and FC = Frontal Cortex

table_2 <- table_2 %>%
  select(Ligand, Model, region, icc) %>%
  ungroup() %>% 
  mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "ICC (Region 1)")) %>% 
  mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "ICC (Region 2)"))
  


#Rename the Ligands and then use Spread() on the regions and the ICC

table_2 <- table_2 %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>% 
  spread(region, icc) %>% 
  mutate("Cor (Region 1)" = Corr_aif$Cor_1) %>% 
  mutate("Cor (Region 2)" = Corr_aif$Cor_2) %>% 
  mutate("Bias (Region 1)" = Bias_aif$Bias_1) %>% 
  mutate("Bias (Region 2)" = Bias_aif$Bias_2) 

# saveRDS(table_2, "tstars_images/aif_tstar_weights_results.rds")

# kable(table_2, format = "latex", booktabs = T, digits = 3) %>%
#   kable_styling(latex_options = c("striped", "condensed", "scale_down"),
#                 full_width = F) 

# save_kable(table_2, 'figs/table_2.pdf')
```

### New combined tables 1 and 2

``` r
combo <- rbind(table_2, table_1) %>% 
mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[$^{11}$C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[$^{11}$C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[$^{11}$C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[$^{11}$C]PBR28")) %>% 
  filter(Model == "MRTM2" | Model == "MA1")
  
  

# kable(combo, format = "latex", booktabs = T,  digits = 2, escape = F,
#              caption = "Agreement between KinfitR and PMOD",
#       col.names = c("Ligand", "Model", "Region 1", "Region 2", "Region 1",
#                     "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
# kable_styling(full_width = F) %>%
# collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
# pack_rows("Invasive", 1, 6) %>%
# pack_rows("Non-Invasive", 7, 12) %>%
# add_header_above(c(" " = 2, "ICC" = 2, "Pearson's r" = 2, "Bias (%)" = 2))

# kable(combo, format = "latex", booktabs = T,  digits = 2, escape = F,
#              caption = "Agreement between KinfitR and PMOD",
#       col.names = c("Ligand", "Model", "Region 1", "Region 2", "Region 1",
#                     "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
# kable_styling(full_width = F, latex_options = "scale_down") %>%
# collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
# pack_rows("Invasive", 1, 6) %>%
# pack_rows("Non-Invasive", 7, 12) %>%
# add_header_above(c(" " = 2, "ICC" = 2, "Pearson's r" = 2, "Bias (%)" = 2)) %>%
# kable_as_image(.)

#save_kable did not work! Perhaps something in the latex code of kable_extra is
#preventing save_kable from working properly?

mix <- kable(combo, format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "Agreement between KinfitR and PMOD",
      col.names = c("Ligand", "Model", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(full_width = F) %>%
collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
pack_rows("Invasive", 1, 2) %>%
 pack_rows("Non-Invasive", 3, 4) %>%
add_header_above(c(" " = 2, "ICC" = 2, "Pearson's r" = 2, "Bias (%)" = 2))


# save_kable(mix, "table_1_tstar_weights.html")
```

Correlation plots showing the relationship between the different models
for the same Ligand.

``` r
col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))

par(mfrow=c(2,2))

cor_az_oc <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "az") %>% 
  filter(region == "OC") %>% 
  filter(software == "kinfitR") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(AZ10419369 ~ OC),
                 mar=c(0,0,1,0))

cor_az_fc_kinfitR <- macro %>% 
  select(PET, region, Ligand, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "az") %>% 
  filter(region == "FC") %>% 
  filter(software == "kinfitR") %>% 
  select(MRTM2:SRTM) %>% 
  cor() %>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(AZ10419369 ~ FC),
                 mar=c(0,0,1,0))

# cor_az_oc_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>% 
#   spread(Model, Value) %>%
#   filter(Ligand == "az") %>%
#   filter(region == "OC") %>% 
#   filter(software == "pmod") %>% 
#   select(MRTM2:SRTM) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(AZ10419369 ~ OC ~ PMOD),
#                  mar=c(0,0,1,0))
# 
# cor_az_fc_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>% 
#   spread(Model, Value) %>% 
#   filter(Ligand == "az") %>%
#   filter(region == "FC") %>% 
#   filter(software == "pmod") %>% 
#   select(MRTM2:SRTM) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(AZ10419369 ~ FC  ~ PMOD),
#                  mar=c(0,0,1,0))

cor_sch_str_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "sch") %>%
  filter(region == "STR") %>% 
  filter(software == "kinfitR") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(SCH23390 ~ STR),
                 mar=c(0,0,1,0))

cor_sch_fc_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>%  
  spread(Model, Value) %>%
  filter(Ligand == "sch") %>%
  filter(region == "FC") %>% 
  filter(software == "kinfitR") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(SCH23390 ~ FC),
                 mar=c(0,0,1,0))
```

![](kinfitr_pmod_comparison_tstar_weights_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# cor_sch_str_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>%  
#   spread(Model, Value) %>% 
#   filter(Ligand == "sch") %>%
#   filter(region == "STR") %>% 
#   filter(software == "pmod") %>% 
#   select(MRTM2:SRTM) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(SCH23390 ~ STR ~ PMOD),
#                  mar=c(0,0,1,0))
# 
# cor_sch_fc_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>%  
#   spread(Model, Value) %>% 
#   filter(Ligand == "sch") %>%
#   filter(region == "FC") %>% 
#   filter(software == "pmod") %>% 
#   select(MRTM2:SRTM) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(SCH23390 ~ FC ~ PMOD),
#                  mar=c(0,0,1,0))
```

``` r
par(mfrow=c(2,2))

cor_pk_tha_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pk") %>% 
  filter(region == "THA") %>% 
  filter(software == "kinfitR") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PK11195 ~ THA),
                 mar=c(0,0,1,0))

cor_pk_fc_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pk") %>% 
  filter(region == "FC") %>% 
  filter(software == "kinfitR") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PK11195 ~ FC),
                 mar=c(0,0,1,0))

# cor_pk_tha_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>% 
#   spread(Model, Value) %>% 
#   filter(Ligand == "pk") %>% 
#   filter(region == "THA") %>% 
#   filter(software == "pmod") %>% 
#   select("2TCM":MA1) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(PK11195 ~ THA ~ PMOD),
#                  mar=c(0,0,1,0))
# 
# cor_pk_fc_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>% 
#   spread(Model, Value) %>% 
#   filter(Ligand == "pk") %>% 
#   filter(region == "FC") %>% 
#   filter(software == "pmod") %>% 
#   select("2TCM":MA1) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(PK11195 ~ FC ~ PMOD),
#                  mar=c(0,0,1,0))

cor_pbr28_tha_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pbr28") %>% 
  filter(region == "THA") %>% 
  filter(software == "kinfitR") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PBR28 ~ THA),
                 mar=c(0,0,1,0))

cor_pbr28_fc_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pbr28") %>% 
  filter(region == "FC") %>% 
  filter(software == "kinfitR") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PBR28 ~ FC),
                 mar=c(0,0,1,0))
```

![](kinfitr_pmod_comparison_tstar_weights_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
# cor_pbr28_tha_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>% 
#   spread(Model, Value) %>% 
#   filter(Ligand == "pbr28") %>% 
#   filter(region == "THA") %>% 
#   filter(software == "pmod") %>% 
#   select("2TCM":MA1) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(PBR28 ~ THA ~ PMOD),
#                  mar=c(0,0,1,0))
# 
# cor_pbr28_fc_pmod <- macro %>% 
#   select(Ligand, PET, region, Model, Value, software) %>% 
#   spread(Model, Value) %>% 
#   filter(Ligand == "pbr28") %>% 
#   filter(region == "FC") %>% 
#   filter(software == "pmod") %>% 
#   select("2TCM":MA1) %>% 
#   cor()%>%
#   corrplot.mixed(lower='ellipse', upper='number', 
#                  lower.col = col2(200), upper.col = col2(200),  diag='n',
#                  number.digits = 2, title=expression(PBR28 ~ FC ~ PMOD),
#                  mar=c(0,0,1,0))
```

\#median of the correlations

``` r
median_invasive_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(software == "pmod") %>% 
  filter(Ligand == "pbr28" | Ligand == "pk") %>% 
  select(Ligand:MA1)

median_invasive_pmod %>% 
  select(`2TCM`:MA1) %>% 
  correlate(quiet = T) %>% 
  shave() %>% 
  stretch() %>% 
  drop_na() %>% 
  summarize(Median = median(r)) %>% 
  kable(digits = 2, caption = "Median correlation of the invasive models with PMOD")
```

<table>
<caption>
Median correlation of the invasive models with PMOD
</caption>
<thead>
<tr>
<th style="text-align:right;">
Median
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.99
</td>
</tr>
</tbody>
</table>

``` r
median_invasive_kinfitr <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(software == "kinfitR") %>% 
  filter(Ligand == "pbr28" | Ligand == "pk") %>% 
  select(Ligand:MA1)

median_invasive_kinfitr %>% 
  select(`2TCM`:MA1) %>% 
  correlate(quiet = T) %>% 
  shave() %>% 
  stretch() %>% 
  drop_na() %>% 
  summarize(Median = median(r)) %>% 
  kable(digits = 2, caption = "Median correlation of the invasive models with kinfitr")
```

<table>
<caption>
Median correlation of the invasive models with kinfitr
</caption>
<thead>
<tr>
<th style="text-align:right;">
Median
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.99
</td>
</tr>
</tbody>
</table>

``` r
median_ref_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(software == "pmod") %>% 
  filter(Ligand == "az" | Ligand == "sch") %>% 
  select(MRTM2:SRTM)

median_ref_pmod %>% 
  correlate(quiet = T) %>% 
  shave() %>% 
  stretch() %>% 
  drop_na() %>% 
  summarize(Median = median(r)) %>% 
  kable(digits = 2, caption = "Median correlation of the non-invasive models with PMOD")
```

<table>
<caption>
Median correlation of the non-invasive models with PMOD
</caption>
<thead>
<tr>
<th style="text-align:right;">
Median
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.99
</td>
</tr>
</tbody>
</table>

``` r
median_ref_kinfitr <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(software == "kinfitR") %>% 
  filter(Ligand == "az" | Ligand == "sch") %>% 
  select(MRTM2:SRTM)

median_ref_kinfitr %>% 
  correlate(quiet = T) %>% 
  shave() %>% 
  stretch() %>% 
  drop_na() %>% 
  summarize(Median = median(r)) %>% 
  kable(digits = 2, caption = "Median correlation of the non-invasive models with kinfitr")
```

<table>
<caption>
Median correlation of the non-invasive models with kinfitr
</caption>
<thead>
<tr>
<th style="text-align:right;">
Median
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.99
</td>
</tr>
</tbody>
</table>

Session info
============

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 18363)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Sweden.1252  LC_CTYPE=English_Sweden.1252   
    ## [3] LC_MONETARY=English_Sweden.1252 LC_NUMERIC=C                   
    ## [5] LC_TIME=English_Sweden.1252    
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] relfeas_0.0.2         granviller_0.0.0.9000 kipettools_0.0.0.9000
    ##  [4] kinfitr_0.5           webshot_0.5.2         magick_2.3           
    ##  [7] corrr_0.4.1           kableExtra_1.1.0      forcats_0.5.0        
    ## [10] dplyr_0.8.5           purrr_0.3.3           readr_1.3.1          
    ## [13] tidyr_1.0.2           tibble_3.0.0          ggplot2_3.3.0        
    ## [16] tidyverse_1.3.0       janitor_2.0.1         viridis_0.5.1        
    ## [19] viridisLite_0.3.0     cowplot_1.0.0         knitr_1.28           
    ## [22] rjags_4-10            coda_0.19-3           lme4_1.1-23          
    ## [25] Matrix_1.2-18         pracma_2.2.9          readxl_1.3.1         
    ## [28] psych_1.9.12.31       RColorBrewer_1.1-2    gridExtra_2.3        
    ## [31] corrplot_0.84         stringr_1.4.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] mcmc_0.9-7         nlme_3.1-142       fs_1.3.1           lubridate_1.7.8   
    ##  [5] httr_1.4.1         tools_3.6.2        backports_1.1.5    R2WinBUGS_2.1-21  
    ##  [9] R6_2.4.1           DBI_1.1.0          colorspace_1.4-1   withr_2.1.2       
    ## [13] tidyselect_1.0.0   mnormt_1.5-6       compiler_3.6.2     quantreg_5.55     
    ## [17] cli_2.0.2          rvest_0.3.5        SparseM_1.78       xml2_1.2.2        
    ## [21] scales_1.1.0       mvtnorm_1.1-0      digest_0.6.25      minqa_1.2.4       
    ## [25] rmarkdown_2.1      MCMCpack_1.4-6     pkgconfig_2.0.3    htmltools_0.4.0   
    ## [29] highr_0.8          dbplyr_1.4.2       rlang_0.4.5        rstudioapi_0.11   
    ## [33] generics_0.0.2     jsonlite_1.6.1     magrittr_1.5       agRee_0.5-3       
    ## [37] Rcpp_1.0.3         munsell_0.5.0      fansi_0.4.1        abind_1.4-5       
    ## [41] lifecycle_0.2.0    stringi_1.4.6      yaml_2.2.1         snakecase_0.11.0  
    ## [45] MASS_7.3-51.5      parallel_3.6.2     crayon_1.3.4       lattice_0.20-38   
    ## [49] haven_2.2.0        splines_3.6.2      hms_0.5.3          pillar_1.4.3      
    ## [53] boot_1.3-23        R2jags_0.5-7       reprex_0.3.0       glue_1.3.1        
    ## [57] evaluate_0.14      modelr_0.1.5       selectr_0.4-2      vctrs_0.2.4       
    ## [61] nloptr_1.2.2.1     MatrixModels_0.4-1 cellranger_1.1.0   gtable_0.3.0      
    ## [65] assertthat_0.2.1   xfun_0.12          broom_0.5.5        miscF_0.1-5       
    ## [69] statmod_1.4.34     ellipsis_0.3.0
