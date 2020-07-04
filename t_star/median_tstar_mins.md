-   [invasive models and reference tissue
    models](#invasive-models-and-reference-tissue-models)

``` r
library(tidyverse)
library(knitr)
library(kableExtra)
```

``` r
pmod_tstar <- readRDS('RawData_PMOD_RDS/pmod_macroparameters_t_star.rds')
```

invasive models and reference tissue models
===========================================

``` r
macro_aif <- pmod_tstar %>% 
  filter(measure == "Vt") %>% 
  filter(model != "two_tcm")

macro_ref <- pmod_tstar %>% 
  filter(measure == "bp") %>% 
  filter(model != "srtm")
```

``` r
macro_aif <- macro_aif %>% 
  select(-c(measure, id, software, Value)) %>% 
  group_by(PET, subjname, PETNo, tracer, region) %>% 
  spread(key = model, value = t_star) %>% 
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
  select(-c(measure, id, software, Value)) %>% 
  group_by(PET, subjname, PETNo, tracer, region) %>% 
  spread(key = model, value = t_star) %>% 
  rename(Region = region, Ligand = tracer) %>% 
  group_by(Region, Ligand) %>% 
  summarize(median_ref_logan = median(ref_logan), iqr_ref_logan = IQR(ref_logan), min_ref_logan = min(ref_logan), max_ref_logan = max(ref_logan), median_mrtm2 = median(mrtm2), iqr_mrtm2= IQR(mrtm2),
            min_mrtm2 = min(mrtm2),max_mrtm2 = max(mrtm2)) %>% 
  select(Ligand, Region, everything()) %>% 
  arrange(Ligand, desc(Region)) 

macro_ref <- kable(macro_ref, format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "t* values fitted by PMOD for the non_invasive models",
      col.names = c("Ligand", "Region", "Median", "IQR", "Min",
                    "Max", "Median", "IQR", "Min",
                    "Max"), align= 'c') %>%
kable_styling(full_width = F) %>%
add_header_above(c(" " = 2, "ref Logan" = 4, "MRTM2" = 4))

save_kable(macro_ref, "t_star_medians_ref.html")
```
