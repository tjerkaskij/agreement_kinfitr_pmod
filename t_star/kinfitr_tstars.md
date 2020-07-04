``` r
library(tidyverse)
library(knitr)
library(kableExtra)
```

``` r
az <- readRDS('tstars_images/az_tstar_kinfitr.rds')

sch <- readRDS('tstars_images/sch_tstar_kinfitr.rds') 

pk <- readRDS('tstars_images/pk11195_tstar_kinfitr.rds')

pbr <- readRDS('tstars_images/pbr28_tstar_kinfitr.rds')
```

``` r
ref <- bind_rows(az,sch) %>% 
  gather(Model, "t*", 'ref Logan', MRTM2)

aif <- bind_rows(pk,pbr) %>% 
  gather(Model, "t*", Logan, MA1)

tstars <- bind_rows(ref, aif)
```

``` r
tstars %>% 
  kable(., format = "html", booktabs = T,  digits = 3, 
             caption = "t* values selected for use in the kinfitr analysis") %>% 
  save_kable(., "tstars_images/kinfitr_tstars.html")
```
