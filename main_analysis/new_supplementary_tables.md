-   [Load packages](#load-packages)
-   [Load data](#load-data)

Load packages
=============

``` r
library(tidyverse)
library(knitr)
library(kableExtra)
```

Load data
=========

``` r
orig_aif <- readRDS("tstars_images/aif_orig_results.rds")

tstar_aif <- readRDS("tstars_images/aif_tstar_results.rds") %>% 
  mutate("t_star" = "Fitted")

orig_ref <- readRDS("tstars_images/ref_orig_results.rds")

tstar_ref <- readRDS("tstars_images/ref_tstar_results.rds") %>% 
  mutate("t_star" = "Fitted")

weights_aif <- readRDS("tstars_images/aif_tstar_weights_results.rds") %>% 
  mutate("Weights" = "Constant")

weights_ref <- readRDS("tstars_images/ref_tstar_weights_results.rds") %>% 
  mutate("Weights" = "Constant")
```

combine datasets

``` r
tstar_aif %>% 
  bind_rows(orig_aif) %>% 
  mutate(t_star = ifelse(is.na(t_star) == TRUE, "Selected", t_star)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, Model, desc(t_star)) %>% 
  select(Ligand, Model, "t*" = t_star, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "2TCM") %>% 
  select(Ligand, Model, "t*", contains("Cor"), everything()) %>% 
  kable(., format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "Agreement between kinfitr and PMOD using fitted or selected t* values for the invasive models",
      col.names = c("Ligand", "Model", "t*", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(full_width = F, latex_options = "scale_down") %>%
collapse_rows(columns = c(1,2), latex_hline = "major", valign = "middle") %>% 
add_header_above(c(" " = 3, "Pearson's r" = 2, "ICC" = 2, "Bias (%)" = 2)) %>% 
  save_kable(., "tstars_images/agreement_aif_tstar.html")
```

``` r
tstar_ref %>% 
  bind_rows(orig_ref) %>% 
  mutate(t_star = ifelse(is.na(t_star) == TRUE, "Selected", t_star)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, desc(Model), desc(t_star)) %>% 
  select(Ligand, Model, "t*" = t_star, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "SRTM") %>% 
  select(Ligand, Model, "t*", contains("Cor"), everything()) %>%
  kable(., format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "Agreement between kinfitr and PMOD using fitted or selected t* values for the non-invasive models",
      col.names = c("Ligand", "Model", "t*", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(full_width = F, latex_options = "scale_down") %>%
collapse_rows(columns = c(1,2), latex_hline = "major", valign = "middle") %>% 
add_header_above(c(" " = 3, "Pearson's r" = 2, "ICC" = 2, "Bias (%)" = 2)) %>%
  save_kable(., "tstars_images/agreement_ref_tstar.html")
```

Improvements

``` r
improv_aif <- tstar_aif %>% 
  bind_rows(orig_aif) %>% 
  mutate(t_star = ifelse(is.na(t_star) == TRUE, "Selected", t_star)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, Model, desc(t_star)) %>% 
  select(Ligand, Model, "t*" = t_star, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "2TCM") %>% 
  select(Ligand, Model, "t*", contains("Cor"), everything())

improv_ref <- tstar_ref %>% 
  bind_rows(orig_ref) %>% 
  mutate(t_star = ifelse(is.na(t_star) == TRUE, "Selected", t_star)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, desc(Model), desc(t_star)) %>% 
  select(Ligand, Model, "t*" = t_star, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "SRTM") %>% 
  select(Ligand, Model, "t*", contains("Cor"), everything())

improv <- bind_rows(improv_aif, improv_ref)

improv_reg1 <- select(improv, Ligand:"t*", contains("Region 1")) %>% 
  mutate(Region = "Region 1")
improv_reg2 <- select(improv, Ligand:"t*", contains("Region 2")) %>% 
  mutate(Region = "Region 2")

colnames(improv_reg1) <- str_remove(colnames(improv_reg1),
                                    " \\(.*")

colnames(improv_reg2) <- str_remove(colnames(improv_reg2),
                                    " \\(.*")

improv_long <- bind_rows(improv_reg1, improv_reg2)

improv_vals <- improv_long %>% 
  group_by(Ligand, Model, Region) %>% 
  summarise(win_cor= Cor[2] - Cor[1],
            win_ICC = ICC[2] - ICC[1],
            win_Bias = abs(Bias[1]) - abs(Bias[2])) %>% 
  ungroup() %>% 
  summarise(win_cor = mean(win_cor),
            win_ICC = mean(win_ICC),
            win_Bias = mean(win_Bias))
  

improv_comp <- improv_long %>% 
  group_by(Ligand, Model, Region) %>% 
  summarise("winner_cor"= `t*`[which.max(Cor)],
            winner_ICC = `t*`[which.max(ICC)],
            winner_Bias = `t*`[which.min(abs(Bias))])

improv_num_reg <- improv_comp %>%
  ungroup() %>% 
  group_by(Ligand) %>% 
  summarise(Fitted_win_cor = sum(winner_cor=="Fitted"),
            Fitted_win_ICC = sum(winner_ICC=="Fitted"),
            Fitted_win_Bias = sum(winner_Bias=="Fitted"))

improv_num <- improv_comp %>%
  ungroup() %>% 
  summarise(Fitted_win_cor = sum(winner_cor=="Fitted"),
            Fitted_win_ICC = sum(winner_ICC=="Fitted"),
            Fitted_win_Bias = sum(winner_Bias=="Fitted"))
```

``` r
weights_ref %>% 
  bind_rows(orig_ref) %>% 
  mutate(Weights = ifelse(is.na(Weights) == TRUE, "Computed", Weights)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, desc(Model), desc(Weights)) %>% 
  select(Ligand, Model, Weights, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "ref Logan") %>% 
  select(Ligand, Model, Weights, contains("Cor"), everything()) %>%
  kable(., format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "Agreement between kinfitr and PMOD using constant or computed weights for the non-invasive models",
      col.names = c("Ligand", "Model", "t*", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(full_width = F, latex_options = "scale_down") %>%
collapse_rows(columns = c(1,2), latex_hline = "major", valign = "middle") %>% 
add_header_above(c(" " = 3, "Pearson's r" = 2, "ICC" = 2, "Bias (%)" = 2)) %>%
  save_kable(., "tstars_images/agreement_ref_weights.html")
```

``` r
weights_aif %>% 
  bind_rows(orig_aif) %>% 
  mutate(Weights = ifelse(is.na(Weights) == TRUE, "Computed", Weights)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, desc(Model), desc(Weights)) %>% 
  select(Ligand, Model, Weights, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "Logan") %>% 
  select(Ligand, Model, Weights, contains("Cor"), everything()) %>%
  kable(., format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "Agreement between kinfitr and PMOD using constant or computed weights for the invasive models",
      col.names = c("Ligand", "Model", "t*", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(full_width = F, latex_options = "scale_down") %>%
collapse_rows(columns = c(1,2), latex_hline = "major", valign = "middle") %>% 
add_header_above(c(" " = 3, "Pearson's r" = 2, "ICC" = 2, "Bias (%)" = 2)) %>%
  save_kable(., "tstars_images/agreement_aif_weights.html")
```

Improvement

``` r
improv_weights_ref <- weights_ref %>% 
  bind_rows(orig_ref) %>% 
  mutate(Weights = ifelse(is.na(Weights) == TRUE, "Computed", Weights)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, desc(Model), desc(Weights)) %>% 
  select(Ligand, Model, Weights, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "ref Logan") %>% 
  select(Ligand, Model, Weights, contains("Cor"), everything())

improv_weights_aif <- weights_aif %>% 
  bind_rows(orig_aif) %>% 
  mutate(Weights = ifelse(is.na(Weights) == TRUE, "Computed", Weights)) %>% 
  group_by(Ligand, Model) %>% 
  arrange(Ligand, desc(Model), desc(Weights)) %>% 
  select(Ligand, Model, Weights, everything()) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  filter(Model != "Logan") %>% 
  select(Ligand, Model, Weights, contains("Cor"), everything())


improv_weights <- bind_rows(improv_weights_aif, improv_weights_ref)

improv_weights_reg1 <- select(improv_weights, Ligand:"Weights", contains("Region 1")) %>% 
  mutate(Region = "Region 1")
improv_weights_reg2 <- select(improv_weights, Ligand:"Weights", contains("Region 2")) %>% 
  mutate(Region = "Region 2")

colnames(improv_weights_reg1) <- str_remove(colnames(improv_weights_reg1),
                                    " \\(.*")

colnames(improv_weights_reg2) <- str_remove(colnames(improv_weights_reg2),
                                    " \\(.*")

improv_weights_long <- bind_rows(improv_weights_reg1, improv_weights_reg2)

improv_weights_vals <- improv_weights_long %>%
  arrange(Weights) %>% 
  group_by(Ligand, Model, Region) %>% 
  summarise(win_cor= Cor[2] - Cor[1],
            win_ICC = ICC[2] - ICC[1],
            win_Bias = abs(Bias[1]) - abs(Bias[2])) %>% 
  ungroup() %>% 
  #group_by(Ligand) %>% 
  summarise(win_cor = mean(win_cor),
            win_ICC = mean(win_ICC),
            win_Bias = mean(win_Bias))
```
