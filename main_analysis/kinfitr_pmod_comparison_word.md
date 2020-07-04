-   [Aims](#aims)
-   [Libraries](#libraries)
    -   [CRAN libraries](#cran-libraries)
    -   [Non-CRAN libraries](#non-cran-libraries)
        -   [Loading Non\_CRAN libraries and setting
            theme](#loading-non_cran-libraries-and-setting-theme)
-   [Macroparameters - Tables 1 and
    2](#macroparameters---tables-1-and-2)
    -   [Load both datasets](#load-both-datasets)
    -   [Combine datasets](#combine-datasets)
    -   [Fig 2, macroparameter dotplot with
        abline](#fig-2-macroparameter-dotplot-with-abline)
    -   [Split into reference tissue models and invasive
        models](#split-into-reference-tissue-models-and-invasive-models)
    -   [Comparison between kinfitr and PMOD - models requiring blood
        sampling](#comparison-between-kinfitr-and-pmod---models-requiring-blood-sampling)
        -   [New combined tables 1 and 2](#new-combined-tables-1-and-2)
        -   [Comparison of the models](#comparison-of-the-models)
-   [Table 3 - consistency
    test\_retest](#table-3---consistency-test_retest)
-   [Spaghetti plot - The occipital
    cortex](#spaghetti-plot---the-occipital-cortex)
-   [Microparameters](#microparameters)
    -   [Load microparameters](#load-microparameters)
    -   [Combine datasets](#combine-datasets-1)
    -   [Fig 2, microparameter dotplot with abline for
        srtm](#fig-2-microparameter-dotplot-with-abline-for-srtm)
    -   [Fig 3, microparameter dotplot with abline for
        2tcm](#fig-3-microparameter-dotplot-with-abline-for-2tcm)
-   [New MRTM2](#new-mrtm2)
-   [Session info](#session-info)

Aims
====

The aim of this assignment is to evaluate the performance of kinfitr and
PMOD

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
kinfitr_macro <- readRDS('Data_PMOD_RDS/kinfitr_macroparameters.rds') 

pmod_macro <- readRDS('Data_PMOD_RDS/pmod_macroparameters.rds')
```

Combine datasets
----------------

``` r
pmod_macro <- pmod_macro %>% 
  mutate(PETNo = as.numeric(PETNo))

macro <- bind_rows(kinfitr_macro, pmod_macro) %>% 
  rename(Ligand = tracer, Model = model) %>% 
  mutate(Model = str_replace(string = Model, pattern = "mrtm2", replacement = "MRTM2")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "srtm", replacement = "SRTM")) %>% mutate(Model = str_replace(string = Model, pattern = "ma1", replacement = "MA1")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "two_tcm", replacement = "2TCM")) %>% mutate(Model = str_replace(string = Model, pattern = "ref_logan", replacement = "ref Logan")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "logan", replacement = "Logan"))
```

\#outlier violin plots

Plots showing the outlier, which was excluded from the final analysis

``` r
outlier <- macro %>% 
  filter(subjname == "rwrd") %>% 
  filter(region == "FC") %>% 
  filter(Ligand == "pbr28") %>% 
  filter(Model == "2TCM") %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD"))  
  
par(mfrow=c(3,1))

macro %>% 
  filter(Ligand == "pbr28") %>% 
  filter(Model == "2TCM") %>%
  filter(region == "FC") %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD"))  %>% 
  ggplot(aes(x = region, y = Value, colour = software, fill = software)) +
  geom_violin(alpha=0.25)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap( ~ software) +
  geom_dotplot(data = outlier, aes(x = region, y = Value, fill = "darkred", color = "darkred"),binaxis='y', stackdir='center', dotsize=1) +
 scale_fill_discrete(name = "Legend:", labels = c("Outlier", "Kinfitr", "PMOD")) +
 scale_color_discrete(name = "Legend:", labels = c("Outlier", "Kinfitr", "PMOD"))+
ggtitle("2TCM for the Frontal Cortex for PBR28")
  
macro %>% 
  filter(Ligand == "pbr28") %>% 
  filter(Model == "2TCM") %>%
  filter(region == "FC") %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD"))  %>% 
  ggplot(aes(x = region, y = Value, colour = software, fill = software)) +
  geom_violin(alpha=0.25)+
  geom_point() +
  facet_wrap( ~ software) +
  geom_point(data = outlier, aes(x = region, y = Value, fill = "darkred", color = "darkred")) +
 scale_fill_discrete(name = "Legend:", labels = c("Outlier", "kinfitr", "PMOD")) +
 scale_color_discrete(name = "Legend:", labels = c("Outlier", "kinfitr", "PMOD"))+
ggtitle("2TCM for the Frontal Cortex for PBR28")

macro %>% 
  filter(Ligand == "pbr28") %>% 
  filter(Model == "2TCM") %>%
  filter(region == "FC") %>%
  filter(subjname != "rwrd") %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD")) %>% 
  ggplot(aes(x = region, y = Value, colour = software, fill = software)) +
  geom_violin(alpha=0.25)+
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  facet_wrap( ~ software) +
  geom_jitter(data = outlier, aes(x = region, y = Value, fill = "darkred", color = "darkred"),shape=16, position=position_jitter(0.2))+
 scale_fill_discrete(name = "Legend:", labels = c("Outlier", "KinfitR", "PMOD")) +
 scale_color_discrete(name = "Legend:", labels = c("Outlier", "KinfitR", "PMOD"))+
ggtitle("2TCM for the Frontal Cortex for PBR28")
```

\#exclusion of the outlier

``` r
macro <- macro %>%
  filter(subjname != "rwrd") %>% 
  select(-id)
```

Violin plots for pbr28 thalamus region

``` r
 macro %>% 
  filter(Ligand == "pbr28") %>% 
  filter(Model == "2TCM") %>%
  filter(region == "THA") %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD"))  %>% 
  ggplot(aes(x = region, y = Value, colour = software, fill = software)) +
  geom_violin(alpha=0.25)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +
  facet_wrap( ~ software)  +
ggtitle("2TCM for the Frontal Cortex for PBR28")
```

``` r
macro_ref <- macro %>% 
  filter(measure == "bp") %>% 
  group_by(PET, subjname, PETNo, Ligand, region, measure, Model) %>% 
  filter(software == "kinfitR") %>% 
  spread(value = Value, key = software) %>% 
  group_by(Ligand, subjname) %>% 
  nest() %>% 
  ungroup() %>%
  group_by(Ligand) %>% 
  mutate(id = row_number()) %>% 
  unnest(cols = c(data)) %>% 
  mutate(id = as_factor(id))

macro_ref <- macro %>% 
  filter(measure == "bp") %>% 
  group_by(PET, subjname, PETNo, Ligand, region, measure, Model) %>% 
  filter(software == "pmod") %>% 
  spread(value = Value, key = software) %>% 
  right_join(macro_ref) %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "SCH23390")) %>% 
  mutate(Ligand = factor(Ligand, levels = c("AZ10419369","SCH23390"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)"), ordered = TRUE)) %>% 
  mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>% 
  mutate(region = str_replace(string = region, pattern = "OC", 
                              replacement = "Occipital cortex")) %>%
  mutate(region = str_replace(string = region, pattern = "STR", 
                              replacement = "Striatum"))

# macro_ref_id <- macro_ref %>% 
#   filter(Ligand == "sch") %>% 
#   group_by(subjname) %>% 
#   nest() %>% 
#   ungroup() %>% 
#   mutate(id = row_number()) %>% 
#   unnest()
# 
# macro_ref_id <- macro_ref %>% 
#   filter(Ligand == "az") %>% 
#   group_by(subjname) %>% 
#   nest() %>% 
#   ungroup() %>% 
#   mutate(id = row_number()) %>% 
#    unnest() %>% 
#    bind_rows(macro_ref_id) %>% 
#     ungroup() 
 
macro_aif <- macro %>% 
  filter(measure == "Vt") %>% 
  group_by(PET, subjname, PETNo, Ligand, region, measure, Model) %>% 
  filter(software == "kinfitR") %>% 
  spread(value = Value, key = software) %>% 
  group_by(Ligand, subjname) %>% 
  nest() %>% 
  ungroup() %>%
  group_by(Ligand) %>% 
  mutate(id = row_number()) %>% 
  unnest(cols = c(data)) %>% 
  mutate(id = as_factor(id))

macro_aif <- macro %>% 
  filter(measure == "Vt") %>% 
  group_by(PET, subjname, PETNo, Ligand, region, measure, Model) %>% 
  filter(software == "pmod") %>% 
  spread(value = Value, key = software) %>% 
  right_join(macro_aif) %>% 
  ungroup() %>% 
   mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>% 
   mutate(Ligand = factor(Ligand, levels = c("PBR28","PK11195"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)"), ordered = TRUE)) %>% 
  mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>% 
  mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus"))


  
# macro_aif_id <- macro_aif %>% 
#   filter(Ligand == "pbr28") %>% 
#   group_by(subjname) %>% 
#   nest() %>% 
#   ungroup() %>% 
#   mutate(id = row_number()) %>% 
#   unnest()
# 
# macro_aif_id <- macro_aif %>% 
#   filter(Ligand == "pk") %>% 
#   group_by(subjname) %>% 
#   nest() %>% 
#   ungroup() %>% 
#   mutate(id = row_number()) %>% 
#    unnest() %>% 
#    bind_rows(macro_aif_id) %>% 
#     ungroup() 
```

Fig 2, macroparameter dotplot with abline
-----------------------------------------

``` r
my_label_parsed <- function (variable, value) {
  if (variable == "Model") {
    return(as.character(value))
  } else {
    lapply(as.character(value), function(x) parse(text = x))    
  }
}
```

``` r
plot_ref1 <-  macro_ref %>% 
  filter(region == "Striatum" | region == "Occipital cortex") %>% 
  filter(str_detect(Ligand, "SCH")) %>% 
  arrange(desc(Model)) %>% 
  mutate(Model = fct_inorder(Model)) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")

plot_ref2 <-  macro_ref %>% 
  filter(region == "Striatum" | region == "Occipital cortex") %>% 
  filter(str_detect(Ligand, "AZ")) %>% 
  arrange(desc(Model)) %>% 
  mutate(Model = fct_inorder(Model)) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")


plot_ref <- cowplot::plot_grid(plot_ref2, plot_ref1, nrow = 2 )

plot_ref
# ggsave(filename = "figs/macro_ref_higher.png", plot = plot_ref, width = 170, height = 180, units = "mm", dpi = 300)
```

``` r
plot_aif1 <-  macro_aif %>% 
  filter(region == "Thalamus") %>% 
  filter(str_detect(Ligand, "PBR28")) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free_y",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr") 

plot_aif2 <-  macro_aif %>% 
  filter(region == "Thalamus") %>% 
  filter(str_detect(Ligand, "PK")) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free_y",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")

plot_aif <- cowplot::plot_grid(plot_aif1, plot_aif2, nrow = 2 )

plot_aif

# ggsave(filename = "figs/macro_aif_higher.png", plot = plot_aif, width = 170, height = 180, units = "mm", dpi = 300)
```

``` r
plot_ref1_fc <-  macro_ref %>% 
    filter(region == "Frontal cortex") %>% 
  filter(str_detect(Ligand, "SCH")) %>% 
  arrange(desc(Model)) %>% 
  mutate(Model = fct_inorder(Model)) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")

plot_ref2_fc <-  macro_ref %>% 
  filter(region == "Frontal cortex") %>% 
  filter(str_detect(Ligand, "AZ")) %>% 
  arrange(desc(Model)) %>% 
  mutate(Model = fct_inorder(Model)) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")


plot_ref_fc <- cowplot::plot_grid(plot_ref2_fc, plot_ref1_fc, nrow = 2 )

plot_ref_fc


# ggsave(filename = "figs/macro_ref_lower.png", plot = plot_ref_fc, width = 170, height = 180, units = "mm", dpi = 300)
```

``` r
plot_aif1_fc <-  macro_aif %>% 
  filter(region == "Frontal cortex") %>% 
  filter(str_detect(Ligand, "PBR28")) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free_y",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr") 

plot_aif2_fc <-  macro_aif %>% 
 filter(region == "Frontal cortex") %>% 
  filter(str_detect(Ligand, "PK")) %>% 
  ggplot(aes(x = kinfitR, y = pmod, color = id)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=id), linetype="dotted") +
  facet_grid(cols = vars(Model), vars(Ligand), 
             scales = "free_y",
             labeller=my_label_parsed)+ theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")

plot_aif_fc <- cowplot::plot_grid(plot_aif1_fc, plot_aif2_fc, nrow = 2 )

plot_aif_fc


# ggsave(filename = "figs/macro_aif_lower.png", plot = plot_aif_fc, width = 170, height = 180, units = "mm", dpi = 300)
```

\#Spaghetti plot for the Frontal Cortex

``` r
spaghetti <-macro %>%
  filter(region == "FC") %>% 
  mutate(PETNo = as.integer(PETNo)) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD")) %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "KinfitR")) %>% 
  mutate(Ligand = factor(Ligand, levels = c("SCH23390","PBR28","PK11195","AZ10419369"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)"), ordered = TRUE))

spaghetti_pk <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>%
    filter(Model == "2TCM") %>% 
 ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid( Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(V[T])) + 
  theme(axis.title.x = element_blank())  

spaghetti_az <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>%
    filter(Model == "SRTM") %>% 
ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid(Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(BP[ND]))+ 
  theme(axis.title.x = element_blank()) 

spaghetti_pbr <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  filter(Model == "2TCM") %>%
 ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid(Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(V[T])) + 
  theme(axis.title.x = element_blank()) 

 

spaghetti_sch <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  filter(Model == "SRTM") %>% 
 ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid(Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(BP[ND])) + 
  theme(axis.title.x = element_blank()) 


spaghetti_plot_two <- plot_grid(spaghetti_sch,spaghetti_az, spaghetti_pk,spaghetti_pbr, nrow = 4, align = "hv")

x.grob <- textGrob("Scan number \n (in chronological order)", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
 
grid.arrange(arrangeGrob(spaghetti_plot_two, bottom = x.grob))

# spag <- grid.arrange(arrangeGrob(spaghetti_plot_two, bottom = x.grob))
# 
# ggsave(filename = "figs/spaghetti_fc.png", plot = spag, width = 170, height = 225, units = "mm", dpi = 300)
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

# kable(table_2, format = "latex", booktabs = T, digits = 3) %>%
#   kable_styling(latex_options = c("striped", "condensed", "scale_down"),
#                 full_width = F) 

# save_kable(table_2, 'figs/table_2.pdf')
```

### New combined tables 1 and 2

``` r
combo <- rbind(table_2, table_1) %>% 
mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28")) %>% 
  select(Ligand, Model, contains("Cor"), everything())
  

kable(combo, format = "html", booktabs = T,  digits = 2, escape = F,
             caption = "Agreement between KinfitR and PMOD",
      col.names = c("Ligand", "Model", "Region 1", "Region 2", "Region 1",
                    "Region 2", "Bias 1", "Bias 2"), align= 'c') %>%
kable_styling(bootstrap_options = c("condensed"),full_width = F) %>%
collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
pack_rows("Invasive", 1, 6) %>%
pack_rows("Non-Invasive", 7, 12) %>%
add_header_above(c(" " = 2, "ICC" = 2, "Pearson's r" = 2, "Bias (%)" = 2))
# 
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
kable_styling(full_width = F, bootstrap_options = "condensed") %>%
collapse_rows(columns = 1, latex_hline = "major", valign = "middle") %>%
pack_rows("Invasive", 1, 6) %>%
pack_rows("Non-Invasive", 7, 12) %>%
add_header_above(c(" " = 2, "Pearson's r" = 2, "ICC" = 2, "Bias (%)" = 2))
# 
# save_kable(mix, 'figs/table_1.html')
```

### Comparison of the models

This figure will be included in the supplementary material.

The first figure compares the reference tissue models

``` r
ref_mod <- ref %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "SCH23390"))%>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD")) %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "KinfitR")) %>% 
  mutate(Ligand = factor(Ligand, levels = c("SCH23390","AZ10419369"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)"), ordered = TRUE)) %>% spread(Model, Value) 

az_srtm_OC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, SRTM) %>% 
  select(- MRTM2, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "OC") %>% 
  mutate(Model = replicate(length(PET), "SRTM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "OC", 
                              replacement = "Occipital cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

az_srtm_FC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, SRTM) %>% 
  select(- MRTM2, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "SRTM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

az_mrtm2_OC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MRTM2) %>% 
  select(- SRTM, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "OC") %>% 
  mutate(Model = replicate(length(PET), "MRTM2")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "OC", 
                              replacement = "Occipital cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

az_mrtm2_FC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MRTM2) %>% 
  select(- SRTM, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "MRTM2")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

az_ref_logan_OC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "ref Logan") %>% 
  select(- MRTM2, -SRTM) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "OC") %>% 
  mutate(Model = replicate(length(PET), "ref Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "OC", 
                              replacement = "Occipital cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

az_ref_logan_FC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "ref Logan") %>% 
  select(- MRTM2, -SRTM) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "ref Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

sch_srtm_STR <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, SRTM) %>% 
  select(- MRTM2, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "STR") %>% 
  mutate(Model = replicate(length(PET), "SRTM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "STR", 
                              replacement = "Striatum")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

sch_srtm_FC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, SRTM) %>% 
  select(- MRTM2, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "SRTM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

sch_mrtm2_STR <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MRTM2) %>% 
  select(- SRTM, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "STR") %>% 
  mutate(Model = replicate(length(PET), "MRTM2")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "STR", 
                              replacement = "Striatum")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

sch_mrtm2_FC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MRTM2) %>% 
  select(- SRTM, -"ref Logan") %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "MRTM2")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

sch_ref_logan_STR <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "ref Logan") %>% 
  select(- MRTM2, -SRTM) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "STR") %>% 
  mutate(Model = replicate(length(PET), "ref Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "STR", 
                              replacement = "Striatum")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

sch_ref_logan_FC <- ref_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "ref Logan") %>% 
  select(- MRTM2, -SRTM) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "ref Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

srtm_plot <- plot_grid(az_srtm_OC,az_srtm_FC,sch_srtm_STR,sch_srtm_FC, ncol =1)

srtm_title <- ggdraw() + draw_label("SRTM", fontface='bold')

srtm_plot <- plot_grid(srtm_title, srtm_plot, ncol = 1, rel_heights=c(0.1, 1))

mrtm2_plot <- plot_grid(az_mrtm2_OC,az_mrtm2_FC,sch_mrtm2_STR,sch_mrtm2_FC, ncol =1)

mrtm2_title <- ggdraw() + draw_label("MRTM2", fontface='bold')

mrtm2_plot <- plot_grid(mrtm2_title, mrtm2_plot, ncol = 1, rel_heights=c(0.1, 1))

ref_plot <- plot_grid(az_ref_logan_OC,az_ref_logan_FC,sch_ref_logan_STR,
                      sch_ref_logan_FC, ncol =1)

ref_title <- ggdraw() + draw_label("Non-invasive Logan", fontface='bold')

ref_plot <- plot_grid(ref_title, ref_plot, ncol = 1, rel_heights=c(0.1, 1))

combined_plot <- plot_grid(srtm_plot, mrtm2_plot, ref_plot, ncol = 3)

title <- ggdraw() + draw_label("Relationship between the same model \nfor PMOD compared to kinfitr", fontface='bold')

plot_grid(title, combined_plot , ncol = 1, rel_heights=c(0.1, 1))  
```

The next figure shows the relationship between kinfitr and pmod for the
models that require arterial blood sampling

``` r
aif_mod <- aif %>% 
  unnest() %>% 
  ungroup() %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD")) %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "KinfitR")) %>% 
  mutate(Ligand = factor(Ligand, levels = c("PBR28","PK11195"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)"), ordered = TRUE)) %>% spread(Model, Value)

pk_2tcm_tha <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "2TCM") %>% 
  select(- Logan, -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "THA") %>% 
  mutate(Model = replicate(length(PET), "2TCM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pk_logan_tha <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, Logan) %>% 
  select(- "2TCM", -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "THA") %>% 
  mutate(Model = replicate(length(PET), "Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pk_ma1_tha <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MA1) %>% 
  select(- "2TCM", -Logan) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "THA") %>% 
  mutate(Model = replicate(length(PET), "MA1")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pk_2tcm_fc <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "2TCM") %>% 
  select(- Logan, -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "2TCM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal Cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pk_logan_fc <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, Logan) %>% 
  select(- "2TCM", -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal Cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pk_ma1_fc <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MA1) %>% 
  select(- "2TCM", -Logan) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "MA1")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal Cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pbr_2tcm_tha <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "2TCM") %>% 
  select(- Logan, -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "THA") %>% 
  mutate(Model = replicate(length(PET), "2TCM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pbr_logan_tha <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, Logan) %>% 
  select(- "2TCM", -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "THA") %>% 
  mutate(Model = replicate(length(PET), "Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pbr_ma1_tha <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MA1) %>% 
  select(- "2TCM", -Logan) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "THA") %>% 
  mutate(Model = replicate(length(PET), "MA1")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "THA", 
                              replacement = "Thalamus")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pbr_2tcm_fc <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, "2TCM") %>% 
  select(- Logan, -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "2TCM")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal Cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pbr_logan_fc <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, Logan) %>% 
  select(- "2TCM", -MA1) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "Logan")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal Cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

pbr_ma1_fc <- aif_mod %>% 
  ungroup() %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  group_by(Ligand, region, PET) %>% 
  spread(software, MA1) %>% 
  select(- "2TCM", -Logan) %>% 
  summarize_at(vars(KinfitR, PMOD), funs(toString(unique(.[!is.na(.)])))) %>% 
  mutate(kinfitr = as.double(KinfitR)) %>% 
  mutate(PMOD = as.double(PMOD)) %>% 
  filter(region == "FC") %>% 
  mutate(Model = replicate(length(PET), "MA1")) %>%
  ungroup() %>% 
   mutate(region = str_replace(string = region, pattern = "FC", 
                              replacement = "Frontal Cortex")) %>%
  ggplot(aes(x = kinfitr, y = PMOD)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red")+
  facet_grid(. ~ Ligand + region, labeller = labeller(Ligand = label_parsed))

two_tcm_plot <- plot_grid(pk_2tcm_tha,pk_2tcm_fc,pbr_2tcm_tha, pbr_2tcm_fc, ncol =1)

two_tcm_title <- ggdraw() + draw_label("2TCM", fontface='bold')

two_tcm_plot <- plot_grid(two_tcm_title, two_tcm_plot, ncol = 1, rel_heights=c(0.1, 1))


ma1_plot <- plot_grid(pk_ma1_tha,pk_ma1_fc,pbr_ma1_tha, pbr_ma1_fc, ncol =1)

ma1_title <- ggdraw() + draw_label("MA1", fontface='bold')

ma1_plot <- plot_grid(ma1_title, ma1_plot, ncol = 1, rel_heights=c(0.1, 1))

logan_plot <- plot_grid(pk_logan_tha,pk_logan_fc,pbr_logan_tha, pbr_logan_fc, ncol =1)

logan_title <- ggdraw() + draw_label("Invasive Logan", fontface='bold')

logan_plot <- plot_grid(logan_title, logan_plot, ncol = 1, rel_heights=c(0.1, 1))

combine_plot <- plot_grid(two_tcm_plot, ma1_plot, logan_plot, ncol = 3)

titles <- ggdraw() + draw_label("Relationship between the same model \nfor PMOD compared to kinfitr", fontface='bold')

plot_grid(titles, combine_plot , ncol = 1, rel_heights=c(0.1, 1))
```

Table 3 - consistency test\_retest
==================================

``` r
macro_sch <- macro %>% 
  filter(Ligand == "sch") %>% 
  filter(region == "STR")

macro_az <- macro %>% 
  filter(Ligand == "az") %>% 
  filter(region == "OC")

macro_aif <- macro %>% 
  filter(measure == "Vt") %>% 
  filter(region == "THA")

macro_sch_2 <- macro %>% 
  filter(Ligand == "sch") %>% 
  filter(region == "FC")

macro_az_2 <- macro %>% 
  filter(Ligand == "az") %>% 
  filter(region == "FC")

macro_aif_2 <- macro %>% 
  filter(measure == "Vt") %>% 
  filter(region == "FC")

  
macro <- bind_rows(macro_sch, macro_aif, macro_az)

macro_2 <- bind_rows(macro_sch_2, macro_aif_2, macro_az_2)


trt <- macro %>% 
  select(software, Ligand, Model, Value, subjname, PETNo ) %>%
  group_by(Ligand, software, Model) %>% 
  nest(.key = "data")

trt <- trt %>% 
  group_by(Ligand, software, Model) %>%
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "subjname",
                                       rater = "PETNo")),
         trt_tidy = map(trt, c("tidy")))

#Note: multiplied "VAR" by 100 so that it is a percentage

trt <- trt %>% 
  select(trt_tidy) %>% 
  unnest(trt_tidy) %>% 
  select(Ligand, Software = software, Model, Mean = mean, "CV" = cv, ICC = icc, "WSCV" = wscv, "AV" = absvar) %>% 
  ungroup() %>% 
  mutate(`AV` = `AV` *100) %>% 
  mutate(`WSCV` = `WSCV` * 100) %>% 
  mutate(`CV` = `CV` * 100) %>% 
  mutate(Software = str_replace( string = Software, pattern = "pmod", replacement = "PMOD")) %>% 
  mutate(Software = str_replace( string = Software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>%
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "ASCH23390")) %>% 
   mutate(Model = str_replace(string = Model, pattern = "SRTM", replacement = "aSRTM")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "ref Logan", replacement = "bref Logan")) %>%
  arrange(desc(Ligand), Model) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "ASCH23390", replacement = "SCH23390")) %>%
  mutate(Model = str_replace(string = Model, pattern = "aSRTM", replacement = "SRTM")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "bref Logan", replacement = "ref Logan")) %>%
mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28"))

kable(trt, format = "html", booktabs = TRUE,digits = c(0,0,0, 2, 1,2,1,1),
      escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD",  col.names = c("Ligand", "Software", "Model", "Mean", "CV (%)", "ICC", "WSCV (%)",
               "AV (%)") , align = 'c') %>%
kable_styling(bootstrap_options = c("condensed"), full_width = F) %>%
collapse_rows(columns = 1, latex_hline = "major" ,valign = "middle") %>%
row_spec(2, extra_latex_after = "\\cline{2-8}") %>%
row_spec(4, extra_latex_after = "\\cline{2-8}") %>%
row_spec(8, extra_latex_after = "\\cline{2-8}") %>%
row_spec(10, extra_latex_after = "\\cline{2-8}") %>%
row_spec(14, extra_latex_after = "\\cline{2-8}") %>%
row_spec(16, extra_latex_after = "\\cline{2-8}") %>%
row_spec(20, extra_latex_after = "\\cline{2-8}") %>%
row_spec(22, extra_latex_after = "\\cline{2-8}") %>%
pack_rows("Invasive", 1, 12) %>%
pack_rows("Non-Invasive", 13, 24)

# kable(trt, format = "latex", booktabs = TRUE,digits = c(0,0,0, 2, 1,2,1,1),
#       escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD", col.names = c("Ligand", "Software", "Model", "Mean", "CV (\\%)", "ICC", "WSCV (\\%)",
#                "AV (\\%)") ,  align = 'c') %>%
# kable_styling(full_width = F, latex_options = "scale_down") %>%
# collapse_rows(columns = 1, latex_hline = "major" ,valign = "middle") %>%
# row_spec(2, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(4, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(8, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(10, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(14, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(16, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(20, extra_latex_after = "\\cline{2-8}") %>%
# row_spec(22, extra_latex_after = "\\cline{2-8}") %>%
# pack_rows("Invasive", 1, 12) %>%
# pack_rows("Non-Invasive", 13, 24) %>%
# save_kable(., file = "figs/table_2.pdf")
```

``` r
trt_aif <- trt %>% 
  filter(str_detect(Ligand, "PBR|PK11195")) %>% 
  select(Ligand, Model, Software, everything())

trt_ref <- trt %>% 
  filter(str_detect(Ligand, "SCH|AZ")) %>% 
  select(Ligand, Model, Software, everything())
```

``` r
kable(trt_aif, format = "html", booktabs = TRUE, digits = c(0,0,0, 2, 1,2,1,1),
      escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD for the invasive models",  col.names = c("Ligand", "Software", "Model", "Mean", "CV (%)", "ICC", "WSCV (%)",
               "AV (%)") , align = 'c') %>%
kable_styling(bootstrap_options = c("condensed"), full_width = F) %>%
collapse_rows(columns = c(1,2), latex_hline = "major" ,valign = "middle") %>% 
  save_kable(. , "figs/trt_aif.html")

kable(trt_ref, format = "html", booktabs = TRUE, digits = c(0,0,0, 2, 1,2,1,1),
      escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD for the non-invasive models",  col.names = c("Ligand", "Software", "Model", "Mean", "CV (%)", "ICC", "WSCV (%)",
               "AV (%)") , align = 'c') %>%
kable_styling(bootstrap_options = c("condensed"), full_width = F) %>%
collapse_rows(columns = c(1,2), latex_hline = "major" ,valign = "middle") %>% 
  save_kable(. , "figs/trt_ref.html")
```

``` r
trt_2 <- macro_2 %>% 
  select(software, Ligand, Model, Value, subjname, PETNo ) %>%
  group_by(Ligand, software, Model) %>% 
  nest(.key = "data")

trt_2 <- trt_2 %>% 
  group_by(Ligand, software, Model) %>%
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "subjname",
                                       rater = "PETNo")),
         trt_tidy = map(trt, c("tidy")))

#Note: multiplied "VAR" by 100 so that it is a percentage

trt_2 <- trt_2 %>% 
  select(trt_tidy) %>% 
  unnest(trt_tidy) %>% 
  select(Ligand, Software = software, Model, Mean = mean, "CV" = cv, ICC = icc, "WSCV" = wscv, "AV" = absvar) %>% 
  ungroup() %>% 
  mutate(`AV` = `AV` *100) %>% 
  mutate(`WSCV` = `WSCV` * 100) %>% 
  mutate(`CV` = `CV` * 100) %>% 
  mutate(Software = str_replace( string = Software, pattern = "pmod", replacement = "PMOD")) %>% 
  mutate(Software = str_replace( string = Software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>%
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "ASCH23390")) %>% 
   mutate(Model = str_replace(string = Model, pattern = "SRTM", replacement = "aSRTM")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "ref Logan", replacement = "bref Logan")) %>%
  arrange(desc(Ligand), Model) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "ASCH23390", replacement = "SCH23390")) %>%
  mutate(Model = str_replace(string = Model, pattern = "aSRTM", replacement = "SRTM")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "bref Logan", replacement = "ref Logan")) %>%
mutate(Ligand = str_replace(string = Ligand, pattern = "AZ10419369", replacement = "[<sup>11</sup>C]AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "SCH23390", replacement = "[<sup>11</sup>C]SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PK11195", replacement = "[<sup>11</sup>C]PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "PBR28", replacement = "[<sup>11</sup>C]PBR28"))

kable(trt_2, format = "html", booktabs = TRUE,digits = c(0,0,0, 2, 1,2,1,1),
      escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD",  col.names = c("Ligand", "Software", "Model", "Mean", "CV (%)", "ICC", "WSCV (%)",
               "AV (%)") , align = 'c') %>%
kable_styling(bootstrap_options = c("condensed"), full_width = F) %>%
collapse_rows(columns = 1, latex_hline = "major" ,valign = "middle") %>%
row_spec(2, extra_latex_after = "\\cline{2-8}") %>%
row_spec(4, extra_latex_after = "\\cline{2-8}") %>%
row_spec(8, extra_latex_after = "\\cline{2-8}") %>%
row_spec(10, extra_latex_after = "\\cline{2-8}") %>%
row_spec(14, extra_latex_after = "\\cline{2-8}") %>%
row_spec(16, extra_latex_after = "\\cline{2-8}") %>%
row_spec(20, extra_latex_after = "\\cline{2-8}") %>%
row_spec(22, extra_latex_after = "\\cline{2-8}") %>%
pack_rows("Invasive", 1, 12) %>%
pack_rows("Non-Invasive", 13, 24)
```

``` r
trt_aif_2 <- trt_2 %>% 
  filter(str_detect(Ligand, "PBR|PK11195")) %>% 
  select(Ligand, Model, Software, everything())

trt_ref_2 <- trt_2 %>% 
  filter(str_detect(Ligand, "SCH|AZ")) %>% 
  select(Ligand, Model, Software, everything())
```

``` r
kable(trt_aif_2, format = "html", booktabs = TRUE, digits = c(0,0,0, 2, 1,2,1,1),
      escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD for the invasive models",  col.names = c("Ligand", "Software", "Model", "Mean", "CV (%)", "ICC", "WSCV (%)",
               "AV (%)") , align = 'c') %>%
kable_styling(bootstrap_options = c("condensed"), full_width = F) %>%
collapse_rows(columns = c(1,2), latex_hline = "major" ,valign = "middle") %>% 
  save_kable(. , "figs/trt_aif_fc.html")

kable(trt_ref_2, format = "html", booktabs = TRUE, digits = c(0,0,0, 2, 1,2,1,1),
      escape = FALSE, caption = "Test-retest reliability of kinfitr and PMOD for the non-invasive models",  col.names = c("Ligand", "Software", "Model", "Mean", "CV (%)", "ICC", "WSCV (%)",
               "AV (%)") , align = 'c') %>%
kable_styling(bootstrap_options = c("condensed"), full_width = F) %>%
collapse_rows(columns = c(1,2), latex_hline = "major" ,valign = "middle") %>% 
  save_kable(. , "figs/trt_ref_fc.html")
```

Spaghetti plot - The occipital cortex
=====================================

``` r
spaghetti <-macro %>%
  filter(region != "FC") %>%
  mutate(PETNo = as.integer(PETNo)) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "SCH23390"))%>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pk", replacement = "PK11195")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "pbr28", replacement = "PBR28")) %>% 
  mutate(software = str_replace(string = software, pattern = "pmod", replacement = "PMOD")) %>% 
  mutate(software = str_replace(string = software, pattern = "kinfitR", replacement = "kinfitr")) %>% 
  mutate(Ligand = factor(Ligand, levels = c("SCH23390","PBR28","PK11195","AZ10419369"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)"), ordered = TRUE))

spaghetti_pk <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)") %>%
    filter(Model == "2TCM") %>% 
 ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid( Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(V[T])) + 
  theme(axis.title.x = element_blank())  

spaghetti_az <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)") %>%
    filter(Model == "SRTM") %>% 
ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid(Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(BP[ND]))+ 
  theme(axis.title.x = element_blank()) 

spaghetti_pbr <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)") %>% 
  filter(Model == "2TCM") %>%
 ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid(Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(V[T])) + 
  theme(axis.title.x = element_blank()) 

 

spaghetti_sch <- spaghetti %>% 
  filter(Ligand == "'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)") %>% 
  filter(Model == "SRTM") %>% 
 ggplot(aes(x = PETNo, y = Value, 
            group = subjname, colour=subjname)) +
  geom_point(size = 3.3) +
  geom_smooth(method = 'lm', se = FALSE) + 
  scale_x_continuous(limits = c(0.75, 2.25), breaks = c(1,2)) +
  facet_grid(Ligand ~ software, scales = "free", labeller=label_parsed) + 
  theme(legend.position = "none")+
  labs(y = expression(BP[ND])) + 
  theme(axis.title.x = element_blank()) 


spaghetti_plot_two <- plot_grid(spaghetti_sch,spaghetti_az, spaghetti_pk,spaghetti_pbr, nrow = 4, align = "hv")

x.grob <- textGrob("Scan number \n (in chronological order)", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))
 
grid.arrange(arrangeGrob(spaghetti_plot_two, bottom = x.grob))

spag <- grid.arrange(arrangeGrob(spaghetti_plot_two, bottom = x.grob))

# ggsave(filename = "spaghetti.pdf", plot = spag, width = 170, height = 225, units = "mm", dpi = 300)

# ggsave(filename = "figs/spaghetti_higher.png", plot = spag, width = 170, height = 225, units = "mm", dpi = 300)
```

Microparameters
===============

Load microparameters
--------------------

``` r
micro_aif_pmod <- readRDS('Data_PMOD_RDS/micro_aif_pmod.rds') %>% 
  mutate(PETNo = as.numeric(PETNo))%>% 
  mutate(subjnumber = as.character(id))%>%
  rename(K1_pmod = K1, k2_pmod = k2, k3_pmod = k3, k4_pmod = k4) %>% 
  select(- software)

micro_ref_pmod <- readRDS('Data_PMOD_RDS/micro_ref_pmod.rds') %>% 
  mutate(PETNo = as.numeric(PETNo))%>% 
  mutate(subjnumber = as.character(id))%>%
  rename(R1_pmod = R1, k2_pmod = k2) %>% 
  select(-software)

micro_aif_kinfitr <- readRDS('Data_PMOD_RDS/micro_aif_kinfitr.rds')%>% 
  mutate(subjnumber = as.character(id))%>%
  rename(K1_kinfitr = K1, k2_kinfitr = k2, k3_kinfitr = k3, k4_kinfitr = k4) %>% 
  select(- software)

micro_ref_kinfitr <- readRDS('Data_PMOD_RDS/micro_ref_kinfitr.rds') %>%
  mutate(subjnumber = as.character(id))%>%
  rename(R1_kinfitr = R1, k2_kinfitr = k2) %>% 
  select(- software)
```

Combine datasets
----------------

``` r
micro_ref <- left_join(micro_ref_pmod, micro_ref_kinfitr)%>% 
  arrange(Ligand) %>% 
  mutate(Ligand = factor(Ligand, levels = c("AZ10419369","SCH23390"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(AZ10419369)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(SCH23390)"), ordered = TRUE)) 

micro_aif <- left_join(micro_aif_pmod, micro_aif_kinfitr) %>% 
  arrange(Ligand) %>% 
  mutate(Ligand = factor(Ligand, levels = c("PBR28","PK11195"), labels = c("'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PBR28)","'[' ~ {}^11 ~ 'C' ~ ']' ~ plain(PK11195)"), ordered = TRUE)) %>% 
  filter(subjname != "rwrd") 
```

Fig 2, microparameter dotplot with abline for srtm
--------------------------------------------------

Using the Occipital cortex ROI

``` r
R1 <- micro_ref %>% 
  ggplot(aes(x = R1_kinfitr, y = R1_pmod, color = subjnumber)) + 
   geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=subjnumber), linetype="dotted") +
  facet_wrap( ~ Ligand, scales = "free", labeller=label_parsed)+ theme(legend.position = "none")+
  ggtitle("R1 (unitless)")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title = element_blank())

k2_ref <- micro_ref %>% 
  ggplot(aes(x = k2_kinfitr, y = k2_pmod, color = subjnumber)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_line(aes(group=subjnumber), linetype="dotted") +
  facet_wrap( ~ Ligand, scales = "free", labeller=label_parsed)+ theme(legend.position = "none") +
  ggtitle(expression(paste (k[2], "(", min^-1, ")", sep = ' ')))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_blank())

plot_ref <- plot_grid(R1, k2_ref, nrow = 2, align = "v")

y.grob <- textGrob("PMOD", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot=90)

x.grob <- textGrob("kinfitr", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))

grid.arrange(arrangeGrob(plot_ref, left = y.grob, bottom = x.grob))

# micro_reference <- grid.arrange(arrangeGrob(plot_ref, left = y.grob, bottom = x.grob))
# 
# ggsave(filename = "micro_ref.pdf", plot = micro_reference, width = 170, height = 180, units = "mm", dpi = 300)

# ggsave(filename = "micro_ref.png", plot = micro_reference, width = 170, height = 180, units = "mm", dpi = 300)
```

Fig 3, microparameter dotplot with abline for 2tcm
--------------------------------------------------

Using the Thalamus ROI

``` r
K1 <- micro_aif %>% 
  ggplot(aes(x = K1_kinfitr, y = K1_pmod, color = subjnumber)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=subjnumber), linetype="dotted") +
  facet_wrap( ~ Ligand, scales = "free", labeller=label_parsed)+ theme(legend.position = "none")+
  ggtitle(expression(paste (K[1], "(",mL[plasma], min^-1, mL[tissue]^-1, ")", sep = ' ')))+
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust=1))

k2_aif <- micro_aif %>% 
  ggplot(aes(x = k2_kinfitr, y = k2_pmod, color = subjnumber)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_line(aes(group=subjnumber), linetype="dotted") +
  facet_wrap( ~ Ligand, scales = "free", labeller=label_parsed)+ theme(legend.position = "none") +
  ggtitle(expression(paste (k[2], "(", min^-1, ")", sep = ' ')))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust=1))

k3 <- micro_aif %>% 
  ggplot(aes(x = k3_kinfitr, y = k3_pmod, color = subjnumber)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_line(aes(group=subjnumber), linetype="dotted") +
  facet_wrap( ~ Ligand, scales = "free", labeller=label_parsed)+ theme(legend.position = "none") +
  ggtitle(expression(paste (k[3], "(", min^-1, ")", sep = ' ')))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust=1))


k4 <- micro_aif %>% 
  ggplot(aes(x = k4_kinfitr, y = k4_pmod, color = subjnumber)) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  geom_line(aes(group=subjnumber), linetype="dotted") +
  facet_wrap( ~ Ligand, scales = "free", labeller=label_parsed)+ theme(legend.position = "none") +
  ggtitle(expression(paste (k[4], "(",  min^-1 , ")", sep = ' ')))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust=1))

plot_aif <- plot_grid(K1, k2_aif, k3, k4)

y.grob <- textGrob("PMOD", 
                   gp=gpar(fontface="bold", col="black", fontsize=12), rot=90)

x.grob <- textGrob("kinfitr", 
                   gp=gpar(fontface="bold", col="black", fontsize=12))

grid.arrange(arrangeGrob(plot_aif, left = y.grob, bottom = x.grob))

# micro_art <- grid.arrange(arrangeGrob(plot_aif, left = y.grob, bottom = x.grob))
 
# # # ggsave(filename = "micro_aif.pdf", plot = micro_art, width = 170, height = 140, units = "mm", dpi = 300)
# # 
# ggsave(filename = "micro_aif.png", plot = micro_art, width = 170, height = 140, units = "mm", dpi = 300)
```

correlations between the microparameters

``` r
get_cors <- function(vals) {
  vals %>% 
    correlate(quiet = T) %>% 
    shave() %>% 
    stretch() %>% 
    drop_na()
}

k1_aif_cor <- micro_aif %>% 
  select(Ligand, PET, contains("K1")) %>% 
  select(-PET) %>% 
  mutate(Ligand = ifelse(str_detect(Ligand,"PBR"), "PBR28", "PK11195")) %>% 
  group_by(Ligand) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  select(r) %>% 
  ungroup() %>% 
  summarize(mean = mean(r)) %>% 
  mutate(microparameter = "K1")

k2_aif_cor <- micro_aif %>% 
  select(Ligand, PET, contains("k2")) %>% 
  select(-PET) %>% 
  mutate(Ligand = ifelse(str_detect(Ligand,"PBR"), "PBR28", "PK11195")) %>% 
  group_by(Ligand) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  select(r) %>% 
  ungroup() %>% 
  summarize(mean = mean(r)) %>% 
  mutate(microparameter = "k2")

k3_aif_cor <- micro_aif %>% 
  select(Ligand, PET, contains("k3")) %>% 
  select(-PET) %>% 
  mutate(Ligand = ifelse(str_detect(Ligand,"PBR"), "PBR28", "PK11195")) %>% 
  group_by(Ligand) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  select(r) %>% 
  ungroup() %>% 
  summarize(mean = mean(r)) %>% 
  mutate(microparameter = "k3")


k4_aif_cor <- micro_aif %>% 
  select(Ligand, PET, contains("k4")) %>% 
  select(-PET) %>% 
  mutate(Ligand = ifelse(str_detect(Ligand,"PBR"), "PBR28", "PK11195")) %>% 
  group_by(Ligand) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  select(r) %>% 
  ungroup() %>% 
  summarize(mean = mean(r)) %>% 
  mutate(microparameter = "k4")

bind_rows(k1_aif_cor, k2_aif_cor, k3_aif_cor, k4_aif_cor) %>% 
  kable(digits = 2, caption = "mean correlations between the microparameters calculated by each tool for the invasive models")
```

correlations between the tools for the microparameters of SRTM

``` r
r1_ref_cor <- micro_ref %>% 
  select(Ligand, PET, contains("R1")) %>% 
  select(-PET) %>% 
  mutate(Ligand = ifelse(str_detect(Ligand,"AZ"), "AZ10419369", 
                         "SCH23390")) %>% 
  group_by(Ligand) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  select(r) %>% 
  ungroup() %>% 
  summarize(mean = mean(r)) %>% 
  mutate(microparameter = "R1")

k2_ref_cor <- micro_ref %>% 
  select(Ligand, PET, contains("k2")) %>% 
  select(-PET) %>% 
  mutate(Ligand = ifelse(str_detect(Ligand,"AZ"), "AZ10419369",
                         "SCH23390")) %>% 
  group_by(Ligand) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  select(r) %>% 
  ungroup() %>% 
  summarize(mean = mean(r)) %>% 
  mutate(microparameter = "k2")

bind_rows(r1_ref_cor, k2_ref_cor) %>% 
  kable(digits = 2, caption = "mean correlations between the microparameters calculated by each tool for the non-invasive models")
```

Correlation plots showing the relationship between the different models
for the same Ligand.

``` r
macro <- bind_rows(kinfitr_macro, pmod_macro) %>% 
  rename(Ligand = tracer, Model = model) %>% 
  mutate(Model = str_replace(string = Model, pattern = "mrtm2", replacement = "MRTM2")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "srtm", replacement = "SRTM")) %>% mutate(Model = str_replace(string = Model, pattern = "ma1", replacement = "MA1")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "two_tcm", replacement = "2TCM")) %>% mutate(Model = str_replace(string = Model, pattern = "ref_logan", replacement = "ref Logan")) %>% 
  mutate(Model = str_replace(string = Model, pattern = "logan", replacement = "Logan")) 

col2 <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))

par(mfrow=c(2,4))

cor_az_oc_kinfitR <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "az") %>% 
  filter(region == "OC") %>% 
  filter(software == "kinfitR") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(AZ10419369 ~ OC ~ kinfitr),
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
                 number.digits = 2, title=expression(AZ10419369 ~ FC  ~ kinfitr),
                 mar=c(0,0,1,0))

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
                 number.digits = 2, title=expression(SCH23390 ~ STR ~ kinfitr),
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
                 number.digits = 2, title=expression(SCH23390 ~ FC ~ kinfitr),
                 mar=c(0,0,1,0))


cor_az_oc_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>%
  filter(Ligand == "az") %>%
  filter(region == "OC") %>% 
  filter(software == "pmod") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(AZ10419369 ~ OC ~ PMOD),
                 mar=c(0,0,1,0))

cor_az_fc_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "az") %>%
  filter(region == "FC") %>% 
  filter(software == "pmod") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(AZ10419369 ~ FC  ~ PMOD),
                 mar=c(0,0,1,0))



cor_sch_str_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>%  
  spread(Model, Value) %>% 
  filter(Ligand == "sch") %>%
  filter(region == "STR") %>% 
  filter(software == "pmod") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(SCH23390 ~ STR ~ PMOD),
                 mar=c(0,0,1,0))

cor_sch_fc_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>%  
  spread(Model, Value) %>% 
  filter(Ligand == "sch") %>%
  filter(region == "FC") %>% 
  filter(software == "pmod") %>% 
  select(MRTM2:SRTM) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(SCH23390 ~ FC ~ PMOD),
                 mar=c(0,0,1,0))
```

``` r
par(mfrow=c(2,4))

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
                 number.digits = 2, title=expression(PK11195 ~ THA ~ kinfitr),
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
                 number.digits = 2, title=expression(PK11195 ~ FC ~ kinfitr),
                 mar=c(0,0,1,0))

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
                 number.digits = 2, title=expression(PBR28 ~ THA ~ kinfitr),
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
                 number.digits = 2, title=expression(PBR28 ~ FC ~ kinfitr),
                 mar=c(0,0,1,0))


cor_pk_tha_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pk") %>% 
  filter(region == "THA") %>% 
  filter(software == "pmod") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PK11195 ~ THA ~ PMOD),
                 mar=c(0,0,1,0))

cor_pk_fc_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pk") %>% 
  filter(region == "FC") %>% 
  filter(software == "pmod") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PK11195 ~ FC ~ PMOD),
                 mar=c(0,0,1,0))



cor_pbr28_tha_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pbr28") %>% 
  filter(region == "THA") %>% 
  filter(software == "pmod") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PBR28 ~ THA ~ PMOD),
                 mar=c(0,0,1,0))

cor_pbr28_fc_pmod <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pbr28") %>% 
  filter(region == "FC") %>% 
  filter(software == "pmod") %>% 
  select("2TCM":MA1) %>% 
  cor()%>%
  corrplot.mixed(lower='ellipse', upper='number', 
                 lower.col = col2(200), upper.col = col2(200),  diag='n',
                 number.digits = 2, title=expression(PBR28 ~ FC ~ PMOD),
                 mar=c(0,0,1,0))
```

\#median of the correlations

``` r
get_cors <- function(vals) {
  vals %>% 
    correlate(quiet = T) %>% 
    shave() %>% 
    stretch() %>% 
    drop_na()
}

invasive_median_r <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "pbr28" | Ligand == "pk") %>% 
  select(Ligand:MA1) %>% 
  select(-PET) %>% 
  group_by(Ligand, region, software) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  group_by(software) %>% 
  summarize(median_r = median(r),
            min = min(r),
            max = max(r),
            IQR = IQR(r))

invasive_median_r %>% 
  kable(digits = 2, caption = "correlations for the invasive models")

noninvasive_median_r <- macro %>% 
  select(Ligand, PET, region, Model, Value, software) %>% 
  spread(Model, Value) %>% 
  filter(Ligand == "az" | Ligand == "sch") %>% 
  select(Ligand:software, MRTM2:SRTM) %>% 
  select(-PET) %>% 
  group_by(Ligand, region, software) %>% 
  nest() %>% 
  mutate(cors = map(data, get_cors)) %>% 
  select(-data) %>% 
  unnest(cors) %>% 
  group_by(software) %>% 
  summarize(median_r = median(r),
            min = min(r),
            max = max(r),
            IQR = IQR(r))


noninvasive_median_r %>% 
  kable(digits = 2, caption = "correlations for the non-invasive models")
```

\#microparameter correlations

``` r
micro_cor_pk <- micro_aif %>% 
  filter(str_detect(Ligand, 'PK')) %>% 
  select(K1_pmod, K1_kinfitr,k2_pmod, k2_kinfitr,k3_pmod, k3_kinfitr,k4_pmod, k4_kinfitr) %>% 
  correlate() %>% 
  shave() %>% 
  stretch() %>% 
  filter(!is.na(r))  %>% 
  slice(-c(2:13,15:22,24:27)) %>% 
  select("Microparameter" = y,r) %>% 
  mutate(Microparameter = str_replace_all(string = Microparameter, pattern = "_kinfitr", replacement = " ")) 

micro_aif %>% 
  filter(str_detect(Ligand, 'PBR')) %>% 
  select(K1_pmod, K1_kinfitr,k2_pmod, k2_kinfitr,k3_pmod, k3_kinfitr,k4_pmod, k4_kinfitr) %>% 
  correlate() %>% 
  shave() %>% 
  stretch() %>% 
  filter(!is.na(r))  %>% 
  slice(-c(2:13,15:22,24:27)) %>% 
  select("Microparameter" = y,r) %>% 
  mutate(Microparameter = str_replace_all(string = Microparameter, pattern = "_kinfitr", replacement = " ")) %>% 
  mutate(r_pk = micro_cor_pk$r) %>% 
  rowwise() %>% mutate(r=mean(c(r, r_pk))) %>% 
  select(- r_pk) %>% 
  kable(digits = 2, caption = "correlation of microparameters between kinfitr and PMOD for the invasive models" )

micro_cor_az <- micro_ref %>% 
  filter(str_detect(Ligand, 'AZ')) %>% 
  select(R1_pmod,R1_kinfitr,k2_pmod, k2_kinfitr) %>% 
  correlate() %>% 
  shave() %>% 
  stretch() %>% 
  filter(!is.na(r)) %>% 
  filter(row_number()==1 | row_number()==n()) %>% 
  select("Microparameter" = y,r) %>% 
  mutate(Microparameter = str_replace_all(string = Microparameter, pattern = "_kinfitr", replacement = " ")) 

micro_ref %>% 
  filter(str_detect(Ligand, 'SCH')) %>% 
  select(R1_pmod,R1_kinfitr,k2_pmod, k2_kinfitr) %>% 
  correlate() %>% 
  shave() %>% 
  stretch() %>% 
  filter(!is.na(r)) %>% 
  filter(row_number()==1 | row_number()==n()) %>% 
  select("Microparameter" = y,r) %>% 
  mutate(Microparameter = str_replace_all(string = Microparameter, pattern = "_kinfitr", replacement = " ")) %>% 
  mutate(r_az = micro_cor_az$r) %>% 
  rowwise() %>% mutate(r=mean(c(r, r_az))) %>% 
  select(- r_az) %>% 
  kable(digits = 2, caption = "correlation of  microparameters between kinfitr and PMOD for the non-invasive models" )
```

New MRTM2
=========

``` r
sch_mrtm2 <- readRDS('RawData_PMOD_RDS/raw_kinfit_SCH23390_mrtm2.rds') %>% 
  mutate(tracer = replicate(length(PET), "sch")) %>%
   mutate(software = replicate(length(PET), "kinfitR"))  %>% 
  rename(region = Region, subjname = Subjname, mrtm2 = bp_MRTM2 )

az_mrtm2 <- readRDS('RawData_PMOD_RDS/raw_kinfit_AZ10419369_mrtm2.rds') %>% 
  mutate(tracer = replicate(length(PET), "az")) %>%
   mutate(software = replicate(length(PET), "kinfitR")) %>% 
  rename(region = Region, subjname = Subjname, mrtm2 = bp_MRTM2 )

mrtm2 <- pmod_macro %>% 
  select(-c(id, measure)) %>% 
  filter(model == "mrtm2") %>% 
  spread(key = model, value = Value) %>% 
  bind_rows(sch_mrtm2, az_mrtm2) %>% 
  mutate(Model = replicate(length(PET), "MRTM2")) %>% 
  rename(Value = mrtm2, Ligand = tracer)
```

``` r
#preparing the data for trt. 

Corr_mrtm2 <- mrtm2 %>%
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
 
Bias_mrtm2 <- mrtm2 %>%
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

mrtm2_ref <- mrtm2 %>%
  select(Ligand, region, Model, PET, software, Value) %>%
  group_by(Ligand, Model, region) %>%
  nest(.key = "data")

#Using trt from the relfeas package to obtain ICC values
#Note: the values were not arranged by the "subjectname+PETNo" or "PET" variable
#Therefore, I had to use the "rater" argument of trt to specify "software", as 
#the default behaviour is to take the value under as the second rater.

mrtm2_table <- mrtm2_ref %>% 
  group_by(Ligand, Model, region) %>% 
  mutate(trt = map(data, ~relfeas::trt(.x, 
                                       values = "Value", 
                                       cases = "PET",
                                       rater = "software")),
         trt_tidy = map(trt, c("tidy"))) %>% 
  select(trt_tidy) %>%  
  unnest() 

#Renaming OC and STR "Region 1" and making FC "Region 2"

mrtm2_table <- mrtm2_table %>%
  select(Ligand, Model, region, icc) %>%
  ungroup() %>% 
  mutate(region = case_when(
    region == "OC" ~ "ICC (Region 1)",
    region == "STR" ~ "ICC (Region 1)",
    TRUE ~ "ICC (Region 2)"
  ))  
 

#Rename the Ligands and then use Spread() on the regions and the ICC

mrtm2_table <- mrtm2_table %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "az", replacement = "AZ10419369")) %>% 
  mutate(Ligand = str_replace(string = Ligand, pattern = "sch", replacement = "SCH23390")) %>% 
  spread(region, icc) %>% 
  mutate("Cor (Region 1)" = Corr_mrtm2$Cor_1) %>% 
  mutate("Cor (Region 2)" = Corr_mrtm2$Cor_2) %>% 
  mutate("Bias (Region 1)" = Bias_mrtm2$Bias_1) %>% 
  mutate("Bias (Region 2)" = Bias_mrtm2$Bias_2) %>% 
  arrange(Ligand, desc(Model))

# kable(mrtm2_table) %>% 
#   save_kable(., "mrtm2.png")
```

``` r
k <- mrtm2 %>% 
  spread(software, Value) %>% 
 ggplot(aes(x = kinfitR, y = pmod, color = subjname)) + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() +
  geom_line(aes(group=subjname), linetype="dotted") + 
  facet_grid(cols = vars(region), vars(Ligand), scales = "fixed") + theme(legend.position = "none") +
  ylab("PMOD") + xlab("kinfitr")
```

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.6.2  magrittr_1.5    tools_3.6.2     htmltools_0.4.0
    ##  [5] yaml_2.2.1      Rcpp_1.0.3      stringi_1.4.6   rmarkdown_2.1  
    ##  [9] knitr_1.28      stringr_1.4.0   xfun_0.12       digest_0.6.25  
    ## [13] rlang_0.4.5     evaluate_0.14
