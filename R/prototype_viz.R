---
title: "Stripped_SUS3"
author: "Jacob Marsh"
date: "01/03/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(root.dir = "C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3", 
                      echo = TRUE)
library(tidyverse)
library(dplyr)
library(viridis)
library(patchwork)
library(ggplot2)
library(scales)
library(ggsci)
library(ggalt)

```


```{r Import}
HF2 <- read.csv("C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3/U_S_haps_fin3.txt") %>% mutate(altalleles = ifelse(altalleles=='',NA,altalleles))

AlleleFile <-read.csv("C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3/U_S_allele_fin3.txt")

PhenoSum <- read.csv("C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3/allpheno_resum3.txt") 

AlleleCounts <- read.csv("C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3/ACAN_tagSNPs.txt")

happhengrp <- read.csv("C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3/happhengrp.txt")

TagPercDiffs <- read.csv("C:/Users/21485753/Desktop/cqProt003/UpSet_Plot/Rinputs/Rinputs2/Rinputs3/percdifftags3.csv") %>%
  subset(select = -X) %>% filter(TAGGING >= 31604127 & TAGGING <= 31777346)

```

```{r Reformat}
#A top plot
grplongHap <- tidyr::gather(HF2, 
                            "dom_status", 
                            "n", 
                            3:6)

#B bot plot
happhengrp

#C left plot
AlPhen <- merge(x = AlleleFile, 
                y = PhenoSum, 
                by.x = "allele", 
                by.y = "Site")

#D right plot
TagProp2 <- AlleleCounts %>% 
  mutate(AltAF = AC / AF) %>% 
  subset(select = -c(AC,AF))
TagAlleles2 <- right_join(TagPercDiffs, 
                          AlleleFile, 
                          by =c("SNP" = "allele"))
TagAlProp2 <- right_join(TagAlleles2, 
                         TagProp2, 
                         by = c("TAGGING" = "SITE")) %>% 
  mutate(percdiff = 100*percdiff) %>% filter(!is.na(SNP))

#E dot plot
HapAlMatrix <- HF2 %>% separate_rows(altalleles,sep=';') %>% 
     mutate(value=1) %>% spread(altalleles,value,fill = 0) %>% select(-`<NA>`)
intersect <- HapAlMatrix %>% 
  as_tibble() %>% 
  gather(allele,
         present,
         7:ncol(.)) %>% 
  mutate(present=as.factor(present)) %>%
  mutate(allele = str_remove(allele, "altalleles_"))
intersect_lines <- intersect %>% 
  filter(present == 1) %>% 
  group_by(hap) %>% 
  mutate(allele = as.numeric(allele)) %>% 
  summarise(max = max(allele), 
            min = min(allele)) %>% 
  mutate(min = as.character(min),
         max = as.character(max))

#colours
npg_col = pal_npg("nrc")(9)
col_list <- c('wt'=npg_col[8],
              'lr' = npg_col[3],
              'oc' =npg_col[2],
              'mc' =npg_col[4])
```


```{r Dotplot (E)}
E <- ggplot() +
  geom_segment(data = intersect_lines, col = 'grey', 
               aes(x = hap, 
                   xend = hap,
                   y = min, 
                   yend = max),
               size = 1.5)+
  geom_point(data = intersect,
             aes(hap,
                 as.character(allele),
                 fill =as.factor(present),
                 size = 2),
             col ='black', 
             pch = 21) +
  scale_fill_manual(values = c('white','black', 'white'))+
  theme_minimal()+
  theme(legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), 
                           "cm"),
        plot.title = element_blank(),
        axis.text.y = element_text(size=10, face='bold', color = c("black", "black", "black", "black", "black", "red")),
        axis.text.x = element_text(size=10, face = 'bold', color = "black")) +
  ylab("Marker group") +
  xlab("Haplotype combination") +
  scale_y_discrete(position = "left", labels = c(paste0("M0",as.character(6:1))))
E
```

```{r Top plot (A)}
A <- ggplot(data = grplongHap %>% 
  mutate(dom_status=factor(dom_status,
                           levels = c('wt','lr','oc','mc'))), 
       aes(fill=dom_status, 
           y=n, 
           x=hap)) + 
  geom_bar(position="stack", 
           stat="identity") +
  theme_minimal() +
  scale_fill_manual("Domestication status",
                    values = col_list,
                    labels = c("Wild-Type", "Landrace", "Old Cultivar", "Modern Cultivar")) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(6, 
                               "mm"),
        axis.text.x = element_text(face = "bold", 
                                   size = 10, color = c("black", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red")),
        plot.margin = unit(c(0,0,0,0), 
                           "cm"),
        legend.position = c(0.8,0.8),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Population size") +
  xlab("Haplotype combination")

A
```

```{r Bot plot (B)}
B <- ggplot(data = happhengrp %>% 
  mutate(Grp=factor(Grp,
                           levels = c('wt','lr','ocult','mcult')))
               ) + 
  geom_jitter(aes(x = Hap, 
                  y = Prot,
                  fill= Grp),
              alpha=0.25, 
              pch=21, 
              width=0.2) +
  scale_fill_manual("Domestication status",values = col_list, guide = "none") +
    geom_crossbar(data= aggregate(happhengrp[, 3], 
                                list(happhengrp$Hap), 
                                median),
                aes(x= as.factor(Group.1),
                    y=x,
                    xmin= as.factor(Group.1) -1,
                    xmax=as.factor(Group.1) +1,
                    ymin=x,
                    ymax=x,
                    colour=x)) +
  scale_colour_gradient('Mean',
                        low='red', 
                        high='green',
                        limits=c(max(top_frac(happhengrp,
                                              -0.2,
                                              Prot)$Prot),
                                 min(top_frac(happhengrp,
                                              0.2,
                                              Prot)$Prot)), 
                        oob = squish) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.margin = unit(c(0,0,0,0), 
                           "cm"),
        axis.text.y = element_text(face = "bold",
                                   size = 10)) +
  ylab("% Seed protein") +
  scale_x_discrete(position = "bottom") +
  scale_y_continuous(position = "left", breaks = scales::pretty_breaks(n = 5))

B
```

```{r Patching AB}

AB <- A + B + plot_layout(design = "AA
                          BB")

AB
```

```{r Left plot (C)}
C <- AlPhen %>% 
  mutate(Type=factor(Type,levels = c("REF","MISS","HET","HETMISS","ALT"))) %>%
  ggplot(aes(x = nInd, 
             y = as.character(allele),
             fill = Type,
             color = Type),
         ylim()) +
  geom_bar(aes(), 
           position = "stack", 
           stat = "identity",
           colour = "black",
           width = 0.8) +
  scale_x_reverse(breaks = scales::pretty_breaks(n = 5), 
                     expand = c(0,0)) +

  theme_void() + 
  theme(axis.text.y = element_blank(), 
        axis.title.x = element_text(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.key.size = unit(5, 
                               "mm"),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0.1), 
                           "cm"),
        plot.title = element_blank()) +
# scale_color_manual(values = c('black', 'black', 'black', 'black', 'black', 'red'), 
#                    guide=F) +
 scale_fill_manual(values = rev(c('#440154FF', "#238A8DFF", "#73D055FF", "#FDE725FF", "#FFFFFF"))) +
  xlab("Allele count") +
  scale_y_discrete(position = "right", labels = c(paste0("M",as.character(20:10)), paste0("M0",as.character(7:1))))
  
C
```

```{r Right plot (D)}
D <- ggplot() +
  geom_jitter(data = TagAlProp2,
              aes(x = abs(percdiff),
                  y = as.character(SNP),
                  fill = AltAF),
              alpha = 0.25,
              pch = 21,
              height = 0.25) +
  scale_fill_gradient('Minor allele frequency',
                      low = 'white',
                      high = '#440154FF') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 10, color = c("black", "black", "black", "black", "black", "red")),
        plot.margin = unit(c(0,0.1,0,0), 
                           "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.position = "bottom",
        legend.key.size = unit(5, 
                               "mm"),
        plot.title = element_blank(),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        axis.title.x = element_text()) +
  xlab("Seed protein association [%]") +
  scale_y_discrete(position = "left", labels = c(paste0(" M0",as.character(6:1))))

D
```

```{r Patching}

AB <- A + B + plot_layout(design = "AA
                          BB")
AB
ggsave('AB3.pdf',
       AB, 
       device = 'pdf', 
       dpi = 1800,
       height = 195,
       width = 110,
       units = "mm")

CD <- C + D + plot_layout(design = "CD
                          CD")
CD

ggsave('CD3.pdf',
       CD, 
       device = 'pdf', 
       dpi = 1800,
       height = 97.875,
       width = 174,
       units = "mm")

E

ggsave('E3.pdf',
       E,
       device = 'pdf',
       dpi = 1800,
       height = 120,
       width = 120,
       units = "mm")

```



```{r BigUpSet}
ABCDE <- A + B + C + D + E + guide_area() + plot_layout(design = "FA#
                                         CED
                                         #B#", guides = "collect")

ABCDE


E2 <- E + theme(#axis.title.y = element_blank(),
          #axis.title.x = element_blank()) + 
  axis.text.x = element_text(face = "bold", 
                                   size = 10, color = c("black", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red")))

ABCDE <- A + B + C + D + E2 + guide_area() + plot_layout(design = "FA#
                                         CED
                                         #B#", guides = "collect")

ggsave('ABCDE4.pdf',
       ABCDE,
       device = 'pdf',
       dpi = 1800,
       height = 240,
       width = 240,
       units = "mm")


ts5 <- fread("C:/Users/21485753/Desktop/cqProt003/mg_locations/tableS5.txt", sep = "\t") %>% as_tibble() %>% 
  mutate(MG=as.numeric(gsub('M0','',MG)))

G <- ggplot() + 
  geom_segment(data=filter(ts5,Type=='nrep'),aes(x=Pos,xend=Pos,y=MG-0.2,yend=MG+0.2),size=0.2)+
  geom_segment(data=filter(ts5,Type=='rep'),aes(x=Pos,xend=Pos,y=MG-0.2,yend=MG+0.2),col='red')+
  geom_point(data=filter(ts5,Type=='rep'),aes(Pos,MG), pch=23,fill='red',size=2)+
  scale_x_continuous(breaks= pretty_breaks(n=3), labels = comma)+
  scale_y_reverse(breaks=1:6, labels=paste('M0',1:6,sep=''), position = "left")+
  labs(x='Position', y='Marker group')+
  theme_minimal()+
  theme(#axis.text.y = element_text(face='bold',color= c("red", "black", "black", "black", "black", "black")),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), 
                           "cm"))+
  expand_limits(x=c(31450000,31900000))+
  NULL

G




ABGDE <- A + B + G + D + E2 + guide_area() + plot_layout(design = "FA#
                                         CED
                                         #B#", guides = "collect")

ABGDE

ggsave('ABGDE7.pdf',
       ABGDE,
       device = 'pdf',
       dpi = 1800,
       height = 240,
       width = 240,
       units = "mm")




FBCDE <- H + B + C + D + E2 + guide_area() + plot_layout(design = "FA#
                                         CED
                                         #B#", guides = "collect")

```
