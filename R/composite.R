---
  title: "CrossHap_clustering"
output: html_document
date: '2022-07-19'
editor_options:
  chunk_output_type: console
---

knitr::opts_chunk$set(echo = TRUE)
library(data.table)
library(tidyverse)
library(factoextra)
library(janitor)
library(dbscan)
#library(fpc)
library(tsne)
library(patchwork)
library(clustree)
library(conflicted)
conflict_prefer_all('dplyr')
setwd("/Users/jmarsh96/Desktop/bash_misc/crosshap_data")


mean_na.rm <- function(x){mean(x,na.rm=T)}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#Read correlation matrix between all SNPs in region
LD <- fread('data/LD_100kb.mtx', nThread = 10) %>%
  as_tibble() %>%
  column_to_rownames("V1")
##Parse a VCF as a matrix and convert alleles to base 3 integers
vcf <- fread('data/impu_100kb_pdh1.vcf', nThread = 10) %>%
  as_tibble() %>%
  select(-c(1,2,4:9)) %>%
  column_to_rownames('ID') %>%
  mutate_all(function(x){ifelse(x=='0|0',0,ifelse(x=='1|0'|x=='0|1',1,ifelse(x=='1|1',2,'failsave')))})

phen_early <- fread('data/early_shatter565.txt') %>% as_tibble() %>% rename('ID' = V1, 'Pheno' = V2)

elon <- seq(1,2,by=0.2)

##Clustering loop

for (arez in elon){
  db40 <- dbscan::dbscan(LD, eps = arez, minPts = 40)
  db40_clust <- tibble(POS=rownames(LD),cluster=db40$cluster)
  assign(paste("export_clusts", arez,sep = ""), db40_clust)
  het_pseudoSNP <- create_pseudoSNP(db40_clust,vcf)
  assign(paste("het_pseudoSNP_MPx", "_E", arez,sep = ""), het_pseudoSNP)
  cnames <- colnames(select(het_pseudoSNP, -ID))
  het_hapCounts <- het_pseudoSNP %>%
    gather(mgs,value,2:ncol(.)) %>%
    group_by(ID) %>%
    mutate(mgs_new=paste(value, collapse = '_')) %>%
    group_by(mgs_new) %>%
    tally() %>% mutate(n=n/(ncol(het_pseudoSNP)-1)) %>%
    separate(col = mgs_new, into = cnames, sep = "_") %>%
    as_tibble() %>%
    arrange(by_group = -n)
  over20_hhCounts <- filter(het_hapCounts, n > 9) %>%
    mutate(hap_eps=LETTERS[1:nrow(.)]) %>%
    rename(!!paste0("hap_eps",arez) := 'hap_eps')
  clustered_hpS <- left_join(het_pseudoSNP, over20_hhCounts) %>%
    mutate_if(is.character,function(x){replace_na(x, '0')}) %>%
    select(1,ncol(.))
  assign(paste("Haplotype_assignments_MP", "_E", arez,sep = ""), clustered_hpS)
}

pre_clustree <- get(paste("Haplotype_assignments_MP","_E",elon[1],sep=""))
for (drez in elon[2:length(elon)]){
  pre_clustree <- pre_clustree %>%
    left_join(get(paste("Haplotype_assignments_MP","_E",drez,sep="")))
}

pre_clustree_phen <- left_join(pre_clustree, phen_early)
#plot
ctree <-
  clustree(pre_clustree_phen,prefix = 'hap_eps',node_colour = 'Pheno',node_colour_aggr = 'mean_na.rm',edge_width = 1,node_alpha = 1)+
  scale_colour_gradient(limits=c(max(top_frac(pre_clustree_phen,-0.1,Pheno)$Pheno),min(top_frac(pre_clustree_phen,0.1,Pheno)$Pheno)),high = "#8ADD81",low = "#6870F6",oob = scales::squish,name = 'Pheno') +
  scale_edge_color_continuous(high = 'black',low = 'grey80') +
  labs(size = 'nIndividuals', edge_alpha = "Proportion") +
  guides(edge_color = "none", size = guide_legend(order = 1))
#create label data
#x and y extractd from clstree object
lbls <- tibble(xval = max(ctree[["data"]][["x"]])*1.1,
               yval=0:(length(elon)-1),
               labelval = paste0("\u03b5"," = ",elon))
#re-plot using label data
labeled_ctree <- ctree +
  geom_text(data = lbls, aes(x=xval, y=yval, label=labelval), hjust = 0)

#save plot
ggsave("clustree_MPx.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')


