---
  title: "CrossHap_clustering"
output: html_document
date: '2022-07-19'
editor_options:
  chunk_output_type: console
---

  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(data.table)
library(tidyverse)
library(factoextra)
library(janitor)
library(dbscan)
library(fpc)
library(tsne)
library(patchwork)
setwd("/mnt/ws/jmarsh/Rscrips/jacob_haplopipe")
```


```{r pre}
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#Read correlation matrix between all SNPs in region
LD <- fread('data/LD_100kb.mtx') %>%
  as_tibble() %>%
  column_to_rownames(var = "V1")
##Parse a VCF as a matrix and convert alleles to base 3 integers
vcf <- fread('data/impu_100kb_pdh1.vcf') %>%
  as_tibble() %>%
  select(-c(1,2,4,5,6,7,8,9)) %>%
  column_to_rownames('ID') %>%
  mutate_all(function(x){ifelse(x=='0|0',0,
                                ifelse(x=='1|0'|x=='0|1',1,
                                       ifelse(x=='1|1',2,'failsave')))})

phen_early <- fread('data/early_shatter565.txt') %>% as_tibble() %>% rename('ID' = V1, 'Pheno' = V2)

elon <- c(0.5,1,1.5,2)
```


```{r cluster loop, echo=FALSE}
for (rez in elon)
{
  db40 <- dbscan(LD, eps = rez, minPts = 40)

  db40_clust <- cbind(LD[, 0], db40$cluster) %>%
    rownames_to_column("POS") %>%
    rename("cluster" = "db40$cluster") %>% as_tibble()
  ##See if can export cluster names from dbscan

  assign(paste("clustree_export_clusts", rez,
               sep = ""), db40_clust)

  write_tsv(db40_clust,paste0("SNP_MGs_MPx_E", rez, '.tsv'),
            quote = 'none')

  ##Step 1: Calculate  alternate vs reference counts for all SNPs in a given cluster

  pseduoNP <- create_pseudo_SNP(db40 = db40, vcf = vcf)

  ##Step 2: Start creating pseudoSNP with all MGs
  for (vel in c(2:max(db40$cluster)))
  {
    db40_cvel <- db40_clust %>%
      filter(cluster == vel) %>%
      as_tibble()

    cvel_vcf <- vcf %>%
      rownames_to_column() %>%
      filter(rowname %in% db40_cvel$POS) %>%
      column_to_rownames()

    i_modes_vel <- apply(cvel_vcf %>%
                           sapply(as.double), 2, mode) %>%
      as_tibble()

    pseudoSNP <- bind_cols(pseudoSNP,i_modes_vel)

    colnames(pseudoSNP)[vel + 1] <- paste0("MG", vel)
  }
  ###################UP TO HERE @JAKOB###################
  #het to full hap
  het_pseudoSNP <- pseudoSNP %>% column_to_rownames("ID") %>%
    mutate_all(function(x) {
      gsub(1, 2, x,
           fixed = T)
    }) %>% rownames_to_column("ID") %>% as_tibble()

  write_tsv(het_pseudoSNP,paste0("ID_MGs_MPx_E", rez, '.tsv'),
            quote = 'none')

  ####NEED TO MAKE IT MUTATE ALL >BUT> ID

  #hapCounts <- pseudoSNP %>%
  # count(MG1, MG2, MG3, MG4, MG5, MG6, MG7) %>%
  #  count() %>%
  #  as_tibble() %>%
  #  arrange(by_group = -n)

  cnames <- colnames(select(pseudoSNP, -ID))

  het_hapCounts <- het_pseudoSNP %>%
    gather(mgs,value,2:ncol(.)) %>%
    group_by(ID) %>%
    mutate(mgs_new=paste(value, collapse = '_')) %>%
    group_by(mgs_new) %>%
    tally() %>% mutate(n=n/(ncol(pseudoSNP)-1)) %>%
    separate(col = mgs_new, into = cnames, sep = "_") %>%
    as_tibble() %>%
    arrange(by_group = -n)

  over20_hhCounts <- filter(het_hapCounts, n > 9) %>%
    mutate(hap_eps=LETTERS[1:nrow(.)]) %>%
    rename(!!paste0("hap_eps",rez) := 'hap_eps')#%>%

  clustered_hpS <- inner_join(het_pseudoSNP, over20_hhCounts[, 1:ncol(over20_hhCounts)],
                              by = (colnames(over20_hhCounts[, 2:ncol(over20_hhCounts)]) =
                                      c(colnames(het_pseudoSNP[, 2:ncol(het_pseudoSNP)]))))

  clustree_expor_hpS <- mutate(clustered_hpS[, c(1, ncol(clustered_hpS))])
  ##multiplex with a few different epsilon parameters
  rexport_hpS <- paste("Haplotype_assignments_MP", "_E", rez, ".txt",
                       sep = "")

  assign(rexport_hpS, clustree_expor_hpS)

  write_tsv(clustree_expor_hpS,paste0("ID_Haps_MPx_E", rez, '.tsv'),
            quote = 'none')
}

pre_clustree <- get(paste("Haplotype_assignments_MP","_E",elon[1],".txt",sep=""))
for (e in elon[2:length(elon)])
{
  pre_clustree <- pre_clustree %>%
    left_join(get(paste("Haplotype_assignments_MP","_E",e,".txt",sep="")), by="ID")
}
pre_clustree[is.na(pre_clustree)] <- "0"
```

```{r}
##CURRENT WORKING
tet <- mutate(over20_hhCounts, LETTERS=LETTERS[1:nrow(over20_hhCounts)])
```

```{r glue into clustree input}
library(clustree)

pre_clustree_phen <- left_join(pre_clustree, phen_early, by = 'ID')

mean_na.rm <- function(x){mean(x,na.rm=T)}



ctree <- clustree(pre_clustree_phen,
                  prefix = 'hap_eps',
                  node_colour = 'Pheno',
                  node_colour_aggr = 'mean_na.rm',
                  edge_width = 1,
                  node_alpha = 1) +
  scale_colour_gradient(limits=c(max(top_frac(pre_clustree_phen,-0.1,Pheno)$Pheno),
                                 min(top_frac(pre_clustree_phen,0.1,Pheno)$Pheno)),
                        high = "#8ADD81",
                        low = "#6870F6",
                        oob = scales::squish,
                        name = 'Pheno') +
  scale_edge_color_continuous(high = 'black',
                              low = 'grey80') +
  labs(size = 'nIndividuals', edge_alpha = "Proportion") +
  guides(edge_color = "none", size = guide_legend(order = 1))

lbls <- tibble(yval=0:(length(elon)-1),
               xval = max(ctree[["data"]][["x"]])*1.1,
               labelval = paste0("\u03b5"," = ",elon))

labeled_ctree <- ctree + geom_text(data = lbls, aes(x=xval, y=yval, label=labelval), hjust = 0)


##0:(length(elon)-1)

ggsave("clustree_MPx.pdf",
       ctree,
       device = 'pdf',
       dpi = 300,
       height = 9,
       width = 16,
       units = 'in')

```


