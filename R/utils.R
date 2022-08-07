#' Title
#'
#' @param x
#'
#' @noRd
#'

library(data.table)
library(tidyverse)
library(factoextra)
library(janitor)
library(dbscan)
library(tsne)
library(patchwork)
library(clustree)
library(conflicted)

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mean_na.rm <- function(x){
  mean(x,na.rm=T)
}

#Read correlation matrix between all SNPs in region
LD <- fread('data/LD_100kb.mtx', nThread = 10) %>%  as_tibble() %>%  column_to_rownames("V1")

##Parse a VCF as a matrix and convert alleles to base 3 integers
vcf <- fread('data/impu_100kb_pdh1.vcf', nThread = 10) %>%  as_tibble() %>%
  select(-c(1,2,4:9)) %>% column_to_rownames('ID') %>%
  mutate_all(function(x){ifelse(x=='0|0',0,ifelse(x=='1|0'|x=='0|1',1,ifelse(x=='1|1',2,'failsave')))})

phen_early <- fread('data/early_shatter565.txt') %>% as_tibble() %>% rename('ID' = V1, 'Pheno' = V2)

elon <- seq(1,2,by=0.2)
