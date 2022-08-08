#' Utility functions and example data
#'
#' @param x Placeholder for mode and mean functions.
#'
#' @noRd
#'
##Dependencies
library(data.table)
library(tidyverse)
library(factoextra)
library(janitor)
library(dbscan)
library(tsne)
library(patchwork)
library(clustree)
library(conflicted)

##Arithmetic mode for calling each individual's allelic states for pseudoSNP marker groups
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

##Arithmetic mean of phenotype scores for each haplotype, removing NAs
mean_na.rm <- function(x){
  mean(x,na.rm=T)
}

setwd("/Users/jmarsh96/Desktop/bash_misc/crosshap_data")

#Read correlation matrix between all SNPs in region
LD <- fread('data/LD_100kb.mtx', nThread = 10) %>%  as_tibble() %>%  column_to_rownames("V1")

##Parse a VCF as a matrix and convert alleles to base 3 integers
vcf <- fread('data/impu_100kb_pdh1.vcf', nThread = 10) %>%  as_tibble() %>%
  select(-c(1,2,4:9)) %>% column_to_rownames('ID') %>%
  mutate_all(function(x){ifelse(x=='0|0',0,ifelse(x=='1|0'|x=='0|1',1,ifelse(x=='1|1',2,'failsave')))})

##Read phenotype data from two column text file without a header (ID | Pheno)
phen_early <- fread('data/early_shatter565.txt') %>% as_tibble() %>%
  rename('ID' = V1, 'Pheno' = V2)

##DBscan epsilon values to be tested and compared using clustree visualization
elon <- seq(1,2,by=0.2)
