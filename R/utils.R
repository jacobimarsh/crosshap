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

#Read correlation matrix between all SNPs in region
read_LD <- function(LDin){
  fread(LDin, nThread = 10) %>%  as_tibble() %>%  column_to_rownames("V1")
  }

##Parse an imputed VCF as a matrix and convert alleles to base 3 integers
read_vcf <- function(VCFin){
  fread(VCFin, nThread = 10) %>%  as_tibble() %>%
    select(-c(1,2,4:9)) %>% column_to_rownames('ID') %>%
    mutate_all(function(x){ifelse(x=='0|0',0,ifelse(x=='1|0'|x=='0|1',1,ifelse(x=='1|1',2,'failsave')))})
}

##Read phenotype data from two column text file without a header (ID | Pheno)
read_pheno <- function(phenoin){
  fread(phenoin) %>% as_tibble() %>%
    rename('ID' = V1, 'Pheno' = V2)
}


