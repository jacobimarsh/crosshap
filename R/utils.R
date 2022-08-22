#' Utility functions and example data
#'
#' @param x Placeholder for mode and mean functions.
#'
#' @noRd
#'
##Dependencies to install
#library(data.table)
#library(tidyverse)
#library(factoextra)
#library(janitor)
#library(dbscan)
#library(tsne)
#library(patchwork)
#library(clustree)

##Arithmetic mode for calling each individual's allelic states for pseudoSNP marker groups
mode <- function(x) {
  ux <- base::unique(x)
  ux[base::which.max(base::tabulate(base::match(x, ux)))]
}

##Arithmetic mean of phenotype scores for each haplotype, removing NAs
mean_na.rm <- function(x){
  base::mean(x,na.rm=T)
}

#Read correlation matrix between all SNPs in region
read_LD <- function(LDin){
  data.table::fread(LDin, nThread = 10) %>%  tibble::as_tibble() %>%  tibble::column_to_rownames("V1")
}

##Parse an imputed VCF as a matrix and convert alleles to base 3 integers
read_vcf <- function(VCFin){
  data.table::fread(VCFin, nThread = 10) %>%  tibble::as_tibble() #%>%
#    dplyr::select(-c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
#    dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,base::ifelse(x=='1|0'|x=='0|1',1,base::ifelse(x=='1|1',2,'failsave')))})
}

position_vcf <- function(VCFin){
  data.table::fread(VCFin, nThread = 10) %>%  tibble::as_tibble() %>%
    dplyr::select(c(2,3))
}

##Read phenotype data from two column text file without a header (ID | Pheno)
read_pheno <- function(phenoin){
  data.table::fread(phenoin) %>% tibble::as_tibble() %>%
    dplyr::rename('Ind' = V1, 'Pheno' = V2)
}


