#' Utility functions and example data
#'
#' @param x Placeholder for mode and mean functions.
#'
#' @export
#'

mode <- function(x) {
  ux <- base::unique(x)
  ux[base::which.max(base::tabulate(base::match(x, ux)))]
}

#' Utility functions and example data for mean
#'
#' @param x Placeholder for mode and mean functions.
#'
#' @export
#'
mean_na.rm <- function(x){
  base::mean(x,na.rm=T)
}



#' Read correlation matrix between all SNPs in region
#'
#' @param LDin Correlation matrix output from PLINK
#' @return A read in table of blah.
#'
#' @export
read_LD <- function(LDin){
  data.table::fread(LDin, nThread = 10) %>%  tibble::as_tibble() %>%  tibble::column_to_rownames("V1")
}

#' Read VCF
#'
#' @param VCFin Input VCF
#' @export
read_vcf <- function(VCFin){
  data.table::fread(VCFin, nThread = 10) %>%  tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),~ base::gsub(":.*","",base::gsub("/","|",.)))) #%>%
#    dplyr::select(-c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
#    dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,base::ifelse(x=='1|0'|x=='0|1',1,base::ifelse(x=='1|1',2,'failsave')))})
}

#' @param VCFin Input VCF
#' @describeIn readVCF read_VCF_pos
#' @export
position_vcf <- function(VCFin){
  data.table::fread(VCFin, nThread = 10) %>%  tibble::as_tibble() %>%
    dplyr::select(c(2,3))
}

#' Read phenotype data from two column text file without a header (ID | Pheno)
#'
#' @param Phenoin Input phenotype file
#' @export
read_pheno <- function(Phenoin){
  data.table::fread(Phenoin) %>% tibble::as_tibble() %>%
    dplyr::rename('Ind' = V1, 'Pheno' = V2)
}

