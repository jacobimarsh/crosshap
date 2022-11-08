#' Utility functions and example data
#'
#' @param x Input vector
#'
#' @export
#'

mode <- function(x) {
  ux <- base::unique(x)
  ux[base::which.max(base::tabulate(base::match(x, ux)))]
}

#' Utility functions and example data for mean
#'
#' @param x Input vector
#'
#' @export
#'

mean_na.rm <- function(x){
  base::mean(x,na.rm=T)
}



#' Read LD correlation matrix to tibble
#'
#' @param LDin Correlation matrix output from PLINK
#'
#' @export
#'
#' @return A tibble.
#'
#' @examples $ plink --r2 square --vcf your_region.vcf
#' read_LD(plink.out)
#'

read_LD <- function(LDin){
  data.table::fread(LDin, nThread = 10) %>%  tibble::as_tibble() %>%  tibble::column_to_rownames("V1")
}

#' Read VCF to tibble
#'
#' @param VCFin Input VCF
#' @export
#'
#' @return A tibble.
#'
#' @examples $ bgzip your.vcf
#' $ tabix -p vcf your.vcf
#' $ tabix your.vcf.gz chr1:15,000,000-15,250,000 > your_region.vcf
#' read_vcf(your_region.vcf)
#'

read_vcf <- function(VCFin){
  data.table::fread(VCFin, nThread = 10) %>%  tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),~ base::gsub(":.*","",base::gsub("/","|",.)))) %>%
    dplyr::mutate(POS = as.numeric(POS))
}

#' Read phenotype data to tibble
#'
#' Requires two column text file without a header (Ind | Pheno)
#'
#' @param Phenoin Input phenotype file
#' @export
#'
#' @return A tibble.
#'
#' @example read_pheno(yield.txt)
#'

read_pheno <- function(Phenoin){
  data.table::fread(Phenoin) %>% tibble::as_tibble() %>%
    dplyr::rename('Ind' = V1, 'Pheno' = V2)
}

#' Read metadata to tibble
#'
#' Requires two column text file without a header (Ind | Metadata)
#'
#' @param Metain Input phenotype file
#' @export
#'
#' @return A tibble.
#'
#' @example read_metadata(country_of_origin.txt)
#'

read_metadata <- function(Metain){
  data.table::fread(Metain, header = F) %>% tibble::as_tibble() %>%
    dplyr::rename('Ind' = V1, 'Metadata' = V2)
}



