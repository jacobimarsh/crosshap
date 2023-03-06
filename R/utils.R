#' Mode utility function
#'
#' @param x Input vector
#'
#' @export
#'

arith_mode <- function(x) {
  ux <- base::unique(x)
  ux[base::which.max(base::tabulate(base::match(x, ux)))]
}

#' Mean utility function
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
#' If your correlation matrix does not have rownames and column names, a VCF
#' will need to be provided so it can be added with read_LD().
#'
#' @param LDin Square correlation matrix
#' @param vcf VCF object created by read_vcf() that can be used to assign column names
#'
#' @export
#'
#' @return A tibble.
#'

read_LD <- function(LDin, vcf = NULL){
  if(is.null(vcf)){
    LD <- data.table::fread(LDin, nThread = 10) %>%  tibble::as_tibble() %>%  tibble::column_to_rownames("V1")
  } else{
    LD <- data.table::fread(LDin, nThread = 10, header = F) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(tempID = vcf$ID) %>%
      tibble::column_to_rownames("tempID")
    colnames(LD) <- vcf$ID
  }
  return(LD)
  }

#' Read VCF to tibble
#'
#' Dashes,'-', in individual names are recoded to '.' for downstream
#' compatability.
#'
#' @param VCFin Input VCF
#' @export
#'
#' @return A tibble.
#'

read_vcf <- function(VCFin){
  vcf <- data.table::fread(VCFin, nThread = 10) %>%  tibble::as_tibble() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),~ base::gsub(":.*","",base::gsub("/","|",.)))) %>%
    dplyr::mutate(POS = as.numeric(POS))

  colnames(vcf) <- gsub('-','.',colnames(vcf))
  return(vcf)
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

read_pheno <- function(Phenoin){
  data.table::fread(Phenoin, header = F) %>% tibble::as_tibble() %>%
    dplyr::filter(V1 !='Ind' & V2 != 'Pheno') %>%
    dplyr::rename('Ind' = V1, 'Pheno' = V2) %>%
    dplyr::mutate(Ind = gsub('-','.',Ind), Pheno = as.double(Pheno))
}

#' Read metadata to tibble
#'
#' Requires two column text file without a header (Ind | Metadata)
#'
#' @param Metain Input phenotype file
#'
#' @export
#'
#' @return A tibble.
#'

read_metadata <- function(Metain){
  data.table::fread(Metain, header = F) %>% tibble::as_tibble() %>%
    dplyr::filter(V1 !='Ind' & V2 != 'Metadata') %>%
    dplyr::rename('Ind' = V1, 'Metadata' = V2) %>%
    dplyr::mutate(Ind = gsub('-','.',Ind), Metadata = gsub('-','.',Metadata))
}



