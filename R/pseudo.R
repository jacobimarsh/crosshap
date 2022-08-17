#' Create pseudoSNPs from SNP marker groups
#'
#' Internal function that calls the most common allelic states for each SNP
#' marker group across individuals, returning dummy SNPs for each marker group
#' mimicking the vcf format. Options to be added for treatment of heterozygotes.
#'
#' @param MGfile SNP marker groups clustered using DBscan.
#' @param vcf Input VCF for region of interest.
#'
#' @return
#' @export
#'
#' @example
#' create_pseudoSNP(MGfile, pdh1_100kb_vcf)
#'
create_pseudoSNP <- function(MGfile, vcf) {

#Extract SNPs in first MG cluster (MG1)
  db40_c1 <- MGfile %>%
    dplyr::filter(cluster == 1) %>% tibble::as_tibble()
  c1_vcf <- vcf %>%
    tibble::rownames_to_column() %>% dplyr::filter(rowname %in% db40_c1$POS) %>% tibble::column_to_rownames()

#Calculate most common (alternate or ref) allele for MG1 across all individuals
  i_modes_1 <- base::apply(c1_vcf %>% base::sapply(as.double), 2, mode) %>% tibble::as_tibble() %>%
    pull(value)

#Build a dummy VCF with one pseudoSNP position (MG1)
  pseudoSNP <- tibble::tibble(ID = base::colnames(c1_vcf), MG1 = i_modes_1)

#Repeat for all other MGs, iteratively adding to pseudoSNP vcf
for (vel in c(2:base::max(MGfile$cluster))) {

    db40_cvel <- MGfile %>%
      dplyr::filter(cluster == vel) %>% tibble::as_tibble()
    cvel_vcf <- vcf %>%
      tibble::rownames_to_column() %>% dplyr::filter(rowname %in% db40_cvel$POS) %>% tibble::column_to_rownames()
    i_modes_vel <-
      base::apply(cvel_vcf %>% base::sapply(as.double), 2, mode) %>% tibble::as_tibble()
    pseudoSNP <- dplyr::bind_cols(pseudoSNP, i_modes_vel)
    base::colnames(pseudoSNP)[vel + 1] <- base::paste0("MG", vel)
  }

#Convert heterozygous marker groups to homozygous alternate (optional)
  het_pseudoSNP <- pseudoSNP %>%
    dplyr::mutate_if(is.numeric,function(x) {base::gsub(1, 2, x,fixed = T)})

  return(het_pseudoSNP)
}
