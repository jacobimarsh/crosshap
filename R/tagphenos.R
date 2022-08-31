#' Calculate SNP phenotypic associations
#'
#' Internal function that splits SNPs by their different allelic states and
#' calculates frequencies across the population. In addition, the mean phenotype
#' is calculated across all individuals that share a given allele for each loci.
#'
#' @param MGfile SNP marker groups clustered using DBscan.
#' @param vcf Input VCF for region of interest.
#' @param pheno Input numeric phenotype data for each individual.
#'
#' @return
#' @export
#'
#' @examples
#'
tagphenos <- function(MGfile, bin_vcf, pheno) {

#Split by allele type
bin_vcf_long <- bin_vcf %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(MGfile, by = "ID") %>%
  tidyr::gather(Ind, key, 2:(base::ncol(.)-3))

#Calculate phenotypic association of each allele type for each SNP
VarFile <- bin_vcf_long %>%
  dplyr::left_join(pheno, by = "Ind") %>%
  dplyr::group_by(ID, MGs, key) %>%
  dplyr::summarize(nInd = dplyr::n(),avPheno=base::mean(Pheno, na.rm = T), .groups = 'keep')

return(VarFile)
}
