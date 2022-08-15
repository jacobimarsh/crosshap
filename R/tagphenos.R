#' Calculate phenotypic associations for SNPs in each cluster
#'
#' @param MGfile SNP marker groups clustered using DBscan.
#' @param vcf Input VCF for region of interest.
#' @param pheno Input numeric phenotype data for colouring clustree.
#'
#' @return
#' @export
#'
#' @examples
#'
run_clustree <- function(MGfile, vcf, pheno) {

vcf_long <- vcf %>%
  rownames_to_column("POS") %>%
  left_join(MGfile) %>%
  gather(ID, key, 2:(ncol(.)-1))

VarFile <- vcf_long %>%
  left_join(phen_early) %>%
  group_by(POS, cluster, key) %>%
  summarize(nInd = n(),avPheno=mean(Pheno, na.rm = T))

return(VarFile)
}
