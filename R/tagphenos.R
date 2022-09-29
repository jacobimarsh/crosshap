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
#' @export
#'
#'
tagphenos <- function(MGfile, bin_vcf, pheno) {

#Split by allele type
bin_vcf_long <- bin_vcf %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(MGfile, by = "ID") %>%
  tidyr::gather(Ind, key, 2:(base::ncol(.)-4))

#Calculate phenotypic association of each allele type for each SNP
preVarfile <- bin_vcf_long %>%
  dplyr::left_join(pheno, by = "Ind") %>%
  dplyr::group_by(ID, MGs, key) %>%
  dplyr::summarize(nInd = dplyr::n(),avPheno=base::mean(Pheno, na.rm = T), .groups = 'keep')

types <- c(ref = '0', het = '1', alt = '2', miss = '<NA>')

noNA_preVarfile <- preVarfile %>%
  dplyr::select(-avPheno) %>%
  tidyr::spread(key, nInd) %>%
  dplyr::rename(dplyr::any_of(types))

if(!("miss" %in% colnames(noNA_preVarfile))){
  noNA_preVarfile <- dplyr::mutate(noNA_preVarfile, miss = 0)
}

if(!("het" %in% colnames(noNA_preVarfile))){
  noNA_preVarfile <- dplyr::mutate(noNA_preVarfile, het = 0)
}

noNA_preVarfile$MGs[is.na(noNA_preVarfile$MGs)] <- "0"
noNA_preVarfile[is.na(noNA_preVarfile)] <- 0

Varfile <- preVarfile %>% dplyr::select(-nInd) %>%
            tidyr::spread(key, avPheno) %>%
            dplyr::rename(dplyr::any_of(types)) %>%
            dplyr::mutate(percdiff = alt - ref) %>%
            dplyr::ungroup() %>%
            dplyr::select(ID, percdiff) %>%
            dplyr::left_join(noNA_preVarfile %>%
                              dplyr::mutate(AltAF = (2*alt+het)/(2*(ref + het + alt))),
                             by = c("ID"))
return(Varfile)
}

