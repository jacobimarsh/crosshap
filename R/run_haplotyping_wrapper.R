#' Cluster SNPs and identify haplotypes
#'
#' Performs density-based clustering of SNPs in region of interest to identify
#' marker groups. Defines haplotype combinations and classifies individuals
#' based on characteristic combinations of marker group alleles. Returns a
#' comprehensive haplotyping object, required to build clustering tree and
#' visualize haplotypes.
#'
#' @param vcf Input VCF for region of interest.
#' @param LD Pairwise correlation matrix of SNPs in region from PLINK.
#' @param pheno Input numeric phenotype data for each individual.
#' @param epsilon Epsilon values for clustering SNPs with DBscan.
#' @param MGmin Minimum SNPs in marker groups, MinPts parameter for DBscan.
#'
#' @return
#' @export
#'
#' @example
#' run_haplotyping(vcf, LD, phen_early, epsilon, MGmin)
#'
run_haplotyping <- function(vcf, LD, pheno, epsilon, MGmin) {
  for (arez in epsilon){
    #Run DBscan on LD matrix
    base::message(paste0("Clustering SNPs into marker groups (eps = ",arez,")"))
    db40 <- dbscan::dbscan(LD, eps = arez, minPts = MGmin)
    MGfile <- tibble::tibble(POS=rownames(LD),cluster=db40$cluster)

    ##Call allelic states for each SNP marker group across individuals
    base::message(paste0("Classifying marker group alleles in individuals (eps = ",arez,")"))
    het_pseudoSNP <- create_pseudoSNP(MGfile = MGfile, vcf = vcf)

    ##Identify haplotype frequencies for different marker group combinations
    base::message(paste0("Calculating haplotype combination frequencies (eps = ",arez,")"))
    clustered_hpS <- pseudo2haps(het_pseudoSNP)

    ##Build summary object with all relevant haplotyping information
    base::message(paste0("Collating haplotype information (eps = ",arez,")"))
    Varfile <- tagphenos(MGfile, vcf, pheno)
    clustered_hpS_obj <-  base::list(epsilon = db40$eps,
                               MGmin = db40$minPts,
                               Hapfile = clustered_hpS$Hapfile,
                               IDfile = dplyr::left_join(clustered_hpS$nophenIDfile,pheno, by = "ID"),
                               Varfile = Varfile,
                               MGfile = MGfile)
    base::assign(paste("Haplotypes_MP", "_E", arez,sep = ""), clustered_hpS_obj, envir = .GlobalEnv)
    base::message(paste0("Done! (object saved as Haplotypes_MP",MGmin,"_E",arez,")",sep = ""))
  }
}
