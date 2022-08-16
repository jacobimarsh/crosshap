#' Call allelic states for each SNP marker group across individuals
#'
#' Have a drescription here.
#' fjdklsa;jfkdl;sa
#' fjdskal;jfkdls;a
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
run_haplotyping <- function(vcf, LD, epsilon, pheno) {
  ##Perform clustering and count haplotype frequencies

  for (arez in epsilon){
    message(paste0("Clustering SNPs into marker groups (eps = ",arez,")"))
    db40 <- dbscan::dbscan(LD, eps = arez, minPts = 40)
    MGfile <- tibble(POS=rownames(LD),cluster=db40$cluster)

    ##Call allelic states for each SNP marker group across individuals
    message(paste0("Classifying marker group alleles in individuals (eps = ",arez,")"))
    het_pseudoSNP <- create_pseudoSNP(MGfile = MGfile, vcf = vcf)

    ##Identify haplotype frequencies for different marker group combinations
    message(paste0("Calculating haplotype combination frequencies (eps = ",arez,")"))
    clustered_hpS <- pseudo2haps(het_pseudoSNP)

    ##Build summary object with all relevant haplotyping information
    message(paste0("Collating haplotype information (eps = ",arez,")"))
    Varfile <- tagphenos(MGfile, vcf, pheno)
    clustered_hpS_obj <-  list(epsilon = db40$eps,
                               minPts = db40$minPts,
                               Hapfile = clustered_hpS$Hapfile,
                               IDfile = left_join(clustered_hpS$nophenIDfile,pheno, by = "ID"),
                               Varfile = Varfile,
                               MGfile = MGfile)
    assign(paste("Haplotypes_MP", "_E", arez,sep = ""), clustered_hpS_obj, envir = .GlobalEnv)
    message(paste0("Done (eps = ",arez"!"))
  }
}
