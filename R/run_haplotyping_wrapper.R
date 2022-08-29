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
#' @param minHap Minimum nIndividuals to keep haplotype combinations
#'
#' @return
#' @export
#'
#' @example
#' run_haplotyping(vcf, LD, phen_early, epsilon, MGmin, minHap)
#'
run_haplotyping <- function(vcf, LD, pheno, epsilon = c(0.5,1,1.5,2,2.5,3), MGmin, minHap) {
    #Reformat VCF
  bin_vcf <- dplyr::select(vcf, -c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
  dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,base::ifelse(x=='1|0'|x=='0|1',1,base::ifelse(x=='1|1',2,'failsave')))})

  for (arez in epsilon){
    #Run DBscan on LD matrix
    base::message(paste0("Clustering SNPs into marker groups (eps = ",arez,")"))
    db40 <- dbscan::dbscan(LD, eps = arez, minPts = MGmin)
    MGfile <- tibble::tibble(ID=rownames(LD),cluster=db40$cluster) %>%
      dplyr::left_join(dplyr::select(vcf, 2:3), by = "ID")

    ##Call allelic states for each SNP marker group across individuals
    base::message(paste0("Classifying marker group alleles in individuals (eps = ",arez,")"))
    het_pseudoSNP <- crosshap::create_pseudoSNP(MGfile = MGfile, bin_vcf = bin_vcf)

    ##Identify haplotype frequencies for different marker group combinations
    base::message(paste0("Calculating haplotype combination frequencies (eps = ",arez,")"))
    clustered_hpS <- crosshap::pseudo2haps(het_pseudoSNP, minHap)

    ##Build summary object with all relevant haplotyping information
    base::message(paste0("Collating haplotype information (eps = ",arez,")"))
    Varfile <- crosshap::tagphenos(MGfile, bin_vcf, pheno)
    clustered_hpS_obj <-  base::list(epsilon = db40$eps,
                               MGmin = db40$minPts,
                               Hapfile = clustered_hpS$Hapfile,
                               Indfile = dplyr::left_join(clustered_hpS$nophenIndfile,pheno, by = "Ind"),
                               Varfile = Varfile,
                               MGfile = MGfile)
    base::assign(paste("Haplotypes_MGmin",MGmin, "_E", arez,sep = ""), clustered_hpS_obj, envir = .GlobalEnv)
    base::message(paste0("Done! (object saved as Haplotypes_MGmin",MGmin,"_E",arez,")",sep = ""))
  }
}
