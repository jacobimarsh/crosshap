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
#' @param minHap Minimum nIndividuals to keep haplotype combinations.
#' @param hetmiss_as Treat heterozygous missing './X' as missing or allele.
#' @param metadata Metadata input
#'
#' @export
#'
#'
run_haplotyping <- function(vcf, LD, pheno, epsilon = c(0.4,0.8,1.2,1.6,2), MGmin, minHap = 5, hetmiss_as = 'allele', metadata = NULL, keep_outliers = F) {
    #Reformat VCF

  cli::cli_progress_bar(total = 4*length(epsilon) + 2
                   )
  step <- "Formatting VCF"

  cli::cli_progress_update()

  bin_vcf <- dplyr::select(vcf, -c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
  dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,
                                             base::ifelse(x=='1|0'|x=='0|1',1,
                                                          base::ifelse(x=='1|1',2,
                                                                       switch(hetmiss_as, "allele" = base::ifelse(x=='1|.'|x=='.|1',1,
                                                                                                                 base::ifelse(x=='0|.'|x=='.|0',0,NA)),
                                                                                         "miss" = NA))))})
  cli::cli_progress_update()
  for (arez in epsilon){

    #Run DBscan on LD matrix
    # base::message(paste0("Clustering SNPs into marker groups (eps = ",arez,")"))
    step <- paste0("eps(",arez,") Clustering SNPs into marker groups")
    cli::cli_progress_update()

    db40 <- dbscan::dbscan(LD, eps = arez, minPts = MGmin)
    preMGfile <- tibble::tibble(ID=rownames(LD),cluster=db40$cluster) %>%
      dplyr::left_join(dplyr::select(vcf, 2:3), by = "ID")


    ##Identify haplotype frequencies for different marker group combinations
    # base::message(paste0("Determining haplotypes from marker group clusters (eps = ",arez,")"))
    step <- paste0("eps(",arez,") Determining haplotypes from marker group clusters")
    cli::cli_progress_update()

    phaps_out <- pseudo_haps(preMGfile = preMGfile, bin_vcf = bin_vcf, minHap = minHap, LD = LD, keep_outliers = keep_outliers)

    ##Build summary object with all relevant haplotyping information
    # base::message(paste0("Collating haplotype information (eps = ",arez,")"))
    step <- paste0("eps(",arez,") Collating haplotype information")
    cli::cli_progress_update()

    Varfile <- tagphenos(MGfile = phaps_out$MGfile, bin_vcf, pheno)
    clustered_hpS_obj <-  base::list(epsilon = db40$eps,
                               MGmin = db40$minPts,
                               Hapfile = phaps_out$Hapfile,
                               Indfile = if (missing(metadata)|is.null(metadata)){
                                 dplyr::left_join(phaps_out$nophenIndfile, pheno, by = "Ind") %>% dplyr::mutate(Metadata = as.character(NA))
                               }else {
                                 dplyr::left_join(phaps_out$nophenIndfile, pheno, by = "Ind") %>% dplyr::left_join(metadata, by = "Ind")
                               },
                               Varfile = Varfile,
                               MGfile = phaps_out$MGfile)
                               #ClusterR2)
    base::assign(paste("Haplotypes_MGmin",MGmin, "_E", arez,sep = ""), clustered_hpS_obj, envir = .GlobalEnv)
    cli::cli_progress_update()
  }
  cli::cli_alert_success(paste0("Haplotyping complete!"))
  base::message(paste0("Info saved in Haplotypes_",MGmin, "_E objects"))
}

