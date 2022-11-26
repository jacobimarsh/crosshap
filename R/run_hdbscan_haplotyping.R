#' Cluster SNPs with HDBSCAN and identify haplotypes
#'
#' run_hdbscan_haplotyping() performs HDBSCAN clustering of SNPs in region of
#' interest to identify marker groups. Individuals are classified by haplotype
#' combination based on shared combinations of marker group alleles. Returns a
#' comprehensive haplotyping object (HapObject), which can be used as input to
#' build clustering tree for epsilon optimization using run_clustree(), and
#' can be visualized with reference to phenotype and metadata using
#' crosshap_viz().
#'
#' @param vcf Input VCF for region of interest.
#' @param LD Pairwise correlation matrix of SNPs in region (e.g. from PLINK).
#' @param pheno Input numeric phenotype data for each individual.
#' @param epsilon Epsilon values for clustering SNPs with DBscan.
#' @param MGmin Minimum SNPs in marker groups, MinPts parameter for DBscan.
#' @param minHap Minimum nIndividuals in a haplotype combination.
#' @param hetmiss_as If hetmiss_as = "allele", heterozygous-missing SNPs './N'
#' are recoded as 'N/N', if hetmiss_as = "miss", the site is recoded as missing.
#' @param metadata Metadata input (optional).
#' @param keep_outliers When FALSE, marker group smoothing is performed to
#' remove outliers.
#'
#' @export
#'
#' @returns A comprehensive haplotyping S3 object (HapObject) for each provided
#' epsilon value, needed for run_clustree() and crosshap_viz().
#'
#' @examples
#'
#' if (FALSE) {
#'     run_hdbscan_haplotyping(vcf, LD, pheno, metadata, MGmin = 30)
#'}
#'

run_hdbscan_haplotyping <- function(vcf, LD, pheno, MGmin, minHap = 5, hetmiss_as = 'allele', metadata = NULL, keep_outliers = F){
  bin_vcf <- dplyr::select(vcf, -c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
    dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,
                                               base::ifelse(x=='1|0'|x=='0|1',1,
                                                            base::ifelse(x=='1|1',2,
                                                                         switch(hetmiss_as, "allele" = base::ifelse(x=='1|.'|x=='.|1',1,
                                                                                                                    base::ifelse(x=='0|.'|x=='.|0',0,NA)),
                                                                                "miss" = NA))))})

    #Run HDBscan on LD matrix
     base::message(paste0("Clustering SNPs into marker groups"))

    db40 <- dbscan::hdbscan(LD, minPts = MGmin)

    preMGfile <- tibble::tibble(ID=rownames(LD),cluster=db40$cluster) %>%
      dplyr::left_join(dplyr::select(vcf, 2:3), by = "ID")


    ##Identify haplotype frequencies for different marker group combinations
    base::message(paste0("Determining haplotypes from marker group clusters"))

    phaps_out <- pseudo_haps(preMGfile = preMGfile, bin_vcf = bin_vcf, minHap = minHap, LD = LD, keep_outliers = keep_outliers)

    ##Build summary object with all relevant haplotyping information
    base::message(paste0("Collating haplotype information"))

    Varfile <- tagphenos(MGfile = phaps_out$MGfile, bin_vcf, pheno)
    clustered_hpS_obj <-  base::list(MGmin = db40$minPts,
                                     Hapfile = phaps_out$Hapfile,
                                     Indfile = if(missing(metadata)|is.null(metadata)){
                                       dplyr::left_join(phaps_out$nophenIndfile, pheno, by = "Ind") %>% dplyr::mutate(Metadata = as.character(NA))
                                     }else {
                                       dplyr::left_join(phaps_out$nophenIndfile, pheno, by = "Ind") %>% dplyr::left_join(metadata, by = "Ind")
                                     },
                                     Varfile = Varfile,
                                     MGfile = phaps_out$MGfile)
    base::assign(paste("Haplotypes_MGmin",MGmin, "_HDBSCAN",sep = ""), clustered_hpS_obj, envir = .GlobalEnv)
  base::message(paste0("Info saved in Haplotypes_",MGmin, "_E_HDBSCAN"))
}

