#' Cluster SNPs and identify haplotypes
#'
#' run_haplotyping() performs density-based clustering of SNPs in region of
#' interest to identify Marker Groups. Individuals are classified by haplotype
#' combination based on shared combinations of Marker Group alleles. Returns a
#' haplotyping object (HapObject), which can be used as input to build clustering
#' tree for epsilon optimization using clustree_viz(), and can be visualized with
#' reference to phenotype and metadata using crosshap_viz().
#'
#' @param vcf Input VCF for region of interest.
#' @param LD Pairwise correlation matrix of SNPs in region (e.g. from PLINK).
#' @param metadata Metadata input (optional).
#' @param pheno Input numeric phenotype data for each individual.
#' @param epsilon Epsilon values for clustering SNPs with DBscan.
#' @param MGmin Minimum SNPs in marker groups, MinPts parameter for DBscan.
#' @param minHap Minimum nIndividuals in a haplotype combination.
#' @param hetmiss_as If hetmiss_as = "allele", heterozygous-missing SNPs './N'
#' are recoded as 'N/N', if hetmiss_as = "miss", the site is recoded as missing.
#' @param het_phenos When FALSE, phenotype associations for SNPs are calculated
#' from reference and alternate allele individuals only, when TRUE, heterozygous
#' individuals are included assuming additive allele effects.
#' @param keep_outliers When FALSE, marker group smoothing is performed to
#' remove outliers.
#'
#' @export
#'
#' @returns A comprehensive haplotyping S3 object (HapObject) for each provided
#' epsilon value, needed for clustree_viz() and crosshap_viz().
#'

run_haplotyping <- function(vcf, LD, pheno, metadata = NULL,
                            epsilon = c(0.2,0.4,0.6,0.8,1),
                            MGmin = 30, minHap = 9, hetmiss_as = 'allele',
                            het_phenos = FALSE, keep_outliers = FALSE){
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
  HapObject <- NULL
  for (arez in epsilon){

    #Run DBscan on LD matrix
    step <- paste0("eps(",arez,") Clustering SNPs into marker groups")
    cli::cli_progress_update()

    db_out <- dbscan::dbscan(LD, eps = arez, minPts = MGmin)
    preMGfile <- tibble::tibble(ID=rownames(LD),cluster=db_out$cluster) %>%
      dplyr::left_join(dplyr::select(vcf, 2:3), by = "ID")


    ##Identify haplotype frequencies for different marker group combinations
    step <- paste0("eps(",arez,") Determining haplotypes from marker group clusters")
    cli::cli_progress_update()
if(length(unique(preMGfile$cluster)) != 1){
    phaps_out <- pseudo_haps(preMGfile = preMGfile, bin_vcf = bin_vcf, minHap = minHap, LD = LD, keep_outliers = keep_outliers)

    ##Build summary object with all relevant haplotyping information
    step <- paste0("eps(",arez,") Collating haplotype information")
    cli::cli_progress_update()

    Varfile <- tagphenos(MGfile = phaps_out$MGfile, bin_vcf, pheno, het_phenos = het_phenos)
    clustered_hap_obj <-  base::list(epsilon = db_out$eps,
                               MGmin = db_out$minPts,
                               Hapfile = phaps_out$Hapfile,
                               Indfile = if (missing(metadata)|is.null(metadata)){
                                 dplyr::left_join(phaps_out$nophenIndfile, pheno, by = "Ind") %>% dplyr::mutate(Metadata = as.character(NA))
                               }else {
                                 dplyr::left_join(phaps_out$nophenIndfile, pheno, by = "Ind") %>% dplyr::left_join(metadata, by = "Ind")
                               },
                               Varfile = Varfile)
    HapObject[[paste("Haplotypes_MGmin",MGmin, "_E", arez,sep = "")]] <- clustered_hap_obj
    } else {
      cli::cli_progress_update()
      base::message(paste0("No Marker Groups identified for Epsilon = ",arez, ""))

      }
    cli::cli_progress_update()
  }

  #Warn users if het_phenos = F and there are heterozygous sites, or phenotype assoc scores not calculated
  if(sum(is.na(HapObject[[1]]$Varfile$phenodiff)) > 0){
    base::message(paste0("WARNING: Phenotype association scores (phenodiff) not calculated for ",
                         sum(is.na(HapObject[[1]]$Varfile$phenodiff)), " of ",length(HapObject[[1]]$Varfile$phenodiff),
                         " sites"))
  }

  if(het_phenos == FALSE && (sum(HapObject[[1]]$Varfile$alt < HapObject[[1]]$Varfile$het) + sum(Varfile$ref < HapObject[[1]]$Varfile$het) > 0)){
    base::message(paste0("NOTE: ",
                         sum(HapObject[[1]]$Varfile$alt < HapObject[[1]]$Varfile$het),
                         " sites have more heterozygous individuals (1/0) than one of the homozygous states (1/1 or 0/0).

Haplotyping was performed with hetphenos = F, meaning phenotype association scores for heterozygous sites was ignored."))
  }

list(paste("Haplotypes_MGmin",MGmin, "_E", epsilon,sep = ""))

  cli::cli_alert_success(paste0("Haplotyping complete!"))
  base::message(paste0("Info saved in Haplotypes_",MGmin, "_E objects"))
return(HapObject)
  }

