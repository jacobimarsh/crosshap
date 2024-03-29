#' Identify haplotypes from clustered SNPs
#'
#'
#' pseudo_haps() calls the most common allelic states for each SNP marker group
#' across individuals, before building dummy SNPs for each marker group that
#' mimic the binary vcf format. This is the step which determines the haplotype
#' combinations, and therefore enables several summaries to be returned - as
#' contained in the $Hapfile and preliminary $Indfile and finalised $MGfile,
#' following marker group smoothing. This is an internal function not intended
#' for external use.
#'
#'
#' @param preMGfile SNP clusters from DBscan.
#' @param bin_vcf Binary VCF for region of interest reformatted by
#' run_haplotyping().
#' @param minHap Minimum size (nIndividuals) to keep haplotype combinations
#' @param LD LD matrix input.
#' @param keep_outliers When FALSE, marker group smoothing is performed to
#' remove outliers.
#'
#' @importFrom rlang ".data"
#' @importFrom rlang ":="
#'
#' @export
#'
#' @return Returns intermediate of haplotype object
#'

pseudo_haps <- function(preMGfile, bin_vcf, minHap, LD, keep_outliers) {

##Call allelic states for each SNP marker group across individuals
#Extract SNPs in first MG cluster (MG1)
  dbscan_cvel <- preMGfile %>%
    dplyr::filter(.data$cluster == 1) %>% tibble::as_tibble()
  c1_vcf <- bin_vcf %>%
    tibble::rownames_to_column() %>% dplyr::filter(.data$rowname %in% dbscan_cvel$ID) %>% tibble::column_to_rownames()

#Calculate most common (alternate or ref) allele for MG1 across all individuals
  i_modes_1 <- base::apply(c1_vcf %>% base::sapply(as.double), 2, crosshap::arith_mode) %>% tibble::as_tibble() %>%
    dplyr::pull(.data$value)

#Build a dummy VCF with one pseudoSNP position (MG1)
  pseudoSNP <- tibble::tibble(Ind = base::colnames(c1_vcf), "1" = i_modes_1)

#Repeat for all other MGs, iteratively adding to pseudoSNP bin_vcf
  if(max(preMGfile$cluster) > 1){
for (vel in c(2:base::max(preMGfile$cluster))) {

    dbscan_cvel <- preMGfile %>%
      dplyr::filter(.data$cluster == vel) %>% tibble::as_tibble()
    cvel_vcf <- bin_vcf %>%
      tibble::rownames_to_column() %>% dplyr::filter(.data$rowname %in% dbscan_cvel$ID) %>% tibble::column_to_rownames()
    i_modes_vel <-
      base::apply(cvel_vcf %>% base::sapply(as.double), 2, crosshap::arith_mode) %>% tibble::as_tibble()
    pseudoSNP <- dplyr::bind_cols(pseudoSNP, i_modes_vel)
    base::colnames(pseudoSNP)[vel + 1] <- base::paste0(vel)
}
}

#Convert heterozygous marker groups to homozygous alternate (optional)
    pseudoSNP <- pseudoSNP %>%
      dplyr::mutate_all((as.character))

##Identify haplotype frequencies from different marker group combinations

  cnames <- base::colnames(dplyr::select(pseudoSNP, -"Ind"))

  hapCounts <- pseudoSNP %>%
    tidyr::gather('mgs','value',2:base::ncol(pseudoSNP)) %>%
    dplyr::group_by(.data$Ind) %>%
    dplyr::mutate(mgs_new=base::paste(.data$value, collapse = '_')) %>%
    dplyr::group_by(.data$mgs_new) %>%
    dplyr::tally() %>% dplyr::mutate(n=.data$n/(base::ncol(pseudoSNP)-1)) %>%
    tidyr::separate(col = .data$mgs_new, into = cnames, sep = "_") %>%
    tibble::as_tibble() %>%
    dplyr::arrange(by_group = -.data$n) %>%
    dplyr::filter(!dplyr::if_any(dplyr::everything(), ~ base::grepl('NA', .)))

#Remove haplotypes with frequency below minHap and label remaining with letters
#Remove any rows with
  overmin_hap_counts <- dplyr::filter(hapCounts, .data$n > minHap) %>%
    dplyr::mutate(hap=LETTERS[1:base::nrow(dplyr::filter(hapCounts, .data$n > minHap))])

#Report what each individual's haplotype letter is based on their pseudoSNP combination
  base::suppressMessages(
    nophenIndfile <- dplyr::left_join(pseudoSNP, overmin_hap_counts, by = c(as.character(1:max(dbscan_cvel$cluster)))) %>%
                           dplyr::mutate_if(is.character,function(x){tidyr::replace_na(x, '0')}) %>%
                           dplyr::select(1,(base::ncol(pseudoSNP)+2))
    )

#Change clusters to as.numeric
   overmin_hap_counts[,1:(base::ncol(overmin_hap_counts)-1)] <-
    base::sapply(overmin_hap_counts[, 1:(base::ncol(overmin_hap_counts)-1)], as.numeric)

#Filter out clusters that don't have an alternate allele in remaining haplotypes
   no0clust <- overmin_hap_counts %>% dplyr::select(dplyr::where(~ base::is.numeric(.x) && base::sum(.x) != 0))

#Reformat as nice hapfile table
   Hapfile <- base::cbind(overmin_hap_counts$hap,no0clust[,base::names(base::sort(base::colSums(no0clust), decreasing = TRUE))]) %>%
     dplyr::rename('hap' = "overmin_hap_counts$hap")

   clust_preMGs <- base::colnames(Hapfile %>% dplyr::select(-c("hap", "n")))

#Add 'MG' prefix to clusters
   for (veal in c(1:(base::length(clust_preMGs)))) {
     base::colnames(Hapfile)[veal+2] <- base::paste0("MG",veal)
     }

#Report which clusters labels from DBSCAN correspond to which final MGs
   cluster2MGs <- base::cbind(clust_preMGs, base::colnames(Hapfile[3:base::ncol(Hapfile)])) %>%
     magrittr::set_colnames(c("cluster", "MGs")) %>%
     tibble::as_tibble() %>%
     dplyr::mutate(cluster = as.numeric(.data$cluster))

#Save MGfile created at this step before smoothing
  unsmoothed_MGfile <- dplyr::left_join(preMGfile, cluster2MGs, by = "cluster")

  unsmoothed_MGfile[is.na(unsmoothed_MGfile)] <- "0"

#Prepare blank tibble for r2 files
  r2prefile <- tibble::tibble(ID = character(), premeanr2 = double(), MGs = character())

#Calculate average R2 of each SNP with other SNPs within same MG and add to r2file
  for (grev in unique(unsmoothed_MGfile$MGs)){
    r2prefile <-  r2prefile %>% rbind(tibble::enframe(colMeans((LD %>%
                                                        dplyr::filter(row.names(LD) %in% dplyr::filter(unsmoothed_MGfile, .data$MGs == grev)$ID))[,dplyr::filter(unsmoothed_MGfile, .data$MGs == grev)$ID])) %>%
                              dplyr::rename('ID' = "name", 'premeanr2' = "value"))
  }

#Add premeanr2 results to unsmoothed MGfile
  unsmoothed_MGfile <- dplyr::left_join(unsmoothed_MGfile, r2prefile, by = "ID")

#Remove SNPs from MGs when they are beyond 2 std devs from median premeanr2 of a given MG
  smoothed_MGfile <- unsmoothed_MGfile %>% dplyr::group_by(.data$MGs) %>%
    dplyr::mutate(MGs = ifelse((abs(.data$premeanr2 - stats::median(.data$premeanr2)) > 2*stats::sd(.data$premeanr2)),0, .data$MGs))

#Re-calculate meanR2 as before, now that outlier SNPs have been trimmed from MGs
  r2file <- tibble::tibble(ID = character(), meanr2 = double(), MGs = character())

  for (grev in unique(smoothed_MGfile$MGs)){
    r2file <-  r2file %>% rbind(tibble::enframe(colMeans((LD %>%
                                                                  dplyr::filter(row.names(LD) %in% dplyr::filter(smoothed_MGfile, .data$MGs == grev)$ID))[,dplyr::filter(smoothed_MGfile, .data$MGs == grev)$ID])) %>%
                                        dplyr::rename(ID = "name", meanr2 = "value"))
  }

#Add information to MGfile ready for export
  MGfile <- dplyr::left_join(smoothed_MGfile, r2file, by = "ID")


  if(keep_outliers == T){
    return(base::list(Hapfile = Hapfile, nophenIndfile = nophenIndfile, MGfile = unsmoothed_MGfile %>% dplyr::rename('meanr2' = 'premeanr2')))
  } else {
    return(base::list(Hapfile = Hapfile, nophenIndfile = nophenIndfile, MGfile = MGfile))
  }
}

