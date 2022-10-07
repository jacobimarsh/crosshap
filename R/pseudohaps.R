#' Identify haplotypes from SNP clusters
#'
#' Internal function that calls the most common allelic states for each SNP
#' marker group across individuals, returning dummy SNPs for each marker group
#' mimicking the binary vcf format. Options to be added for treatment of heterozygotes.
#'
#' @param preMGfile SNP clusters from DBscan.
#' @param bin_vcf Input binary VCF for region of interest.
#' @param minHap Minimum size (nIndividuals) to keep haplotype combinations
#' @param LD LD matrix input
#' @param keep_outliers whether to keep outlier SNPs in MGs
#'
#' @export
#'
#'
pseudo_haps <- function(preMGfile, bin_vcf, minHap, LD, keep_outliers) {

##Call allelic states for each SNP marker group across individuals
#Extract SNPs in first MG cluster (MG1)
  db40_c1 <- preMGfile %>%
    dplyr::filter(cluster == 1) %>% tibble::as_tibble()
  c1_vcf <- bin_vcf %>%
    tibble::rownames_to_column() %>% dplyr::filter(rowname %in% db40_c1$ID) %>% tibble::column_to_rownames()

#Calculate most common (alternate or ref) allele for MG1 across all individuals
  i_modes_1 <- base::apply(c1_vcf %>% base::sapply(as.double), 2, crosshap::mode) %>% tibble::as_tibble() %>%
    dplyr::pull(value)

#Build a dummy VCF with one pseudoSNP position (MG1)
  pseudoSNP <- tibble::tibble(Ind = base::colnames(c1_vcf), "1" = i_modes_1)

#Repeat for all other MGs, iteratively adding to pseudoSNP bin_vcf
for (vel in c(2:base::max(preMGfile$cluster))) {

    db40_cvel <- preMGfile %>%
      dplyr::filter(cluster == vel) %>% tibble::as_tibble()
    cvel_vcf <- bin_vcf %>%
      tibble::rownames_to_column() %>% dplyr::filter(rowname %in% db40_cvel$ID) %>% tibble::column_to_rownames()
    i_modes_vel <-
      base::apply(cvel_vcf %>% base::sapply(as.double), 2, mode) %>% tibble::as_tibble()
    pseudoSNP <- dplyr::bind_cols(pseudoSNP, i_modes_vel)
    base::colnames(pseudoSNP)[vel + 1] <- base::paste0(vel)
  }

#Convert heterozygous marker groups to homozygous alternate (optional)
  het_pseudoSNP <- pseudoSNP %>%
    dplyr::mutate_if(is.numeric,function(x) {base::gsub(1, 2, x,fixed = T, )})

##Identify haplotype frequencies from different marker group combinations

  cnames <- base::colnames(dplyr::select(het_pseudoSNP, -Ind))

  het_hapCounts <- het_pseudoSNP %>%
    tidyr::gather(mgs,value,2:base::ncol(.)) %>%
    dplyr::group_by(Ind) %>%
    dplyr::mutate(mgs_new=base::paste(value, collapse = '_')) %>%
    dplyr::group_by(mgs_new) %>%
    dplyr::tally() %>% dplyr::mutate(n=n/(base::ncol(het_pseudoSNP)-1)) %>%
    tidyr::separate(col = mgs_new, into = cnames, sep = "_") %>%
    tibble::as_tibble() %>%
    dplyr::arrange(by_group = -n)

  over20_hhCounts <- dplyr::filter(het_hapCounts, n > minHap) %>%
    dplyr::mutate(hap=LETTERS[1:base::nrow(.)])

  base::suppressMessages(
    clustered_hpS <- dplyr::left_join(het_pseudoSNP, over20_hhCounts, by = c(as.character(1:max(db40_cvel$cluster)))) %>%
                           dplyr::mutate_if(is.character,function(x){tidyr::replace_na(x, '0')}) %>%
                           dplyr::select(1,base::ncol(.))
    )

#Change clusters to as.numeric
   over20_hhCounts[,1:(base::ncol(over20_hhCounts)-1)] <-
    base::sapply(over20_hhCounts[, 1:(base::ncol(over20_hhCounts)-1)], as.numeric)

   no0clust <- over20_hhCounts %>% dplyr::select(where(~ base::is.numeric(.x) && base::sum(.x) != 0))

   dat1 <- base::cbind(over20_hhCounts$hap,no0clust[,base::names(base::sort(base::colSums(no0clust), decreasing = TRUE))]) %>%
     dplyr::rename(hap = "over20_hhCounts$hap")

   clust_preMGs <- base::colnames(dat1 %>% dplyr::select(-c(hap, n)))

   for (veal in c(1:(base::length(clust_preMGs)))) {
     base::colnames(dat1)[veal+2] <- base::paste0("MG",veal)
     }

   cluster2MGs <- base::cbind(clust_preMGs, base::colnames(dat1[3:base::ncol(dat1)])) %>%
     magrittr::set_colnames(c("cluster", "MGs")) %>%
     tibble::as_tibble() %>%
     dplyr::mutate(cluster = as.numeric(cluster))

  unsmoothed_MGfile <- dplyr::left_join(preMGfile, cluster2MGs, by = "cluster")

  unsmoothed_MGfile[is.na(unsmoothed_MGfile)] <- "0"

  r2prefile <- tibble::tibble(ID = character(), premeanr2 = double(), MGs = character())

  for (grev in unique(unsmoothed_MGfile$MGs)){
    r2prefile <-  r2prefile %>% rbind(tibble::enframe(colMeans((LD %>%
                                                        dplyr::filter(row.names(LD) %in% dplyr::filter(unsmoothed_MGfile, MGs == grev)$ID))[,dplyr::filter(unsmoothed_MGfile, MGs == grev)$ID])) %>%
                              dplyr::rename(ID = "name", premeanr2 = "value"))
  }

  unsmoothed_MGfile <- dplyr::left_join(unsmoothed_MGfile, r2prefile, by = "ID")

  smoothed_MGfile <- unsmoothed_MGfile %>% group_by(MGs) %>%
    dplyr::mutate(MGs = ifelse((abs(premeanr2 - median(premeanr2)) > 2*sd(premeanr2)),0, MGs))

  r2file <- tibble::tibble(ID = character(), meanr2 = double(), MGs = character())

  for (grev in unique(smoothed_MGfile$MGs)){
    r2file <-  r2file %>% rbind(tibble::enframe(colMeans((LD %>%
                                                                  dplyr::filter(row.names(LD) %in% dplyr::filter(smoothed_MGfile, MGs == grev)$ID))[,dplyr::filter(smoothed_MGfile, MGs == grev)$ID])) %>%
                                        dplyr::rename(ID = "name", meanr2 = "value"))
  }

  MGfile <- dplyr::left_join(smoothed_MGfile, r2file, by = "ID")


  if(keep_outliers == T){
    return(base::list(Hapfile = dat1, nophenIndfile = clustered_hpS, MGfile = unsmoothed_MGfile %>% dplyr::rename(meanr2 = premeanr2)))
  } else {
    return(base::list(Hapfile = dat1, nophenIndfile = clustered_hpS, MGfile = MGfile))
  }
}

