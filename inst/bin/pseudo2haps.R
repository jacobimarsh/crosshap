#' Classify haplotype combinations and populations
#'
#' Internal function that reads pseudoSNP data to identify and call the
#' frequency for characteristic marker group combinations. Defines common marker
#' group allele combinations as haplotypes and returns the haplotype of each
#' individual.
#'
#' @param pSNP 'pseudoSNPs' created by pseudo.R
#' @param minHap Minimum size (nIndividuals) to keep haplotype combinations
#'
#' @return
#' @export
#'
#' @examples
#' pseudo2haps(het_pseudoSNP, 9)
#'
pseudo2haps <- function(pSNP, minHap) {
  cnames <- base::colnames(dplyr::select(pSNP, -Ind))

  het_hapCounts <- pSNP %>%
  tidyr::gather(mgs,value,2:base::ncol(.)) %>%
  dplyr::group_by(Ind) %>%
  dplyr::mutate(mgs_new=base::paste(value, collapse = '_')) %>%
  dplyr::group_by(mgs_new) %>%
  dplyr::tally() %>% dplyr::mutate(n=n/(base::ncol(pSNP)-1)) %>%
  tidyr::separate(col = mgs_new, into = cnames, sep = "_") %>%
  tibble::as_tibble() %>%
  dplyr::arrange(by_group = -n)

  over20_hhCounts <- dplyr::filter(het_hapCounts, n > minHap) %>%
  dplyr::mutate(hap=LETTERS[1:base::nrow(.)])

  base::suppressMessages(clustered_hpS <- dplyr::left_join(pSNP, over20_hhCounts) %>%
  dplyr::mutate_if(is.character,function(x){tidyr::replace_na(x, '0')}) %>%
  dplyr::select(1,base::ncol(.)))

  return(base::list(Hapfile = over20_hhCounts, nophenIndfile = clustered_hpS))
}

