#' Identify haplotype frequencies for different marker group combinations
#'
#' @param elon Epsilon values for DBscan.
#' @param vcf Input VCF for region of interest.
#' @param LD Pairwise R squared correlation matrix between all SNP loci.
#'
#' @return
#' @export
#'
#' @examples
#'
pseudo2haps <- function(pSNP) {
  cnames <- colnames(select(pSNP, -ID))

  het_hapCounts <- pSNP %>%
  gather(mgs,value,2:ncol(.)) %>%
  group_by(ID) %>%
  mutate(mgs_new=paste(value, collapse = '_')) %>%
  group_by(mgs_new) %>%
  tally() %>% mutate(n=n/(ncol(pSNP)-1)) %>%
  separate(col = mgs_new, into = cnames, sep = "_") %>%
  as_tibble() %>%
  arrange(by_group = -n)

  over20_hhCounts <- dplyr::filter(het_hapCounts, n > 9) %>%
  mutate(hap_eps=LETTERS[1:nrow(.)]) %>%
  rename(!!paste0("hap_eps",arez) := 'hap_eps')

  clustered_hpS <- left_join(pSNP, over20_hhCounts) %>%
  mutate_if(is.character,function(x){replace_na(x, '0')}) %>%
  select(1,ncol(.))

  return(clustered_hpS)
}
