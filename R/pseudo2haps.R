#' Run Clustering on Data
#'
#' @param elon Epsilon values to search.
#' @param vcf Input vcf of region.
#' @param LD Pairwise R squared correlation matrix.
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

over20_hhCounts <- filter(het_hapCounts, n > 9) %>%
  mutate(hap_eps=LETTERS[1:nrow(.)]) %>%
  rename(!!paste0("hap_eps",arez) := 'hap_eps')

clustered_hpS <- left_join(pSNP, over20_hhCounts) %>%
  mutate_if(is.character,function(x){replace_na(x, '0')}) %>%
  select(1,ncol(.))

return(clustered_hpS)
}
