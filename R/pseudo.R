#' Call allelic states for each SNP marker group across individuals
#'
#' @param db40_clust SNP marker groups clustered using DBscan.
#' @param vcf Input VCF for region of interest.
#'
#' @return
#' @export
#'
#' @examples
create_pseudoSNP <- function(db40_clust, vcf) {

  db40_c1 <- db40_clust %>%
    dplyr::filter(cluster == 1) %>% tibble::as_tibble()

  c1_vcf <- vcf %>%
    rownames_to_column() %>% dplyr::filter(rowname %in% db40_c1$POS) %>% column_to_rownames()

  i_modes_1 <- apply(c1_vcf %>% sapply(as.double), 2, crosshap::mode) %>% as_tibble() %>%
    pull(value)

  pseudoSNP <- tibble(ID = colnames(c1_vcf), MG1 = i_modes_1)

for (vel in c(2:max(db40$cluster))) {

    db40_cvel <- db40_clust %>%
      dplyr::filter(cluster == vel) %>% as_tibble()

    cvel_vcf <- vcf %>%
      rownames_to_column() %>% dplyr::filter(rowname %in% db40_cvel$POS) %>% column_to_rownames()

    i_modes_vel <-
      apply(cvel_vcf %>% sapply(as.double), 2, mode) %>% as_tibble()

    pseudoSNP <- bind_cols(pseudoSNP, i_modes_vel)

    colnames(pseudoSNP)[vel + 1] <- paste0("MG", vel)
  }

  het_pseudoSNP <- pseudoSNP %>%
    mutate_if(is.numeric,function(x) {gsub(1, 2, x,fixed = T)})

  return(het_pseudoSNP)
}
