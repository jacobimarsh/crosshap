#' Create pseudo snips!
#'
#' @param db40
#' @param vcf
#'
#' @return
#' @export
#'
#' @examples
create_pseudo_SNP <- function(db40, vcf) {
  db40_c1 <- db40 %>%
    dplyr::filter(cluster == 1) %>%
    tibble::as_tibble()

  c1_vcf <- vcf %>%
    rownames_to_column() %>%
    filter(rowname %in% db40_c1$POS) %>%
    column_to_rownames()

  i_modes_1 <- apply(c1_vcf %>% sapply(as.double), 2, crosshap::mode) %>%
    as_tibble() %>% pull(value)

  pseudoSNP <- tibble(ID = colnames(c1_vcf),
                      MG1 = i_modes_1)

  pseudoSNP
}
