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
hap_run_cluster <- function(vcf, LD, elon =  c(0.5, 1, 1.5, 2)) {
  for (rez in elon) {
    db40 <- dbscan(LD, eps = rez, minPts = 40)

    db40_clust <- cbind(LD[, 0], db40$cluster) %>%
      rownames_to_column("POS") %>%
      rename("cluster" = "db40$cluster") %>% as_tibble()
    ##See if can export cluster names from dbscan

    assign(paste("clustree_export_clusts", rez,
                 sep = ""), db40_clust)

    write_tsv(db40_clust, paste0("SNP_MGs_MPx_E", rez, '.tsv'),
              quote = 'none')

    ##Step 1: Calculate  alternate vs reference counts for all SNPs in a given cluster

    db40_c1 <- db40_clust %>%
      filter(cluster == 1) %>%
      as_tibble()

    c1_vcf <- vcf %>%
      rownames_to_column() %>%
      filter(rowname %in% db40_c1$POS) %>%
      column_to_rownames()

    i_modes_1 <- apply(c1_vcf %>% sapply(as.double), 2, mode) %>%
      as_tibble() %>% pull(value)

    pseudoSNP <- tibble(ID = colnames(c1_vcf),
                        MG1 = i_modes_1)

    ##Step 2: Start creating pseudoSNP with all MGs
    for (vel in c(2:max(db40$cluster)))
    {
      db40_cvel <- db40_clust %>%
        filter(cluster == vel) %>%
        as_tibble()

      cvel_vcf <- vcf %>%
        rownames_to_column() %>%
        filter(rowname %in% db40_cvel$POS) %>%
        column_to_rownames()

      i_modes_vel <- apply(cvel_vcf %>%
                             sapply(as.double), 2, mode) %>%
        as_tibble()

      pseudoSNP <- bind_cols(pseudoSNP, i_modes_vel)

      colnames(pseudoSNP)[vel + 1] <- paste0("MG", vel)
    }
    ###################UP TO HERE @JAKOB###################
    #het to full hap
    het_pseudoSNP <- pseudoSNP %>% column_to_rownames("ID") %>%
      mutate_all(function(x) {
        gsub(1, 2, x,
             fixed = T)
      }) %>% rownames_to_column("ID") %>% as_tibble()

    write_tsv(het_pseudoSNP, paste0("ID_MGs_MPx_E", rez, '.tsv'),
              quote = 'none')

    ####NEED TO MAKE IT MUTATE ALL >BUT> ID

    cnames <- colnames(select(pseudoSNP,-ID))

    het_hapCounts <- het_pseudoSNP %>%
      gather(mgs, value, 2:ncol(.)) %>%
      group_by(ID) %>%
      mutate(mgs_new = paste(value, collapse = '_')) %>%
      group_by(mgs_new) %>%
      tally() %>% mutate(n = n / (ncol(pseudoSNP) - 1)) %>%
      separate(col = mgs_new,
               into = cnames,
               sep = "_") %>%
      as_tibble() %>%
      arrange(by_group = -n)

    over20_hhCounts <- filter(het_hapCounts, n > 9) %>%
      mutate(hap_eps = LETTERS[1:nrow(.)]) %>%
      rename(!!paste0("hap_eps", rez) := 'hap_eps')#%>%

    clustered_hpS <-
      inner_join(het_pseudoSNP, over20_hhCounts[, 1:ncol(over20_hhCounts)],
                 by = (colnames(over20_hhCounts[, 2:ncol(over20_hhCounts)]) =
                         c(colnames(het_pseudoSNP[, 2:ncol(het_pseudoSNP)]))))

    clustree_expor_hpS <-
      mutate(clustered_hpS[, c(1, ncol(clustered_hpS))])
    ##multiplex with a few different epsilon parameters
    rexport_hpS <-
      paste("Haplotype_assignments_MP", "_E", rez, ".txt",
            sep = "")

    assign(rexport_hpS, clustree_expor_hpS)

    write_tsv(clustree_expor_hpS,
              paste0("ID_Haps_MPx_E", rez, '.tsv'),
              quote = 'none')
  }

}
