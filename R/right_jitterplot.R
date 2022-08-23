#' Right SNP-Pheno jitterplot
#'
#' Internal function that creates a horizontal jitterplot displaying the
#' difference in phenotype means between alternate and reference alleles for
#' each SNP loci, grouped by marker group. Makes use of the $Varfile phenotypic
#' information for each allele, first calculating the difference between alt/ref
#' before plotting. May be missing some axis to allow for 'crosshap' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @return
#' @export
#'
#' @example
#' build_right_jitterplot(Haplotypes_MP_E2)
#'

build_right_jitterplot <- function(HapObject) {
  pdiff_altAF <- HapObject$Varfile %>%
  select(-nInd) %>%
  spread(key, avPheno) %>%
  rename(ref = '0', het = '1', alt = '2') %>%
  mutate(percdiff = alt - ref) %>%
  select(cluster, ID, percdiff) %>%
  left_join(HapObject$Varfile %>%
              select(-avPheno) %>%
              spread(key, nInd) %>%
              rename(ref = '0', het = '1', alt = '2') %>%
              mutate(AltAF = alt/(ref + het + alt)) %>%
              select(cluster, ID, AltAF), by = c("ID","cluster"))

#D right plot

right_jitterplot <- ggplot() +
  geom_jitter(data = pdiff_altAF %>% dplyr::filter(cluster > 0),
              aes(x = abs(percdiff), y = as.factor(cluster), fill = AltAF),
              alpha = 0.25, pch = 21, height = 0.25) +
  scale_fill_gradient('Alt allele frequency', low = 'white', high = '#440154FF') +
                        scale_x_continuous(breaks = scales::pretty_breaks()) +
                        theme_minimal() +
                        theme(axis.title.y = element_blank(),
                              axis.text.y = element_text(face = "bold", size = 10, color = "black"),
                              plot.margin = unit(c(0,0,0,0), "cm"),
                              legend.title = element_text(size = 7),
                              legend.text = element_text(size = 5),
                              legend.position = "none",
                              legend.key.size = unit(5, "mm"),
                              plot.title = element_blank(),
                              axis.text.x = element_text(face = "bold", size = 10),
                              axis.title.x = element_text()) +
  xlab("Pheno association") +
  scale_y_discrete(limits = rev, position = "left",
                   labels = c(paste0("MG",as.character(max(HapObject$Varfile$cluster):1))))

return(right_jitterplot)
}
