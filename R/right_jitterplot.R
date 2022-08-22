#' Right SNP-Pheno jitterplot
#'
#' DESCRIPTION
#'
#' @param vcf Input VCF for region of interest.
#' @param LD Pairwise correlation matrix of SNPs in region from PLINK.
#' @param pheno Input numeric phenotype data for each individual.
#' @param epsilon Epsilon values for clustering SNPs with DBscan.
#' @param MGmin Minimum SNPs in marker groups, MinPts parameter for DBscan.
#'
#' @return
#' @export
#'
#' @example
#' run_haplotyping(vcf, LD, phen_early, epsilon, MGmin)
#'


MP_E1.5_pdiff_altAF <- Haplotypes_MP_E1.5$Varfile %>%
  select(-nInd) %>%
  spread(key, avPheno) %>%
  rename(ref = '0', het = '1', alt = '2') %>%
  mutate(percdiff = alt - ref) %>% select(ID, percdiff) %>%
  left_join(Haplotypes_MP_E1.5$Varfile %>%
              select(-avPheno) %>%
              spread(key, nInd) %>%
              rename(ref = '0', het = '1', alt = '2') %>%
              mutate(AltAF = alt/(ref + het + alt)) %>%
              select(ID, AltAF))


#D right plot

Dplot <- ggplot() +
  geom_jitter(data = MP_E1.5_pdiff_altAF %>% dplyr::filter(cluster > 0),
              aes(x = abs(percdiff), y = as.character(cluster), fill = AltAF),
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
                   labels = c(paste0("MG",as.character(max(MP_E1.69_pdiff_altAF$cluster):1))))


#top mid bot right
layout <-
"A#
 BD
 C#"

Aplot + Eplot + Bplot + Dplot + plot_layout(design = layout)
