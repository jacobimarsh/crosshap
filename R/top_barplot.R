#' Middle MG/Hap dotplot
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


#A top plot
#Haplotypes_MP_E1.69$IDfile %>% group_by(hap) %>% tally()

Aplot <- ggplot(Haplotypes_MP_E1.5$Hapfile,
       aes(y = n, x = hap)) +
  geom_bar(position="stack", stat = "identity", fill = "black") +
  theme_minimal() +
  scale_fill_manual("Metadata") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(6, "mm"),
        axis.text.x = element_text(face = "bold", size = 10),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_text(face = "bold",
                                   size = 10),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Population size") +
  xlab("Haplotype combination")


#top and mid
Aplot + Eplot + plot_layout(ncol = 1)
