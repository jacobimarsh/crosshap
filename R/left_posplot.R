#' Left SNP posplot
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

IQRs <- c(((max(Haplotypes_MP_E1.5$MGfile$POS) - min(Haplotypes_MP_E1.5$MGfile$POS))*0.1 + min(Haplotypes_MP_E1.5$MGfile$POS)),
          ((max(Haplotypes_MP_E1.5$MGfile$POS) - min(Haplotypes_MP_E1.5$MGfile$POS))*0.5 + min(Haplotypes_MP_E1.5$MGfile$POS)),
          ((max(Haplotypes_MP_E1.5$MGfile$POS) - min(Haplotypes_MP_E1.5$MGfile$POS))*0.9 + min(Haplotypes_MP_E1.5$MGfile$POS)))


Fplot <- filter(Haplotypes_MP_E1.5$MGfile) %>% dplyr::filter(cluster > 0) %>% ggplot() +
  geom_segment(size = 0.2,
    aes(x = POS, xend = POS, y = cluster-0.2, yend = cluster+0.2)) +
  scale_x_continuous(breaks = IQRs) +
  scale_y_reverse(breaks = 1:max(Haplotypes_MP_E1.5$MGfile$cluster),
    labels = c(paste0("MG",as.character(max(Haplotypes_MP_E1.5$MGfile$cluster):1))),
    position = "right") +
  labs(x = "Position") +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 10, color = "black"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.x = element_text(face = "bold", size = 10, color = "black"))

#top mid bot right left
layout <-
  "#A#
 EBD
 #C#"

Aplot + Eplot + Bplot + Dplot + Fplot + plot_layout(design = layout)

