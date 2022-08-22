#' Bot Hap-Pheno jitterplot
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


#B bot plot

Bplot <- ggplot(data = Haplotypes_MP_E1.5$Indfile) +
  geom_jitter(aes(x = hap, y = Pheno), alpha = 0.25, pch = 21, width = 0.2) +
  geom_crossbar(data = aggregate(Haplotypes_MP_E1.5$Indfile$Pheno,
                                 list(Haplotypes_MP_E1.5$Indfile$hap), mean, na.rm = TRUE),
                aes(x = as.factor(Group.1),
                    y = x,
                    xmin= as.factor(Group.1),
                    xmax=as.factor(Group.1),
                    ymin=x,
                    ymax=x,
                    colour=x)) +
  scale_colour_gradient('Mean',
                        low='red',
                        high='green',
                        limits=c(max(top_frac(Haplotypes_MP_E1.5$Indfile,
                                              -0.05,
                                              Pheno)$Pheno),
                                 min(top_frac(Haplotypes_MP_E1.5$Indfile,
                                              0.05,
                                              Pheno)$Pheno)),
                        oob = squish) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black")) +
  ylab("Pheno") +
  scale_y_continuous(position = "left", breaks = scales::pretty_breaks())


#top mid bot
Aplot + Eplot + Bplot + plot_layout(ncol = 1)
