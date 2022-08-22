#' Left SNP-info alleleplot
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

pre_leftplotdata <- Haplotypes_MP_E1.5$Varfile %>% dplyr::filter(cluster > 0) %>%
  select(-avPheno) %>%
  spread(key, nInd) %>%
  rename(ref = '0', het = '1', alt = '2')

leftplot_data <- aggregate(pre_leftplotdata$ref,
          list(pre_leftplotdata$cluster),
          mean) %>% rename(cluster = 'Group.1', ref = 'x') %>%
  left_join(aggregate(pre_leftplotdata$het,
                      list(pre_leftplotdata$cluster),
                      mean) %>% rename(cluster = 'Group.1', het = 'x')) %>%
  left_join(aggregate(pre_leftplotdata$alt,
                      list(pre_leftplotdata$cluster),
                      mean) %>% rename(cluster = 'Group.1', alt = 'x'))

Cplot <- ggplot(leftplot_data %>% gather("Type", "nInd", 2:4) %>%
           mutate(Type = factor(Type, levels = c("ref", "het", "alt"))),
         aes(x = nInd,
             y = as.character(cluster),
             fill = Type,
             color = Type),
         ylim()) +
  geom_bar(aes(), position = 'stack', stat = "identity", width = 0.8, color = "black") +
  scale_x_reverse(breaks = scales::pretty_breaks(),
                  expand = c(0,0)) +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.key.size = unit(5,
                               "mm"),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.1),
                           "cm"),
        plot.title = element_blank()) +
  scale_fill_manual(values = c("#FFFFFF", "#73D055FF", '#440154FF')) +
  xlab("Allele count") +
  scale_y_discrete(position = "right", limits = rev,
                   labels = c(paste0("MG",as.character(max(Haplotypes_MP_E1.5$Varfile$cluster):1))))


#top mid bot right left
layout <-
  "#A#
 EBD
 #C#"

Aplot + Eplot + Bplot + Dplot + Cplot + plot_layout(design = layout)
