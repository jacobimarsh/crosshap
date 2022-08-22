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

intersect <- Haplotypes_MP_E1.5$Hapfile %>%
  gather(MG, present, 1:(ncol(.)-2)) %>%
  mutate(present = as.factor(present)) %>%
  mutate(MG = as.numeric(str_remove(MG, "MG")))

intersect_lines <- intersect %>%
  dplyr::filter(present == 2) %>%
  group_by(hap) %>%
  summarise(max = max(MG), min = min(MG)) %>%
  mutate(min = as.character(min), max = as.character(max))

Eplot <- ggplot() +
  geom_segment(data = intersect_lines, col = "grey", size = 1.5,
               aes(x = hap, xend = hap, y = min, yend = max)) +
  geom_point(data = intersect, col = 'black', pch = 21,
             aes(hap, as.character(MG), fill = present, size= 2)) +
  scale_fill_manual(values = c('white', 'black')) +
  theme_minimal() +
  theme(legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), "cm"),
        plot.title = element_blank(),
        axis.text.x = element_text(size = 10, face = 'bold', color = 'black'),
        axis.text.y = element_text(size = 10, face = 'bold', color = 'black')) +
  ylab("Marker Group") +
  xlab("Haplotype combination") +
  scale_y_discrete(limits = rev, position = "left",
                   labels = c(paste0("MG",as.character(max(intersect$MG):1))))
