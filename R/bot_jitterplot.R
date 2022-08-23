#' Bot Hap-Pheno jitterplot
#'
#' Internal function that creates a vertical jitterplot displaying the
#' phenotypic scores for each individual, grouped by haplotype. Makes use of the
#' $Indfile information. May be missing some axis to allow for 'crosshap'
#' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @return
#' @export
#'
#' @example
#' build_bot_jitterplot(Haplotypes_MP_E2)
#'

build_bot_jitterplot <- function(HapObject) {
bot_jitterplot <- ggplot(data = HapObject$Indfile) +
  geom_jitter(aes(x = hap, y = Pheno), alpha = 0.25, pch = 21, width = 0.2) +
  geom_crossbar(data = aggregate(HapObject$Indfile$Pheno,
                                 list(HapObject$Indfile$hap), mean, na.rm = TRUE),
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
                        limits=c(max(top_frac(HapObject$Indfile,
                                              -0.05,
                                              Pheno)$Pheno),
                                 min(top_frac(HapObject$Indfile,
                                              0.05,
                                              Pheno)$Pheno)),
                        oob = scales::squish) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_text(face = "bold", size = 10, color = "black")) +
  ylab("Pheno") +
  scale_y_continuous(position = "left", breaks = scales::pretty_breaks())

return(bot_jitterplot)
}
