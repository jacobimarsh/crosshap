#' Bot Hap-Pheno jitterplot
#'
#' Internal function that creates a vertical jitterplot displaying the
#' phenotypic scores for each individual, grouped by haplotype. Makes use of the
#' $Indfile information. May be missing some axis to allow for 'crosshap'
#' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#' @param hide_labels
#'
#' @return
#' @export
#'
#' @example
#' build_bot_jitterplot(Haplotypes_MP_E2)
#'

build_bot_jitterplot <- function(HapObject, hide_labels) {

no0data <- HapObject$Indfile %>% dplyr::filter(hap !=0)

bot_jitterplot <-
  ggplot2::ggplot(data = no0data) +
  ggplot2::geom_jitter(ggplot2::aes(x = hap, y = Pheno, fill = Metadata),
                       alpha = 0.3, pch = 21, width = 0.2)  +
  ggplot2::scale_fill_brewer(palette = "Dark2") +
  ggplot2::geom_crossbar(data = stats::aggregate(no0data$Pheno,
                                 base::list(no0data$hap), mean, na.rm = TRUE),
                         ggplot2::aes(x = as.factor(Group.1),
                    y = x,
                    xmin = base::as.factor(Group.1),
                    xmax = base::as.factor(Group.1),
                    ymin = x,
                    ymax = x,
                    colour = x),
                    show.legend = F) +
  ggplot2::scale_colour_gradient('Mean',
                        low='red',
                        high='green',
                        limits=c(base::max(dplyr::top_frac(no0data,
                                              -0.05,
                                              Pheno)$Pheno),
                                 base::min(dplyr::top_frac(no0data,
                                              0.05,
                                              Pheno)$Pheno)),
                        oob = scales::squish) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        axis.text.y = ggplot2::element_text(face = "bold", size = 10, color = "black")) +
  ggplot2::ylab("Pheno") +
  ggplot2::scale_y_continuous(position = "left", breaks = scales::pretty_breaks())

if(hide_labels == T){
  return(bot_jitterplot + ggplot2::theme(legend.position = "none"))
} else {
  return(bot_jitterplot)
}
}
