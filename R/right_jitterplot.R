#' Right SNP-Pheno jitterplot
#'
#' Internal function that creates a horizontal jitterplot displaying the
#' difference in phenotype means between alternate and reference alleles for
#' each SNP loci, grouped by marker group. Makes use of the $Varfile phenotypic
#' information for each allele, first calculating the difference between alt/ref
#' before plotting. May be missing some axis to allow for 'crosshap' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#' @param hide_labels
#'
#' @return
#' @export
#'
#' @example
#' build_right_jitterplot(Haplotypes_MP_E2)
#'

build_right_jitterplot <- function(HapObject, hide_labels) {

#D right plot

right_jitterplot <- ggplot2::ggplot() +
  ggplot2::geom_jitter(data = HapObject$Varfile %>% dplyr::filter(MGs != 0),
                       ggplot2::aes(x = base::abs(percdiff), y = base::as.factor(MGs), fill = AltAF),
              alpha = 0.25, pch = 21, height = 0.25) +
  ggplot2::scale_fill_gradient('Alt allele frequency', low = 'white', high = '#440154FF') +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                              axis.text.y = ggplot2::element_text(face = "bold", size = 10, color = "black"),
                              plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
                              legend.title = ggplot2::element_text(size = 7),
                              legend.text = ggplot2::element_text(size = 5),
                              legend.key.size = ggplot2::unit(5, "mm"),
                              legend.direction = "horizontal",
                              plot.title = ggplot2::element_blank(),
                              axis.text.x = ggplot2::element_text(face = "bold", size = 10),
                              axis.title.x = ggplot2::element_text()) +
  ggplot2::xlab("Pheno association") +
  ggplot2::scale_y_discrete(limits = rev, position = "left",
                   labels = c(paste0("MG", base::as.character((base::length(unique(HapObject$Varfile$MGs))-1):1))))

if(hide_labels == T){
  return(right_jitterplot + ggplot2::theme(legend.position = "none"))
} else {
  return(right_jitterplot)
}
}
