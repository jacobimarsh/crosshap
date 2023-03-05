#' Right SNP-pheno phenoplot
#'
#' build_right_phenoplot() builds a horizontal plot displaying the mean
#' difference in phenotype score between individuals with the alternate vs
#' reference alleles for each SNP locus, grouped by marker group, coloured by
#' the alternate allele frequency of each SNP. Makes use of the $Varfile
#' phenotypic information from haplotyping object. It is an internal function
#' called by crosshap_viz(), though can be called separately to build a
#' stand-alone plot.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param hide_labels If TRUE, legend is hidden.
#'
#' @export
#'
#' @return A ggplot2 object.
#'
#' @examples
#'
#' if (FALSE) {
#'     build_right_phenoplot(Haplotypes_MGmin30_E0.6, hide_labels = F)
#' }
#'

build_right_phenoplot <- function(HapObject, hide_labels) {

right_phenoplot <- ggplot2::ggplot() +
  ggplot2::geom_jitter(data = HapObject$Varfile %>% dplyr::filter(MGs != 0),
                       ggplot2::aes(x = base::abs(phenodiff), y = base::as.factor(MGs), fill = AltAF),
              alpha = 0.25, pch = 21, height = 0.25) +
  ggplot2::scale_fill_gradient('Alt allele frequency', low = 'white', high = '#440154FF') +
  ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                 axis.text.y = ggplot2::element_text(face = "bold", size = 10, color = "black"),
                 plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
                 legend.title = ggplot2::element_text(size = 10),
                 legend.text = ggplot2::element_text(size = 7),
                 legend.key.size = ggplot2::unit(5, "mm"),
                 legend.direction = "horizontal",
                 plot.title = ggplot2::element_blank(),
                 axis.text.x = ggplot2::element_text(face = "bold", size = 10),
                 axis.title.x = ggplot2::element_text()) +
  ggplot2::xlab("Pheno association") +
  ggplot2::scale_y_discrete(limits = rev, position = "left",
                   labels = c(paste0("MG", base::as.character((base::length(unique(HapObject$Varfile$MGs))-1):1))))

if(hide_labels == T){
  return(right_phenoplot + ggplot2::theme(legend.position = "none"))
} else {
  return(right_phenoplot)
}
}
