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
#' @param epsilon Epsilon to visualize haplotyping results for.
#' @param hide_labels If TRUE, legend is hidden.
#'
#' @importFrom rlang ".data"
#'
#' @export
#'
#' @return A ggplot2 object.
#'
#' @examples
#'     build_right_phenoplot(HapObject, epsilon = 0.6, hide_labels = FALSE)
#'

build_right_phenoplot <- function(HapObject, epsilon, hide_labels = TRUE) {

  #Extract haplotype results for given epsilon
  for (x in 1:length(HapObject)){
    if(HapObject[[x]]$epsilon == epsilon){
      HapObject_eps <- HapObject[[x]]
    }
  }

right_phenoplot <- ggplot2::ggplot() +
  ggplot2::geom_jitter(data = HapObject_eps$Varfile %>% dplyr::filter(.data$MGs != 0),
                       ggplot2::aes(x = base::abs(.data$phenodiff),
                                    y = base::factor(.data$MGs,
                                                     levels = base::paste0("MG",
                                                                           1:(base::length(base::unique(c(HapObject_eps$Varfile$MGs)))-1))),
                                    fill = .data$AltAF),
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
                   labels = c(paste0("MG", base::as.character((base::length(unique(HapObject_eps$Varfile$MGs))-1):1))))

if(hide_labels == T){
  return(right_phenoplot + ggplot2::theme(legend.position = "none"))
} else {
  return(right_phenoplot)
}
}
