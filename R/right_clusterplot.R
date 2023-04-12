#' Right intra-cluster linkage plot
#'
#'
#' build_right_jitterplot() builds a horizontal plot displaying the mean
#' pairwise R^2 linkage between each SNP and all other SNPs in its marker group,
#' grouped by marker group, coloured by alternate allele frequency. Makes use of
#' the $Varfile information from haplotyping object. It is an internal function
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
#' build_right_clusterplot(HapObject, epsilon = 0.6, hide_labels = FALSE)
#'

build_right_clusterplot <- function(HapObject, epsilon, hide_labels = FALSE) {

  #Extract haplotype results for given epsilon
  for (x in 1:length(HapObject)){
    if(HapObject[[x]]$epsilon == epsilon){
      HapObject_eps <- HapObject[[x]]
    }
  }

  right_clusterplot <-   ggplot2::ggplot() +
    ggplot2::geom_boxplot(data = HapObject_eps$Varfile %>% dplyr::filter(.data$MGs != 0),
                         ggplot2::aes(x = .data$meanr2, y = .data$MGs),
                         alpha = 0.25, pch = 21, #height = 0.25,
                         coef = 2) +
  ggplot2::scale_fill_gradient('Alt allele frequency', low = 'white', high = '#440154FF')  +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()
                                ,limits = c(0.5,1)
                               ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "bold", size = 10, color = "black"),
                   plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = ggplot2::unit(7, "mm"),
                   legend.direction = "horizontal",
                   plot.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(face = "bold", size = 10),
                   axis.title.x = ggplot2::element_text()) +
    ggplot2::xlab("Mean intra-cluster R^2") +
    ggplot2::scale_y_discrete(limits = rev, position = "left",
                              labels = c(paste0("MG", base::as.character((base::length(unique(HapObject_eps$Varfile$MGs))-1):1))))

  if(hide_labels == T){
    return(right_clusterplot + ggplot2::theme(legend.position = "none"))
  } else {
    return(right_clusterplot)
  }
}
