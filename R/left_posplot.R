#' Left SNP-position plot
#'
#' build_left_alleleplot() builds a horizontal plot displaying the chromosomal
#' position of each SNP locus, grouped by marker group. Makes use of the $Varfile
#' file from haplotype object. It is an internal function called by
#' crosshap_viz(), though can be called separately to build a stand-alone plot.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param hide_labels If TRUE, legend is hidden.
#'
#' @export
#'
#' @return A ggplot2 object.
#'
#' @examples
#' build_left_posplot(Haplotypes_MGmin30_E0.6, hide_labels = FALSE)
#'

build_left_posplot <- function(HapObject, hide_labels = T) {

#Find chromosomal positions at 10%, 50% and 90% intervals across region of interest
#These act as the x-axis labels
IQRs <- base::as.numeric(HapObject$Varfile$POS) %>% {c(((base::max(.) - base::min(.))*0.1 + base::min(.)),
                                        ((base::max(.) - base::min(.))*0.5 + base::min(.)),
                                        ((base::max(.) - base::min(.))*0.9 + base::min(.)))} %>%
  plyr::round_any(1000)

left_posplot <- HapObject$Varfile %>% dplyr::filter(MGs != 0) %>% dplyr::mutate(MGs = as.numeric(stringr::str_remove(MGs,"MG"))) %>%
  ggplot2::ggplot() +
  ggplot2::geom_segment(size = 0.2,
                        ggplot2::aes(x = POS, xend = POS, y = MGs-0.2, yend = MGs+0.2)) +
  ggplot2::scale_x_continuous(breaks = IQRs) +
  ggplot2::scale_y_reverse(breaks = (base::length(unique(HapObject$Varfile$MGs))-1):1,
    labels = c(paste0("MG", base::as.character((base::length(unique(HapObject$Varfile$MGs))-1):1))),
    position = "right") +
  ggplot2::labs(x = "Position") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                 legend.title = ggplot2::element_text(size = 7),
                 legend.text = ggplot2::element_text(size = 5),
                 legend.key.size = ggplot2::unit(5, "mm"),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold", size = 10, color = "black"),
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10, color = "black"))

if(hide_labels == T){
  return(left_posplot + ggplot2::theme(legend.position = "none"))
} else {
  return(left_posplot)
}
}
