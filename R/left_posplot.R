#' Left SNP posplot
#'
#' Internal function that plots the position of each SNP in each marker group
#' horizontally using small vertical lines. Makes use of the $MGfile position
#' information for each allele. The first snippet of code determines the labels
#' to use, as genomic position can be exceptionally high numbers that end up
#' being cut-off when default labeling is used. May be missing some axis to
#' allow for 'crosshap' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#' @param hide_labels
#'
#' @return
#' @export
#'
#' @example
#' build_left_posplot(Haplotypes_MP_E2)
#'

build_left_posplot <- function(HapObject, hide_labels) {

IQRs <- base::as.numeric(HapObject$MGfile$POS) %>% {c(((base::max(.) - base::min(.))*0.1 + base::min(.)),
                                        ((base::max(.) - base::min(.))*0.5 + base::min(.)),
                                        ((base::max(.) - base::min(.))*0.9 + base::min(.)))}

left_posplot <- HapObject$MGfile %>% dplyr::filter(MGs != 0) %>% dplyr::mutate(MGs = as.numeric(stringr::str_remove(MGs,"MG"))) %>%
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
