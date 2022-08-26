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
#'
#' @return
#' @export
#'
#' @example
#' build_left_posplot(Haplotypes_MP_E2)
#'

build_left_posplot <- function(HapObject) {
IQRs <- c(((base::max(HapObject$MGfile$POS) - base::min(HapObject$MGfile$POS))*0.1 + base::min(HapObject$MGfile$POS)),
          ((base::max(HapObject$MGfile$POS) - base::min(HapObject$MGfile$POS))*0.5 + base::min(HapObject$MGfile$POS)),
          ((base::max(HapObject$MGfile$POS) - base::min(HapObject$MGfile$POS))*0.9 + base::min(HapObject$MGfile$POS)))


left_posplot <- dplyr::filter(HapObject$MGfile) %>% dplyr::filter(cluster > 0) %>% ggplot2::ggplot() +
  ggplot2::geom_segment(size = 0.2,
                        ggplot2::aes(x = POS, xend = POS, y = cluster-0.2, yend = cluster+0.2)) +
  ggplot2::scale_x_continuous(breaks = IQRs) +
  ggplot2::scale_y_reverse(breaks = 1:base::max(HapObject$MGfile$cluster),
    labels = c(paste0("MG", base::as.character(max(HapObject$MGfile$cluster):1))),
    position = "right") +
  ggplot2::labs(x = "Position") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(face = "bold", size = 10, color = "black"),
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        axis.title.x = ggplot2::element_text(face = "bold", size = 10, color = "black"))

return(left_posplot)
}
