#' Middle MG/Hap dotplot
#'
#' Internal function that creates central dotplot displaying the relationship
#' between haplotype combinations and the characteristic marker groups they
#' possess. Makes use of the $Hapfile 'pseudoSNP' information. Required as an
#' anchor for all peripheral plots in the 'crosshap' visualization.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @return
#' @export
#'
#' @example
#' build_mid_dotplot(Haplotypes_MP_E2)
#'

build_mid_dotplot <- function(HapObject) {
intersect <- HapObject$Hapfile %>%
  tidyr::gather(MG, present, 1:(base::ncol(.)-2)) %>%
  dplyr::mutate(present = base::as.factor(present)) %>%
  dplyr::mutate(MG = base::as.numeric(str_remove(MG, "MG")))

intersect_lines <- intersect %>%
  dplyr::filter(present == 2) %>%
  dplyr::group_by(hap) %>%
  dplyr::summarise(max = base::max(MG), min = base::min(MG)) %>%
  dplyr::mutate(min = base::as.character(min), max = base::as.character(max))

mid_dotplot <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = intersect_lines, col = "grey", size = 1.5,
               ggplot2::aes(x = hap, xend = hap, y = min, yend = max)) +
  ggplot2::geom_point(data = intersect, col = 'black', pch = 21,
             ggplot2::aes(hap, base::as.character(MG), fill = present, size= 2)) +
  ggplot2::scale_fill_manual(values = c('white', 'black')) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = 'none',
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        plot.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 10, face = 'bold', color = 'black'),
        axis.text.y = ggplot2::element_text(size = 10, face = 'bold', color = 'black')) +
  ggplot2::ylab("Marker Group") +
  ggplot2::xlab("Haplotype combination") +
  ggplot2::scale_y_discrete(limits = rev, position = "left",
                   labels = c(paste0("MG",base::as.character(base::max(intersect$MG):1))))

return(mid_dotplot)
}
