#' Middle MG/hap dot plot
#'
#' build_mid_dotplot() builds a central dot plot displaying the relationship
#' between haplotype combinations and the characteristic marker group alleles
#' that define them. Makes use of the $Hapfile information from a haplotype
#' object. This is an internal function called by crosshap_viz(), though can
#' be called separately to build a stand-alone plot (can be useful when patched
#' to a peripheral plot).
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#' @param hide_labels If TRUE, legend is hidden.
#'
#' @importFrom rlang ".data"
#'
#' @export
#'
#' @return A ggplot2 object.
#'
#' @examples
#' build_mid_dotplot(Haplotypes_MGmin30_E0.6, hide_labels = FALSE)
#'

build_mid_dotplot <- function(HapObject, hide_labels) {

#Recode hapfile to long format, with 0 as REF, 1 as HET and 2 as ALT (dots in plot)
intersect <- HapObject$Hapfile %>% tibble::as_tibble() %>%
  tidyr::gather("MG", "present", 3:(base::ncol(HapObject$Hapfile))) %>%
  dplyr::mutate(present = base::as.factor(.data$present)) %>%
  dplyr::mutate(MG = base::as.numeric(stringr::str_remove(.data$MG, "MG"))) %>%
  dplyr::mutate(present = gsub(as.factor(2),"ALT",
                               gsub(as.factor(1),"HET",
                                    gsub(as.factor(0),"REF",.data$present)))) %>%
  dplyr::mutate(Allele = factor(.data$present, levels = c("REF", "HET", "ALT")))

#Report min and max MG that each hap has an ALT allele for (edges in plot)
intersect_lines <- intersect %>%
  dplyr::filter(.data$Allele == "ALT") %>%
  dplyr::group_by(.data$hap) %>%
  dplyr::summarise(max = base::max(.data$MG), min = base::min(.data$MG)) %>%
  dplyr::mutate(min = base::as.character(min), max = base::as.character(max))

mid_dotplot <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = intersect_lines, col = "grey", linewidth = 1.5,
               ggplot2::aes(x = .data$hap, xend = .data$hap, y = min, yend = max)) +
  ggplot2::geom_point(data = intersect, col = 'black', pch = 21,
             ggplot2::aes(.data$hap, base::as.character(.data$MG), fill = .data$Allele, size= 2)) +
  ggplot2::scale_fill_manual(values = c('white','black','grey')) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.direction = 'horizontal',
                 legend.justification = "left",
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        plot.title = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 10, face = 'bold', color = 'black'),
        axis.text.y = ggplot2::element_text(size = 10, face = 'bold', color = 'black'),
        legend.title = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 7),
        legend.key.size = ggplot2::unit(7, "mm"),
 #       axis.title = ggplot2::element_blank()
 ) +
  ggplot2::ylab("Marker Group") +
  ggplot2::xlab("Haplotype combination") +
  ggplot2::guides('size' = "none",
                  fill = ggplot2::guide_legend(override.aes = list(size = 5), title = "Allele")) +
  ggplot2::scale_y_discrete(limits = rev, position = "left",
                   labels = c(paste0("MG",base::as.character(base::max(intersect$MG):1))))

if(hide_labels == T){
  return(mid_dotplot + ggplot2::theme(legend.position = "none"))
} else {
  return(mid_dotplot)
}
}
