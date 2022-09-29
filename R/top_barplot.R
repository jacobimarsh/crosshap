#' Top Metadata-hap Plot
#'
#' Internal function that creates a vertical stacked barplot representing the
#' frequency of each haplotype combination, broken down by each categorical
#' metadata variable. Makes use of the $Hapfile pseudoSNP count information. May
#' be missing some axis to allow for 'crosshap' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @export
#'
#'

build_top_metaplot <- function(HapObject, hide_labels) {


topplot_data <- suppressMessages(HapObject$Indfile %>%
  dplyr::group_by(hap, Metadata) %>%
  dplyr::summarise(counts = length(Metadata)) %>%
  dplyr::filter(hap != 0))


top_metaplot <- ggplot2::ggplot(topplot_data,
       ggplot2::aes(y = counts, x = hap, fill = base::factor(Metadata, levels = c(NA, sort(unique(HapObject$Indfile$Metadata))), exclude = NULL))) +
  #ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_bar(position="stack", stat = "identity") +
  ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer("Metadata", palette = "Dark2", na.value = "grey20") +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8),
                   legend.direction = "horizontal",
        legend.text = ggplot2::element_text(size = 7),
        legend.key.size = ggplot2::unit(6, "mm"),
        axis.text.x = ggplot2::element_text(face = "bold", size = 10, colour = "black"),
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        axis.text.y = ggplot2::element_text(face = "bold",
                                   size = 10, colour = "black"),
        axis.title.x = ggplot2::element_blank()) +
  ggplot2::scale_y_continuous(expand = c(0,0)) +
  ggplot2::ylab("Population size") +
  ggplot2::xlab("Haplotype combination")

if(hide_labels == T){
  return(top_metaplot + ggplot2::theme(legend.position = "none"))
} else {
  return(top_metaplot)
}
}
