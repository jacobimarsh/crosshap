#' Top metadata-hap bar plot
#'
#' build_top_metaplot() builds a vertical stacked bar plot displaying the
#' frequency of each haplotype combination, broken down by each categorical
#' metadata variable provided. Makes use of the $Indfile information from a
#' haplotype object. This is an in internal function called by crosshap_viz(),
#' though can be called separately to build a stand-alone plot
#'
#' @param HapObject Haplotype object created by run_haplotyping()
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
#' build_top_metaplot(HapObject, epsilon = 0.6, hide_labels = FALSE)
#'

build_top_metaplot <- function(HapObject, epsilon, hide_labels = FALSE) {

#Extract haplotype results for given epsilon
  for (x in 1:length(HapObject)){
    if(HapObject[[x]]$epsilon == epsilon){
      HapObject_eps <- HapObject[[x]]
    }
  }

#Remove unassigned individuals and convert Indfile to long format
topplot_data <- suppressMessages(HapObject_eps$Indfile %>%
  dplyr::group_by(.data$hap, .data$Metadata) %>%
  dplyr::summarise(counts = length(.data$Metadata)) %>%
  dplyr::filter(.data$hap != 0))


top_metaplot <- ggplot2::ggplot(topplot_data,
       ggplot2::aes(y = .data$counts, x = .data$hap, fill = base::factor(.data$Metadata, levels = c(NA, sort(unique(HapObject_eps$Indfile$Metadata))), exclude = NULL))) +
  ggplot2::geom_bar(position="stack", stat = "identity") +
  ggplot2::theme_minimal() +
    ggplot2::scale_fill_brewer("Metadata", palette = "Dark2", na.value = "grey20") +
    ggplot2::theme(legend.direction = "horizontal",
                   legend.title = ggplot2::element_text(size = 10),
                   legend.text = ggplot2::element_text(size = 7),
                   legend.key.size = ggplot2::unit(7, "mm"),
        axis.text.x = ggplot2::element_text(face = "bold", size = 10, colour = "black"),
        plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
        axis.text.y = ggplot2::element_text(face = "bold",
                                   size = 10, colour = "black"),
        axis.title.x = ggplot2::element_blank()) +
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5), title = "Metadata")) +
  ggplot2::scale_y_continuous(expand = c(0,0)) +
  ggplot2::ylab("Population size") +
  ggplot2::xlab("Haplotype combination")

if(hide_labels == T){
  return(top_metaplot + ggplot2::theme(legend.position = "none"))
} else {
  return(top_metaplot)
}
}
