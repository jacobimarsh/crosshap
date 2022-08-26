#' Top Metadata-hap Plot
#'
#' Internal function that creates a vertical stacked barplot representing the
#' frequency of each haplotype combination, broken down by each categorical
#' metadata variable. Makes use of the $Hapfile pseudoSNP count information. May
#' be missing some axis to allow for 'crosshap' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @return
#' @export
#'
#' @example
#' build_top_metaplot(Haplotypes_MP_E2)
#'

build_top_metaplot <- function(HapObject) {
top_metaplot <- ggplot2::ggplot(HapObject$Hapfile,
       ggplot2::aes(y = n, x = hap)) +
  ggplot2::geom_bar(position="stack", stat = "identity", fill = "black") +
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual("Metadata") +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 8),
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

return(top_metaplot)
}
