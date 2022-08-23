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
top_metaplot <- ggplot(HapObject$Hapfile,
       aes(y = n, x = hap)) +
  geom_bar(position="stack", stat = "identity", fill = "black") +
  theme_minimal() +
  scale_fill_manual("Metadata") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(6, "mm"),
        axis.text.x = element_text(face = "bold", size = 10, colour = "black"),
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.text.y = element_text(face = "bold",
                                   size = 10, colour = "black"),
        axis.title.x = element_blank()) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Population size") +
  xlab("Haplotype combination")

return(top_metaplot)
}
