#' Left SNP-info allele plot
#'
#' build_left_alleleplot() builds a horizontal plot displaying mean allelic
#' frequencies (reference/alternate/missing/heterozygous) of all SNP loci,
#' grouped by marker group. Makes use of $Varfile information from a HapObject
#' created by run_haplotyping(). This is an internal function called by
#' crosshap_viz(), though can be called separately to build a stand-alone plot.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param hide_labels If TRUE, legend is hidden.
#'
#' @export
#'
#'

build_left_alleleplot <- function(HapObject, hide_labels = T) {
pre_leftplotdata <- HapObject$Varfile %>% dplyr::filter(MGs != 0)

#pre_leftplotdata[base::is.na(pre_leftplotdata)] <- 0

leftplot_data <- stats::aggregate(pre_leftplotdata$ref,
          base::list(pre_leftplotdata$MGs),
          mean) %>% dplyr::rename(MGs = 'Group.1', ref = 'x') %>%
  dplyr::left_join(stats::aggregate(pre_leftplotdata$miss,
                                    base::list(pre_leftplotdata$MGs),
                                    mean) %>% dplyr::rename(MGs = 'Group.1', miss = 'x'), by = "MGs") %>%
  dplyr::left_join(stats::aggregate(pre_leftplotdata$het,
                      base::list(pre_leftplotdata$MGs),
                      mean) %>% dplyr::rename(MGs = 'Group.1', het = 'x'), by = "MGs") %>%
  dplyr::left_join(stats::aggregate(pre_leftplotdata$alt,
                      base::list(pre_leftplotdata$MGs),
                      mean) %>% dplyr::rename(MGs = 'Group.1', alt = 'x'), by = "MGs") %>%
  dplyr::rename("REF" = "ref", "MISS" = "miss", "HET" = "het", "ALT" = "alt")

left_alleleplot <- ggplot2::ggplot(leftplot_data %>% tidyr::gather("Type", "nInd", 2:5) %>%
           dplyr::mutate(Type = base::factor(Type, levels = c("REF", "MISS", "HET", "ALT"))),
           ggplot2::aes(x = nInd,
             y = base::as.character(MGs),
             fill = Type,
             color = Type),
           ggplot2::ylim()) +
  ggplot2::geom_bar(ggplot2::aes(), position = 'stack', stat = "identity", width = 0.8, color = "black") +
  ggplot2::scale_x_reverse(breaks = scales::pretty_breaks(),
                  expand = c(0,0)) +
  ggplot2::theme_void() +
  ggplot2::theme(axis.text.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(),
        axis.title.y = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(size = 7),
        legend.text = ggplot2::element_text(size = 5),
        legend.key.size = ggplot2::unit(5,
                               "mm"),
        legend.direction = "horizontal",
        axis.text.x = ggplot2::element_text(face = "bold",
                                   size = 10),
        plot.margin = ggplot2::unit(c(0,0,0,0.1),
                           "cm"),
        plot.title = ggplot2::element_blank()) +
  ggplot2::scale_fill_manual(values = c("#FFFFFF", '#FDE725FF', "#73D055FF", '#440154FF')) +
  ggplot2::xlab("Allele count") +
  ggplot2::scale_y_discrete(position = "right", limits = rev,
                   labels = c(paste0("MG",base::as.character((base::length(unique(HapObject$Varfile$MGs))-1):1))))


if(hide_labels == T){
  return(left_alleleplot + ggplot2::theme(legend.position = "none"))
} else {
  return(left_alleleplot)
}
}
