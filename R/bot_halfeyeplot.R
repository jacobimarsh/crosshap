#' Bot hap-pheno raincloud plot
#'
#' build_bot_halfeyeplot() builds a vertical plot displaying the
#' phenotypic scores for each individual, grouped by haplotype, coloured by
#' metadata variable. Makes use of the $Indfile information from haplotype
#' object. It is an internal function called by crosshap_viz(), though can be
#' called separately to build a stand-alone plot.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param isolate_groups Only show specific metadata groups on bot plot.
#' @param hide_labels If TRUE, legend is hidden.
#'
#' @export
#'
#' @return A ggplot2 object.
#'
#' @examples
#'
#' if (FALSE) {
#'      build_bot_halfeyeplot(Haplotypes_MGmin30_E0.6, hide_labels = F)
#'}
#'

build_bot_halfeyeplot <- function(HapObject, hide_labels = T, isolate_groups = NULL) {

no0data <- tidyr::drop_na(if(is.null(isolate_groups)){
  HapObject$Indfile
} else {
  HapObject$Indfile %>% dplyr::filter(Metadata %in% isolate_groups)
}, Pheno) %>% dplyr::filter(hap !=0 )

bot_halfeyeplot <-
  ggplot2::ggplot(data = no0data, ggplot2::aes(x = hap, y = Pheno#, fill = Metadata
  )) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) +
  ggplot2::geom_jitter(alpha = 0.3, pch = 21, width = 0.1,
                       ggplot2::aes(fill = Metadata)) +
  ggplot2::scale_fill_brewer(palette = "Dark2") +
  ggplot2::theme_minimal() +
  ggplot2::ylab("Pheno")+
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks()) +
  ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                 axis.text.x  = ggplot2::element_blank(),
                 legend.key.size = ggplot2::unit(7,
                                                 "mm"),
                 legend.title = ggplot2::element_text(size = 10),
                 legend.text = ggplot2::element_text(size = 7),
                 plot.margin = ggplot2::unit(c(0,0,0,0), "cm"),
                 axis.text.y = ggplot2::element_text(face = "bold", size = 10, color = "black")) +
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5, alpha = 0.6), title = "Metadata"))

if(hide_labels == T){
  return(bot_halfeyeplot + ggplot2::theme(legend.position = "none"))
} else {
  return(bot_halfeyeplot)
}
}
