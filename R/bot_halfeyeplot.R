#' Bot hap-pheno raincloud plot
#'
#' build_bot_halfeyeplot() builds a vertical plot displaying the
#' phenotypic scores for each individual, grouped by haplotype, coloured by
#' metadata variable. Makes use of the $Indfile information from haplotype
#' object. It is an internal function called by crosshap_viz(), though can be
#' called separately to build a stand-alone plot.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param hide_labels If TRUE, legend is hidden.
#' @param isolate_group If a Metadata group is provided, all other Metadata
#' groups will be masked from the plot. NOTE: it does change the summary tables
#' or marker group phenotype scores.
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

build_bot_halfeyeplot <- function(HapObject, hide_labels = T, isolate_group = "none") {

no0data <- tidyr::drop_na(HapObject$Indfile, Pheno) %>% dplyr::filter(hap !=0 )

halfeyedat <-  if(is.na(isolate_group)) {
  no0data } else {
  tidyr::drop_na(HapObject$Indfile %>% dplyr::mutate(Pheno = ifelse(Metadata != isolate_group, NA, Pheno)),
                                 Pheno) %>% dplyr::filter(hap !=0 )
}

jitterdat <- if(is.na(isolate_group)) {
  no0data %>% dplyr::mutate(keep = "A")
} else {
  no0data %>% dplyr::mutate(keep = ifelse(Metadata != isolate_group, NA, "A"))
}

if(length(setdiff(HapObject$Indfile$hap,c(unique(halfeyedat$hap),0))) > 0){
jitterdat <- jitterdat %>% tibble::add_row(Ind = "Dummy",
                hap = setdiff(HapObject$Indfile$hap,c(unique(halfeyedat$hap),0)),
                Pheno = crosshap::mean_na.rm(HapObject$Indfile$Pheno),
                Metadata = HapObject$Indfile$Metadata[[1]],
                keep = NA)
}

bot_halfeyeplot <-
  ggplot2::ggplot(data = no0data, ggplot2::aes(x = hap, y = Pheno)) +
  ggdist::stat_halfeye(data = halfeyedat, adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) +
  ggplot2::geom_jitter(data = jitterdat,
                         pch = 21, width = 0.1,
                       ggplot2::aes(fill = Metadata, alpha = keep)) +
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
  ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5, alpha = 0.6),
                                               title = "Metadata"),
                  alpha = "none") +
  ggplot2::scale_alpha_manual(values = 0.3, na.value = 0)


if(!is.na(isolate_group)) {
phen <- halfeyedat %>%
  dplyr::filter(hap != 0) %>%
  dplyr::group_by(hap) %>%
  dplyr::summarise(phenav = mean_na.rm(Pheno)) %>%
  tidyr::spread(hap, phenav) %>%
  tibble::as_tibble() %>%
  dplyr::mutate_if(is.double, function(x){signif(x, 3)}) %>%
  dplyr::mutate_if(is.double, as.character) %>%
  dplyr::mutate(rname = "Pheno") %>%
  tibble::column_to_rownames("rname")
message(paste("Haplotype phenotype averages of", isolate_group ,"individuals only:"))
message(paste(capture.output({phen}), collapse = "\n"))
}

if(hide_labels == T){
  return(bot_halfeyeplot + ggplot2::theme(legend.position = "none"))
} else {
  return(bot_halfeyeplot)
}
}
