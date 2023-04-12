#' Visualize haplotypes
#'
#' crosshap_viz() builds five individual plots using various elements of a
#' HapObject created by run_haplotyping(). The central dotplot displays
#' relationship between clusters of linked SNPs (Marker Groups), and distinct
#' haplotypes present within the population. Vertical plots (top/bottom)
#' visualize individuals and populations, grouped by haplotype. Horizontal plots
#' (left/right) visualize SNP information, grouped by Marker Group cluster.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param epsilon Epsilon to visualize haplotyping results for.
#' @param plot_left When plot_left = "allele", SNP allele frequency information
#' is displayed, when plot_left = "pos", SNP position information is displayed.
#' @param plot_right When plot_right = "pheno", phenotype associations for SNPs
#' are displayed, when plot_right = "cluster", internal marker group linkage is
#' displayed.
#' @param hide_labels When TRUE, legends from plots are hidden.
#' @param isolate_group If one or more Metadata groups are provided, all other
#' Metadata groups will be masked from the plot. NOTE: it does change the
#' summary tables or marker group phenotype scores.
#'
#' @export
#'
#' @return A patchwork object.
#'
crosshap_viz <- function(HapObject, epsilon, plot_left = "allele", plot_right = "pheno",
                         hide_labels = FALSE, isolate_group = NA) {

#Extract haplotype results for given epsilon
  for (x in 1:length(HapObject)){
    if(HapObject[[x]]$epsilon == epsilon){
      HapObject_eps <- HapObject[[x]]
    }
  }

  if(sum(is.na(HapObject_eps$Indfile$Pheno)) > 0){
    base::message(paste0(sum(is.na(HapObject_eps$Indfile$Pheno)), " individuals without phenotype information"))
  }

  mid <- build_mid_dotplot(HapObject, epsilon = epsilon, hide_labels = hide_labels)

  top <- build_top_metaplot(HapObject, epsilon = epsilon, hide_labels =  hide_labels)

  bot <- build_bot_halfeyeplot(HapObject, epsilon = epsilon, hide_labels =  T, isolate_group = isolate_group)

  left <- switch(plot_left,
                 "allele" = build_left_alleleplot(HapObject, epsilon = epsilon, hide_labels = hide_labels),
                 "pos" = build_left_posplot(HapObject, epsilon = epsilon, hide_labels = hide_labels))


  right <- switch(plot_right,
                  "pheno" = build_right_phenoplot(HapObject, epsilon = epsilon, hide_labels = hide_labels),
                  "cluster" = build_right_clusterplot(HapObject, epsilon = epsilon, hide_labels = hide_labels))

  tables <- build_summary_tables(HapObject, epsilon = epsilon)
  MGtable <- tables[[1]]
  haptable <- tables[[2]]

  layout <- "#BF
  DAE
  HCG"

  crosshap_stitched <-
    patchwork::wrap_plots(mid, top, bot , left, right, MGtable, haptable) +
    patchwork::guide_area() +
    patchwork::plot_layout(design = layout, guides = "collect")

  return(crosshap_stitched)
}

