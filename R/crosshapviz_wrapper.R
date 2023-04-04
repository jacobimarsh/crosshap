#' Visualize haplotypes
#'
#' crosshap_viz() builds five individual plots using various elements of a
#' HapObject created by run_haplotyping(). The central dotplot displays
#' relationship between clusters of linked SNPs (Marker Groups), and distinct
#' haplotypes present within the population. Vertical plots (top/bottom)
#' visualize individuals and populations, grouped by haplotype. Horizontal plots
#' (left/right) visualize SNP information, grouped by Marker Group cluster.
#'
#' @param plot_left When plot_left = "allele", SNP allele frequency information
#' is displayed, when plot_left = "pos", SNP position information is displayed.
#' @param plot_right When plot_right = "pheno", phenotype associations for SNPs
#' are displayed, when plot_right = "cluster", internal marker group linkage is
#' displayed.
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param hide_labels When TRUE, legends from plots are hidden.
#' @param isolate_group If one or more Metadata groups are provided, all other
#' Metadata groups will be masked from the plot. NOTE: it does change the
#' summary tables or marker group phenotype scores.
#'
#' @export
#'
#' @return A patchwork object.
#'
crosshap_viz <- function(HapObject, plot_left = "allele", plot_right = "pheno",
                         hide_labels = F, isolate_group = NA) {

  if(sum(is.na(HapObject$Indfile$Pheno)) > 0){
    base::message(paste0(sum(is.na(HapObject$Indfile$Pheno)), " individuals without phenotype information"))
  }

  mid <- build_mid_dotplot(HapObject, hide_labels)

  top <- build_top_metaplot(HapObject, hide_labels)

  bot <- build_bot_halfeyeplot(HapObject, hide_labels = T, isolate_group = isolate_group)

  left <- switch(plot_left,
                 "allele" = build_left_alleleplot(HapObject, hide_labels),
                 "pos" = build_left_posplot(HapObject, hide_labels))


  right <- switch(plot_right,
                  "pheno" = build_right_phenoplot(HapObject, hide_labels),
                  "cluster" = build_right_clusterplot(HapObject, hide_labels))

  tables <- build_summary_tables(HapObject)
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

