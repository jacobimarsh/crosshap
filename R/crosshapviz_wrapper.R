#' Visualize haplotypes
#'
#' crosshap_viz() builds five individual plots using various elements of a
#' HapObject created by run_haplotyping(). The central dotplot displays
#' relationship between clusters of linked SNPs (marker groups), and distinct
#' haplotypes present within the population. Vertical plots (top/bottom)
#' visualize individuals and populations, grouped by haplotype. Horizontal plots
#' (left/right) visualize SNP information, grouped by marker group cluster.
#'
#' @param plot_left When plot_left = "allele", SNP allele frequency information
#' is displayed, when plot_left = "pos", SNP position information is displayed.
#' @param plot_right When plot_right = "pheno", phenotype associations for SNPs
#' are displayed, when plot_right = "cluster", internal marker group linkage is
#' displayed.
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param hide_labels When TRUE, legends from plots are hidden.
#'
#' @export
#'
#' @return A patchwork object.
#'
#' @examples
#'
#' if (FALSE) {
#'      crosshap_viz(Haplotypes_MGmin30_E0.6)
#'}
#'
crosshap_viz <- function(HapObject, plot_left = "allele", plot_right = "pheno", hide_labels = F) {
 # base::message(paste0("Building Mid Dot plot"))
  mid <- build_mid_dotplot(HapObject, hide_labels)

 # base::message(paste0("Building Top Metadata-Hap plot"))
  top <- build_top_metaplot(HapObject, hide_labels)

 # base::message(paste0("Building Bottom Hap-Pheno plot"))
  bot <- build_bot_halfeyeplot(HapObject)

#  base::message(paste0("Building Left SNP info plot"))

  left <- switch(plot_left,
                 "allele" = build_left_alleleplot(HapObject, hide_labels),
                 "pos" = build_left_posplot(HapObject, hide_labels))

 # base::message(paste0("Building Right SNP-Pheno plot"))
  right <- switch(plot_right,
                  "pheno" = build_right_phenoplot(HapObject, hide_labels),
                  "cluster" = build_right_clusterplot(HapObject, hide_labels))

  tables <- build_summary_tables(HapObject)
  MGtable <- tables[[1]]
  haptable <- tables[[2]]

  layout <- "#BF
  DAE
  HCG"

  #base::message(paste0("Stitching plots"))
  crosshap_stitched <-
    patchwork::wrap_plots(mid, top, bot , left, right, MGtable, haptable) +
    patchwork::guide_area() +
    patchwork::plot_layout(design = layout, guides = "collect")

  if(sum(is.na(HapObject$Indfile$Pheno)) > 0){
  base::message(paste0(sum(is.na(HapObject$Indfile$Pheno)), " individuals without phenotypes omitted from botplot"))
}
  return(crosshap_stitched)
}

