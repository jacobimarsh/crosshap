#' Build 'crosshap' visualization
#'
#' Builds five individual plots using various elements of a HapObject created by
#' crosshap::run_haplotyping. The central dotplot displays relationship between
#' clusters of linked SNPs (marker groups), and distinct haplotypes present
#' within the population. Vertical plots (top/bottom) visualize individuals and
#' populations, grouped by haplotype. Horizontal plots (left/right) visualize
#' SNP information, grouped by marker group cluster.
#'
#' @param plot_left
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @return Returns a crosshap visualization
#' @export
#'
#' @example
#' crosshap_viz(Haplotypes_MP_E2)
#'
crosshap_viz <- function(HapObject, plot_left = "allele") {
  base::message(paste0("Building Mid Dot plot"))
  mid <- build_mid_dotplot(HapObject)

  base::message(paste0("Building Top Metadata-Hap plot"))
  top <- build_top_metaplot(HapObject)

  base::message(paste0("Building Bottom Hap-Pheno plot"))
  bot <- build_bot_jitterplot(HapObject)

  base::message(paste0("Building Left SNP info plot"))
  # if (left_plot == "allele") {
  #   left <- build_left_alleleplot(HapObject)
  # } else if (left_plot == "pos") {
  #   left <- build_left_posplot(HapObject)
  # }

  left <- switch(plot_left,
                 "allele" = build_left_alleleplot(HapObject),
                 "pos" = build_left_posplot(HapObject))

  base::message(paste0("Building Right SNP-Pheno plot"))
  right <- build_right_jitterplot(HapObject)

  layout <- "#B#
  DAE
  #CF"

  base::message(paste0("Stitching plots"))
  crosshap_stitched <- mid + top + bot + left + right + patchwork::guide_area() + patchwork::plot_layout(design = layout)

  base::message(paste0("Done!"))
  return(crosshap_stitched)
}
