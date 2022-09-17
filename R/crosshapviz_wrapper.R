#' Build 'crosshap' visualization
#'
#' Builds five individual plots using various elements of a HapObject created by
#' crosshap::run_haplotyping. The central dotplot displays relationship between
#' clusters of linked SNPs (marker groups), and distinct haplotypes present
#' within the population. Vertical plots (top/bottom) visualize individuals and
#' populations, grouped by haplotype. Horizontal plots (left/right) visualize
#' SNP information, grouped by marker group cluster. ABDCASCACSACS
#'
#' @param plot_left
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#' @param hide_labels
#'
#' @return Returns a crosshap visualization
#' @export
#'
#' @example
#' crosshap_viz(Haplotypes_MP_E2)
#'
crosshap_viz <- function(HapObject, plot_left = "allele", hide_labels = F) {
  base::message(paste0("Building Mid Dot plot"))
  mid <- build_mid_dotplot(HapObject)

  base::message(paste0("Building Top Metadata-Hap plot"))
  top <- build_top_metaplot(HapObject, hide_labels)

  base::message(paste0("Building Bottom Hap-Pheno plot"))
  bot <- build_bot_jitterplot(HapObject, hide_labels)

  base::message(paste0("Building Left SNP info plot"))

  left <- switch(plot_left,
                 "allele" = build_left_alleleplot(HapObject, hide_labels),
                 "pos" = build_left_posplot(HapObject, hide_labels))

  base::message(paste0("Building Right SNP-Pheno plot"))
  right <- build_right_jitterplot(HapObject, hide_labels)

  layout <- "#B#
  DAE
  #CF"

  base::message(paste0("Stitching plots"))
  crosshap_stitched <-
    patchwork::wrap_plots(mid, top, bot, left, right) +
    patchwork::guide_area() +
    patchwork::plot_layout(design = layout, guides = "collect")

  base::message(paste0("Done!"))
  return(crosshap_stitched)
}

