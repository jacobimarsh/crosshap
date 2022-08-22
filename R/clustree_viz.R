#' Visualize clusters
#'
#' This function reads haplotype objects calculated for a range of epsilon
#' values. It calculates the mean phenotypic value for each haplotype group
#' and rearranges the data into a single sheet. Visualization is conducted using
#' the clustree package, with labeling manually added afterward.
#'
#' @param epsilon Epsilon values to search.
#' @param pheno Input numeric phenotype data for each individual.
#'
#' @return
#' @export
#'
#' @examples
#'
run_clustree <- function(epsilon, pheno) {

#Extract ID file first epsilon value and change column name to hap_epsXX
pre_clustree <- base::get(base::paste("Haplotypes_MP_E",epsilon[1],sep=""))[["Indfile"]] %>%
    dplyr::rename(!!base::paste0("hap_eps",epsilon[1]) := 'hap')

#Iterate over all other epsilon values, adding hap_epsXX columns to tibble
for (drez in epsilon[2:base::length(epsilon)]){
  pre_clustree <- pre_clustree %>%
    dplyr::left_join(base::get(base::paste("Haplotypes_MP_E",drez,sep=""))[["Indfile"]] %>%
                dplyr::rename(!!base::paste0("hap_eps",drez) := 'hap'))
}

pre_clustree_phen <- dplyr::left_join(pre_clustree, pheno)

#Plot with clustree
ctree <-
  clustree::clustree(pre_clustree_phen, prefix = 'hap_eps', node_colour = 'Pheno',node_colour_aggr = 'mean_na.rm', edge_width = 1, node_alpha = 1)+
  ggplot2::scale_colour_gradient(limits=c(base::max(dplyr::top_frac(pre_clustree_phen,-0.1,Pheno)$Pheno),
                                 base::min(dplyr::top_frac(pre_clustree_phen,0.1,Pheno)$Pheno)),
                        high = "#8ADD81",low = "#6870F6",oob = scales::squish,name = 'Pheno') +
  ggraph::scale_edge_color_continuous(high = 'black',low = 'grey80') +
  ggplot2::labs(size = 'nIndividuals', edge_alpha = "Proportion") +
  ggplot2::guides(edge_color = "none", size = ggplot2::guide_legend(order = 1))

#Create epsilon label data
#Extract x and ay coordinates from clustree object and build labels
lbls <-
  tibble::tibble(xval = base::max(ctree[["data"]][["x"]])*1.1,
         yval=0:(base::length(epsilon)-1),
         labelval = base::paste0("\u03b5"," = ",epsilon))

#Re-plot with epsilon label data added
labeled_ctree <- ctree +
  ggplot2::geom_text(data = lbls, ggplot2::aes(x=xval, y=yval, label=labelval), hjust = 0)

return(labeled_ctree)
}

