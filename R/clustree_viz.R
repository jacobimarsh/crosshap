#' Visualize clusters
#'
#' This function reads haplotype objects calculated for a range of epsilon
#' values. It calculates the mean phenotypic value for each haplotype group
#' and rearranges the data into a single sheet. Visualization is conducted using
#' the clustree package, with labeling manually added afterward.
#'
#' @param epsilon Epsilon values to search.
#' @param pheno Input numeric phenotype data for each individual.
#' @param type Hap or MG clustree
#'
#' @return
#' @export
#'
#' @examples
#'
run_clustree <- function(epsilon = c(0.4,0.8,1.2,1.6,2), MGmin, pheno, type = "MG") {
#Extract ID file first epsilon value and change column name to hap_epsXX
pre_clustree <- base::get(base::paste("Haplotypes_MGmin",MGmin,"_E",epsilon[1], sep=""))[["Indfile"]] %>%
    dplyr::rename(!!base::paste0("hap_eps",epsilon[1]) := 'hap')

#Iterate over all other epsilon values, adding hap_epsXX columns to tibble
for (drez in epsilon[2:base::length(epsilon)]){
  pre_clustree <- pre_clustree %>%
    dplyr::left_join(base::get(base::paste("Haplotypes_MGmin",MGmin,"_E",drez,sep=""))[["Indfile"]] %>%
                dplyr::rename(!!base::paste0("hap_eps",drez) := 'hap'), by = c("Ind", "Pheno"))
}

pre_clustree_phen <- dplyr::left_join(pre_clustree, pheno, by = c("Ind", "Pheno"))

#Plot with clustree
haptree <- base::suppressMessages(
  clustree::clustree(pre_clustree_phen, prefix = 'hap_eps', node_colour = 'Pheno',node_colour_aggr = "mean_na.rm", edge_width = 1, node_alpha = 1)+
  ggplot2::scale_colour_gradient(limits=c(base::max(dplyr::top_frac(pre_clustree_phen,-0.1,Pheno)$Pheno),
                                 base::min(dplyr::top_frac(pre_clustree_phen,0.1,Pheno)$Pheno)),
                        high = "#8ADD81",low = "#6870F6",oob = scales::squish,name = 'Pheno') +
  ggraph::scale_edge_color_continuous(high = 'black',low = 'grey80') +
  ggplot2::labs(size = 'nIndividuals', edge_alpha = "Proportion") +
  ggplot2::guides(edge_color = "none", size = ggplot2::guide_legend(order = 1))
)
#Create epsilon label data
#Extract x and ay coordinates from clustree object and build labels
haplbls <-
  tibble::tibble(xval = base::max(haptree[["data"]][["x"]])*1.1,
         yval=0:(base::length(epsilon)-1),
         labelval = base::paste0("\u03b5"," = ",epsilon))

#Re-plot with epsilon label data added
labeled_haptree <- haptree +
  ggplot2::geom_text(data = haplbls, ggplot2::aes(x=xval, y=yval, label=labelval), hjust = 0)

pre_MGtree <- base::get(base::paste("Haplotypes_MGmin",MGmin,"_E",epsilon[1], sep=""))[["Varfile"]] %>%
  dplyr::rename(!!base::paste0("MGs_eps",epsilon[1]) := 'MGs')

for (drez in epsilon[2:base::length(epsilon)]){
  pre_MGtree <- dplyr::left_join(pre_MGtree,
                                 base::get(base::paste("Haplotypes_MGmin",MGmin,"_E",drez,sep=""))[["Varfile"]] %>%
                                   dplyr::select(ID, MGs)%>%
                                   dplyr::rename(!!base::paste0("MGs_eps",drez) := 'MGs'), by = "ID")
}

MGtree <- base::suppressMessages(
clustree::clustree(pre_MGtree, prefix = 'MGs_eps', node_colour = 'percdiff',node_colour_aggr = "mean_na.rm", edge_width = 1, node_alpha = 1) +
  ggplot2::scale_colour_gradient(limits=c(base::max(pre_MGtree$percdiff),
                                          base::min(pre_MGtree$percdiff)),
                                 high = "#8ADD81",low = "#6870F6",oob = scales::squish,name = 'percdiff') +
  ggraph::scale_edge_color_continuous(high = 'black',low = 'grey80') +
  ggplot2::labs(size = 'nSNPs', edge_alpha = "Proportion") +
  ggplot2::guides(edge_color = "none", size = ggplot2::guide_legend(order = 1))
)
#Create epsilon label data
#Extract x and ay coordinates from clustree object and build labels
MGlbls <-
  tibble::tibble(xval = base::max(MGtree[["data"]][["x"]])*1.1,
                 yval=0:(base::length(epsilon)-1),
                 labelval = base::paste0("\u03b5"," = ",epsilon))

#Re-plot with epsilon label data added
labeled_MGtree <- MGtree +
  ggplot2::geom_text(data = MGlbls, ggplot2::aes(x=xval, y=yval, label=labelval), hjust = 0)

return(switch(type,
       "hap" = labeled_haptree,
       "MG" = labeled_MGtree))
}

