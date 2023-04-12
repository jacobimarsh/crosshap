#' Clustering tree
#'
#' clustree_viz() builds a clustering tree displaying changes in haplotype
#' assignment between individuals or changes in Marker Group assignment for
#' SNPs, across different epsilon values. This function is a `clustree` wrapper.
#'
#' @param HapObject A haplotyping object with a range of results from different
#' epsilons created by run_haplotyping()
#' @param type When type = "hap", nodes represent haplotype populations, when
#' type = "MG", nodes represent marker groups.
#'
#' @importFrom rlang ".data"
#' @importFrom rlang ":="
#'
#' @export
#'
#' @return A ggplot2 object.


clustree_viz <- function(HapObject, type = "MG") {
#Prepare epsilon values
epsilon <- NULL

for (x in 1:length(HapObject)){
epsilon <- base::append(epsilon, HapObject[[x]][["epsilon"]])
}

#Extract ID file first epsilon value and change column name to hap_epsXX
pre_clustree <- HapObject[[1]]$Indfile %>%
  dplyr::rename(!!base::paste0("hap_eps", HapObject[[1]]$epsilon) := 'hap') %>%
  dplyr::select(-("Metadata"))

#Iterate over all other epsilon values, adding hap_epsXX columns to tibble
for (ecks in 2:length(HapObject)){
  pre_clustree <- dplyr::left_join(pre_clustree,
                     HapObject[[ecks]][["Indfile"]] %>%
    dplyr::select(-("Metadata")) %>%
    dplyr::rename(!!base::paste0("hap_eps", HapObject[[ecks]]$epsilon) := 'hap'), by = c("Ind", "Pheno")
)
}

#Plot with clustree
haptree <- base::suppressMessages(
  clustree::clustree(pre_clustree, prefix = 'hap_eps', node_colour = 'Pheno',node_colour_aggr = "mean_na.rm", edge_width = 1, node_alpha = 1)+
  ggplot2::scale_colour_gradient(limits=c(base::max(dplyr::top_frac(pre_clustree,-0.1,.data$Pheno)$Pheno),
                                 base::min(dplyr::top_frac(pre_clustree,0.1,.data$Pheno)$Pheno)),
                        high = "#8ADD81",low = "#6870F6", oob = scales::squish,name = 'Pheno') +
  ggplot2::scale_colour_gradient(high = 'black', low = 'grey80', guide = "edge_colourbar", aesthetics = "edge_colour") +
  ggplot2::labs(size = 'nIndividuals', edge_alpha = "Proportion") +
  ggplot2::guides(edge_color = "none", size = ggplot2::guide_legend(order = 1))
)
#Create epsilon label data
#Extract x and y coordinates from clustree object and build labels
haplbls <-
  tibble::tibble(xval = base::max(haptree[["data"]][["x"]])*1.1,
         yval=(base::length(epsilon)-1):0,
         labelval = base::paste0("\u03b5"," = ",epsilon))

#Re-plot with epsilon label data added
labeled_haptree <- haptree +
  ggplot2::geom_text(data = haplbls, ggplot2::aes(x=.data$xval, y=.data$yval, label=.data$labelval), hjust = 0)

#Repeat with MGs rather than haplotype groups

pre_MGtree <- HapObject[[1]]$Varfile %>%
  dplyr::select(c("ID", "phenodiff", "MGs")) %>%
  dplyr::rename(!!base::paste0("MGs_eps", HapObject[[1]]$epsilon) := 'MGs')

for (ecks in 2:length(HapObject)){
  pre_MGtree <- dplyr::left_join(pre_MGtree,
                                 HapObject[[ecks]][["Varfile"]] %>%
                                   dplyr::select("ID", "MGs")%>%
                                   dplyr::rename(!!base::paste0("MGs_eps",HapObject[[ecks]]$epsilon) := 'MGs'), by = "ID")
}

MGtree <- base::suppressMessages(
clustree::clustree(pre_MGtree, prefix = 'MGs_eps', node_colour = 'phenodiff',node_colour_aggr = "mean_na.rm", edge_width = 1, node_alpha = 1) +
  ggplot2::scale_colour_gradient(limits=c(base::max(pre_MGtree$phenodiff),
                                          base::min(pre_MGtree$phenodiff)),
                                 high = "#8ADD81",low = "#6870F6",oob = scales::squish,name = 'phenodiff') +
  ggplot2::scale_colour_gradient(high = 'black', low = 'grey80', guide = "edge_colourbar", aesthetics = "edge_colour") +
  ggplot2::labs(size = 'nSNPs', edge_alpha = "Proportion") +
  ggplot2::guides(edge_color = "none", size = ggplot2::guide_legend(order = 1))
)
#Create epsilon label data
#Extract x and ay coordinates from clustree object and build labels
MGlbls <-
  tibble::tibble(xval = base::max(MGtree[["data"]][["x"]])*1.1,
                 yval=(base::length(epsilon)-1):0,
                 labelval = base::paste0("\u03b5"," = ",epsilon))

#Re-plot with epsilon label data added
labeled_MGtree <- MGtree +
  ggplot2::geom_text(data = MGlbls, ggplot2::aes(x=.data$xval, y=.data$yval, label=.data$labelval), hjust = 0)

return(switch(type,
       "hap" = labeled_haptree,
       "MG" = labeled_MGtree))
}

