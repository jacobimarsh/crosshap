#' Run clustering on data
#'
#' @param elon Epsilon values to search.
#' @param pheno Input numeric phenotype data for colouring clustree.
#'
#' @return
#' @export
#'
#' @examples
#'
run_clustree <- function(elon, pheno) {
pre_clustree <- get(paste("Haplotype_assignments_MP","_E",elon[1],sep=""))

for (drez in elon[2:length(elon)]){
pre_clustree <- pre_clustree %>%
    left_join(get(paste("Haplotype_assignments_MP","_E",drez,sep="")))
}

pre_clustree_phen <- left_join(pre_clustree, pheno)

#plot
ctree <-
  clustree(pre_clustree_phen, prefix = 'hap_eps', node_colour = 'Pheno',node_colour_aggr = 'mean_na.rm', edge_width = 1, node_alpha = 1)+
  scale_colour_gradient(limits=c(max(top_frac(pre_clustree_phen,-0.1,Pheno)$Pheno),
                                 min(top_frac(pre_clustree_phen,0.1,Pheno)$Pheno)),
                        high = "#8ADD81",low = "#6870F6",oob = scales::squish,name = 'Pheno') +
  scale_edge_color_continuous(high = 'black',low = 'grey80') +
  labs(size = 'nIndividuals', edge_alpha = "Proportion") +
  guides(edge_color = "none", size = guide_legend(order = 1))

#create label data
#x and y extractd from clustree object
lbls <-
  tibble(xval = max(ctree[["data"]][["x"]])*1.1,
         yval=0:(length(elon)-1),
         labelval = paste0("\u03b5"," = ",elon))

#re-plot using label data
labeled_ctree <- ctree +
  geom_text(data = lbls, aes(x=xval, y=yval, label=labelval), hjust = 0)

return(labeled_ctree)
}
