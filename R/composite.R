##Perform clustering and count haplotype frequencies
for (arez in elon){
  db40 <- dbscan::dbscan(LD, eps = arez, minPts = 40)
  db40_clust <- tibble(POS=rownames(LD),cluster=db40$cluster)
  assign(paste("export_clusts", arez,sep = ""), db40_clust)

##Call allelic states for each SNP marker group across individuals
  het_pseudoSNP <- create_pseudoSNP(db40_clust,vcf)
  assign(paste("het_pseudoSNP_MPx", "_E", arez,sep = ""), het_pseudoSNP)

##Identify haplotype frequencies for different marker group combinations
  clustered_hpS <- pseudo2haps(het_pseudoSNP)
  assign(paste("Haplotype_assignments_MP", "_E", arez,sep = ""), clustered_hpS)
}

##Visualize differences between clusters based on epsilon value input to DBscan
labeled_ctree <- run_clustree(elon, phen_early)
ggsave("clustree_MPx.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')


