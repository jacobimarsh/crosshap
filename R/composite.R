---
  title: "CrossHap_clustering"
date: '2022-07-19'
---

#conflict_prefer_all('dplyr')
setwd("/Users/jmarsh96/Desktop/bash_misc/crosshap_data")

#Read correlation matrix between all SNPs in region
LD <- fread('data/LD_100kb.mtx', nThread = 10) %>%  as_tibble() %>%  column_to_rownames("V1")

##Parse a VCF as a matrix and convert alleles to base 3 integers
vcf <- fread('data/impu_100kb_pdh1.vcf', nThread = 10) %>%  as_tibble() %>%
  select(-c(1,2,4:9)) %>% column_to_rownames('ID') %>%
  mutate_all(function(x){ifelse(x=='0|0',0,ifelse(x=='1|0'|x=='0|1',1,ifelse(x=='1|1',2,'failsave')))})

phen_early <- fread('data/early_shatter565.txt') %>% as_tibble() %>% rename('ID' = V1, 'Pheno' = V2)

elon <- seq(0.9,1.5,by=0.2)

##Clustering loop

for (arez in elon){
  db40 <- dbscan::dbscan(LD, eps = arez, minPts = 40)
  db40_clust <- tibble(POS=rownames(LD),cluster=db40$cluster)
  assign(paste("export_clusts", arez,sep = ""), db40_clust)

  het_pseudoSNP <- create_pseudoSNP(db40_clust,vcf)
  assign(paste("het_pseudoSNP_MPx", "_E", arez,sep = ""), het_pseudoSNP)

  clustered_hpS <- pseudo2haps(het_pseudoSNP)
  assign(paste("Haplotype_assignments_MP", "_E", arez,sep = ""), clustered_hpS)
}

#clustree
labeled_ctree <- run_clustree(elon, phen_early)
ggsave("clustree_MPx.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')


