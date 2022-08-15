##Example

setwd("/Users/jmarsh96/Desktop/bash_misc/crosshap_data")

LD <- read_LD('data/LD_100kb.mtx')
vcf <- read_vcf('data/impu_100kb_pdh1.vcf')
phen_early <- read_pheno('data/early_shatter565.txt')
prot <- read_pheno('data/prot_phen.txt')

##DBscan epsilon values to be tested and compared using clustree visualization
epsilon <- seq(1,3,by=0.5)

#run_haplotyping <- function(vcf, LD, pheno, epsilon) {
##Perform clustering and count haplotype frequencies
for (arez in epsilon){
  db40 <- dbscan::dbscan(LD, eps = arez, minPts = 40)
  MGfile <- tibble(POS=rownames(LD),cluster=db40$cluster)

##Call allelic states for each SNP marker group across individuals
  het_pseudoSNP <- create_pseudoSNP(MGfile,vcf)

##Identify haplotype frequencies for different marker group combinations
  clustered_hpS <- pseudo2haps(het_pseudoSNP)
  Varfile <- run_clustree(MGfile, vcf, phen_early)
  clustered_hpS_obj <-  list(epsilon = db40$eps,
                             minPts = db40$minPts,
                             Hapfile = clustered_hpS$Hapfile,
                             IDfile = left_join(clustered_hpS$nophenIDfile,phen_early),
                             Varfile = Varfile,
                             MGfile = MGfile)
  assign(paste("Haplotypes_MP", "_E", arez,sep = ""), clustered_hpS_obj)
}

#  return(get(paste("Haplotypes_MP", "_E", arez,sep = "")))

##Visualize differences between clusters based on epsilon value input to DBscan
labeled_ctree <- run_clustree(epsilon,phen_early)
ggsave("clustree_MPx2.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')

labeled_ctree
