##Example

setwd("/Users/jmarsh96/Desktop/bash_misc/crosshap_data")

LD <- read_LD('data/LD_100kb.mtx')
vcf <- read_vcf('data/impu_100kb_pdh1.vcf')
phen_early <- read_pheno('data/early_shatter565.txt')
prot <- read_pheno('data/prot_phen.txt')

##DBscan epsilon values to be tested and compared using clustree visualization
epsilon <- seq(0.69,2.69,by=0.5)
MGmin <- 40

run_haplotyping(vcf, LD, phen_early, epsilon, MGmin)

##Visualize differences between clusters based on epsilon value input to DBscan
labeled_ctree <- run_clustree(epsilon,phen_early)
ggsave("clustree_MPx2.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')

labeled_ctree
