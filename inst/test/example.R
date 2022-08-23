##Example
library(crosshap)

"hello" %>% print()

LD <- read_LD('~/Desktop/bash_misc/crosshap_data/data/LD_100kb.mtx')

vcf <- read_vcf('~/Desktop/bash_misc/crosshap_data/data/impu_100kb_pdh1.vcf')
phen_early <- read_pheno('~/Desktop/bash_misc/crosshap_data/data/early_shatter565.txt')
prot <- read_pheno('~/Desktop/bash_misc/crosshap_data/data/prot_phen.txt')

##DBscan epsilon values to be tested and compared using clustree visualization
epsilon <- seq(2.5,3,by=0.5)
MGmin <- 10
minHap <- 9

run_haplotyping(vcf, LD, phen_early, epsilon, MGmin, minHap)

##Visualize differences between clusters based on epsilon value input to DBscan
labeled_ctree <- run_clustree(epsilon, MGmin, phen_early)

labeled_ctree

Hap2_viz <- crosshap_viz(Haplotypes_MGmin10_E3)

ggsave("clustree_MPx2.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')

