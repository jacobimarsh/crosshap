##Example
library(crosshap)
library(patchwork)

"hello" %>% print()

LD <- crosshap::read_LD('~/Desktop/bash_misc/crosshap_data/data/LD_100kb.mtx')

vcf <- crosshap::read_vcf('~/Desktop/bash_misc/crosshap_data/data/impu_100kb_pdh1.vcf')
phen_early <- crosshap::read_pheno('~/Desktop/bash_misc/crosshap_data/data/early_shatter565.txt')
prot <- crosshap::read_pheno('~/Desktop/bash_misc/crosshap_data/data/prot_phen.txt')

##DBscan epsilon values to be tested and compared using clustree visualization
epsilon <- seq(1,3,by=0.5)
MGmin <- 40
minHap <- 9

crosshap::run_haplotyping(vcf = vcf, LD = LD, pheno = prot, MGmin = MGmin, minHap = minHap, epsilon = epsilon)
##Visualize differences between clusters based on epsilon value input to DBscan
labeled_ctree <- crosshap::run_clustree(epsilon, MGmin, prot)

labeled_ctree

Hap2_viz <- crosshap::crosshap_viz(Haplotypes_MGmin40_E2.5)

ggsave("clustree_MPx2.pdf",labeled_ctree,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')

