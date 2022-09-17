##Example
library(crosshap)
#library(patchwork)

"hello" %>% print()

LD <- crosshap::read_LD('~/Desktop/bash_misc/crosshap_data/data/LD_100kb.mtx')

vcf <- crosshap::read_vcf('~/Desktop/bash_misc/crosshap_data/data/impu_100kb_pdh1.vcf')
phen_early <- crosshap::read_pheno('~/Desktop/bash_misc/crosshap_data/data/early_shatter565.txt')
prot <- crosshap::read_pheno('~/Desktop/bash_misc/crosshap_data/data/prot_phen.txt')

##DBscan epsilon values to be tested and compared using clustree visualization
epsilon <- seq(1.5,2,by=0.5)
MGmin <- 30
minHap <- 9

crosshap::run_haplotyping(vcf = vcf, LD = LD, pheno = phen_early, MGmin = MGmin, minHap = minHap, hetmiss_as = "miss")
##Visualize differences between clusters based on epsilon value input to DBscan
labeled_haptree <- crosshap::run_clustree(MGmin = MGmin, pheno = phen_early, type = "hap")
labeled_MGtree <- crosshap::run_clustree(MGmin = MGmin, pheno = phen_early)

Hap2_viz <- crosshap_viz(Haplotypes_MGmin30_E2, plot_left = "pos")
Hap2_labs_viz <- crosshap_viz(Haplotypes_MGmin30_E2, hide_labels = F)

ggplot2::ggsave("botrightlabs.pdf",crosshap_stitched,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')


##nonImpu

nonimpu_vcf <- read_vcf('~/Desktop/bash_misc/crosshap_data/data/headed_100kb_pdh1.vcf')



vcf <- read_vcf('~/Desktop/bash_misc/crosshap_data/data/dummy_test.vcf')


#prot

protLD <- crosshap::read_LD("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/LD_173kb.mtx")
prot_phen <- crosshap::read_pheno("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/prot_phen.txt")
prot_vcf <- crosshap::read_vcf("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/fin_b51_173kb_only.vcf")
metadata <- crosshap::read_metadata('/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/namepopfile.txt')

eps <- seq(.2,1,by=.2)

crosshap::run_haplotyping(vcf = prot_vcf,
                LD = protLD,
                pheno = prot_phen,
                metadata = metadata,
                epsilon = eps,
                MGmin = 30,
                minHap = 9)

prot_clustree <- crosshap::run_clustree(epsilon = eps,
                              MGmin = 30,
                              pheno = prot_phen)

prot_viz <- crosshap::crosshap_viz(Haplotypes_MGmin30_E0.6, hide_labels = F)

#tprot <- tsne(protLD)
#
#pca_protLD <- prcomp(protLD)
#
#tprot_labeled <- tprot %>% rename("X" = "V3", "Y" = "V4","ID" = "POS") %>% mutate(Y = as.numeric(Y), X = as.numeric(X))
#
#tsne_prot_MG40_E2 <- Haplotypes_MGmin40_E2$MGfile %>% select(-POS) %>% left_join(tprot_labeled)
#
#ggplot(Haplotypes_MGmin40_E2.5$MGfile %>% select(-POS) %>%
#         left_join(tprot_labeled), aes(X, Y)) +
#  geom_point(aes(colour = factor(cluster)))

