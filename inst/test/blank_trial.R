library(tidyverse)
library(patchwork)
library(crosshap)

##Pod shatter

pod_vcf <- read_vcf("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/impu_100kb_pdh1.vcf")
podLD <- read_LD("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/LD_100kb.mtx")
pod_phen <- read_pheno("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/early_shatter565.txt")

run_haplotyping(vcf = pod_vcf,
                LD = podLD,
                pheno = pod_phen,
                epsilon = eps,
                MGmin = 30,
                minHap = 9)

pod_clustree <- run_clustree(epsilon = eps,
                              MGmin = 30,
                              pheno = pod_phen)

pod_viz <- crosshap_viz(Haplotypes_MGmin30_E2)

tpod <- tsne(podLD)

tpod_labeled <- tpod %>% as_tibble() %>% cbind(rownames(podLD)) %>% rename("X" = "V1", "Y" = "V2","ID" = "rownames(podLD)") %>% mutate(Y = as.numeric(Y), X = as.numeric(X))

tsne_pod_MG40_E2 <- Haplotypes_MGmin30_E2$MGfile %>% select(-POS) %>% left_join(tpod_labeled)

ggplot(Haplotypes_MGmin30_E1.5$MGfile %>% select(-POS) %>%
         left_join(tpod_labeled, by = "ID"), aes(X, Y)) +
  geom_point(aes(colour = factor(cluster)))

##Protein
protLD <- read_LD("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/LD_173kb.mtx")

prot_phen <- read_pheno("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/prot_phen.txt")
prot_vcf <- read_vcf("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/fin_b51_173kb_only.vcf")

eps <- seq(0.5,3,by=0.5)

run_haplotyping(vcf = prot_vcf,
                LD = protLD,
                pheno = prot_phen,
                epsilon = eps,
                MGmin = 30,
                minHap = 9)

prot_clustree <- run_clustree(epsilon = eps,
                              MGmin = 30,
                              pheno = prot_phen)

prot_viz <- crosshap_viz(Haplotypes_MGmin30_E0.5, plot_left = "allele")

tprot <- tsne(protLD)

pca_protLD <- prcomp(protLD)

tprot_labeled <- tprot %>% rename("X" = "V3", "Y" = "V4","ID" = "POS") %>% mutate(Y = as.numeric(Y), X = as.numeric(X))

tsne_prot_MG40_E2 <- Haplotypes_MGmin40_E2$MGfile %>% select(-POS) %>% left_join(tprot_labeled)

ggplot(Haplotypes_MGmin40_E2.5$MGfile %>% select(-POS) %>%
         left_join(tprot_labeled), aes(X, Y)) +
  geom_point(aes(colour = factor(cluster)))
