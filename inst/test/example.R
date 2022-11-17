##Example
library(crosshap)
#library(patchwork)

"hello" %>% print()

LD <- crosshap::read_LD('~/Desktop/bash_misc/crosshap_data/data/LD_100kb.mtx')

vcf <- crosshap::read_vcf('~/Desktop/bash_misc/crosshap_data/data/impu_100kb_pdh1.vcf')
phen_early <- crosshap::read_pheno('~/Desktop/bash_misc/crosshap_data/data/early_shatter565.txt')
prot <- crosshap::read_pheno('~/Desktop/bash_misc/crosshap_data/data/prot_phen.txt')

##DBscan epsilon values to be tested and compared using clustree visualization
epsilon <- seq(.5,1,by=0.5)
MGmin <- 30
minHap <- 9



crosshap::run_haplotyping(epsilon = epsilon, vcf = vcf, LD = LD, pheno = phen_early, MGmin = MGmin, minHap = minHap, )
##Visualize differences between clusters based on epsilon value input to DBscan
  labeled_haptree <- crosshap::run_clustree(MGmin = MGmin, pheno = pheno, type = "hap")
labeled_MGtree <- crosshap::run_clustree(MGmin = MGmin, pheno = phen_early)

Hap2_viz <- crosshap::crosshap_viz(Haplotypes_MGmin30_E1, plot_left = "pos")
Hap2_labs_viz <- crosshap::crosshap_viz(Haplotypes_MGmin30_E1, hide_labels = F)

#ggplot2::ggsave("botrightlabs.pdf",crosshap_stitched,device = 'pdf',dpi = 300,height = 9,width = 16,units = 'in')


##nonImpu

nonimpu_vcf <- crosshap::read_vcf('~/Desktop/bash_misc/crosshap_data/data/headed_100kb_pdh1.vcf')



vcf <- crosshap::read_vcf('~/Desktop/bash_misc/crosshap_data/data/dummy_test.vcf')


#prot

LD <- crosshap::read_LD("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/LD_173kb.mtx")
pheno <- crosshap::read_pheno("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/prot_phen.txt")
vcf <- crosshap::read_vcf("/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/fin_b51_173kb_only.vcf")
metadata <- crosshap::read_metadata('/Users/jmarsh96/Desktop/bash_misc/crosshap_data/data/labmeeting/namepopfile.txt')

eps <- seq(.2,1,by=.2)

crosshap::run_haplotyping(vcf = vcf,
                LD = LD,
                metadata = metadata,
                pheno = pheno,MGmin = 30, minHap = 9)
                #


prot_clustree <- run_clustree(epsilon = eps,
                              MGmin = 30,
                              pheno = prot_phen,
                              type = 'hap')

prot_viz <- crosshap::crosshap_viz(Haplotypes_MGmin30_E0.6, hide_labels = F, plot_right = "cluster")

posplot_prot_viz <- crosshap_viz(Haplotypes_MGmin30_E1, hide_labels = F, plot_left = "pos", plot_right = "cluster")

hdbposplot_prot_viz <- crosshap_viz(Haplotypes_MGmin30_EX, hide_labels = F, plot_left = "pos", plot_right = "cluster")

run_hdbscan_haplotyping(vcf = prot_vcf,
                        LD = protLD,
                        pheno = prot_phen,
                        metadata = metadata,
                        MGmin = 30,
                        minHap =9,
                        keep_outliers = F
)

hdbposplot_prot_viz <- crosshap_viz(Haplotypes_MGmin30_HDBSCAN, hide_labels = F, plot_left = "allele", plot_right = "cluster")


run_haplotyping(vcf = prot_vcf,
                LD = protLD,
                pheno = prot_phen,
                metadata = metadata,
                epsilon = eps,
                MGmin = 30,
                minHap = 9,
)




posplot_prot_viz <- crosshap_viz(Haplotypes_MGmin30_E1, hide_labels = F, plot_left = "allele", plot_right = "cluster")

tprot <- tsne(protLD)

pca_protLD <- prcomp(protLD)

tprot_labeled <- tprot %>% rename("X" = "V3", "Y" = "V4","ID" = "POS") %>% mutate(Y = as.numeric(Y), X = as.numeric(X))

tsne_prot_MG40_E2 <- Haplotypes_MGmin30_E0.6$MGfile %>% select(-POS) %>% left_join(tprot_labeled)

ggplot(Haplotypes_MGmin40_E2.5$MGfile %>% select(-POS) %>%
         left_join(tprot_labeled), aes(X, Y)) +
  geom_point(aes(colour = factor(cluster)))





dbE1 <- build_right_clusterplot(Haplotypes_MGmin30_E0.6, hide_labels = T)

Haplotypes_MGmin30_E1$MGfile %>% group_by(MGs) %>% summarise(groupvar = var(meanr2),
                                                             count = length(x = POS))

hdb <- build_right_clusterplot(Haplotypes_MGmin30_HDBSCAN, hide_labels = T)

Haplotypes_MGmin30_HDBSCAN$MGfile %>% group_by(MGs) %>% summarise(groupvar = var(meanr2),
                                                             count = length(x = POS))





Haplotypes_MGmin30_E1$MGfile %>% group_by(MGs) %>%
#  filter(!(abs(meanr2 - median(meanr2)) > 2*sd(meanr2))) %>%
  summarise_each(mean, meanr2)

Haplotypes_MGmin30_E0.4$MGfile %>% group_by(MGs) %>%
#  filter(!(abs(meanr2 - median(meanr2)) > 2*sd(meanr2))) %>%
  summarise_each(mean, meanr2)


df1 = df %>%
  group_by(element) %>%
  filter(!(abs(value - median(value)) > 2*sd(value))) %>%
  summarise_each(funs(mean), value)


dbE1 <- build_right_clusterplot(Haplotypes_MGmin30_E0.4, hide_labels = T)

smoothed <- build_right_clusterplot(Haplotypes_MGmin30_E0.4, hide_labels = T)

smoothed <- build_right_clusterplot(Haplotypes_MGmin30_HDBSCAN, hide_labels = T)

outs_jit

noouts_jit


umap_in2 <- umap::umap(LD, min_dist = 2, spread = 2.5, n_neighbors = 30)

pre_anim_gg <- prepare_hap_umap(umap_in2,
                            HapObject = Haplotypes_MGmin30_E0.6,
                            vcf = vcf,
                            nsamples = 25)

hap_gganim <- pre_anim_gg +
  ggplot2::facet_wrap(~hap)+
  gganimate::transition_states(Frame,
                    transition_length = 0,
                    state_length = 1)

gganimate::animate(
  hap_gganim,
  renderer = gganimate::gifski_renderer(),
  fps = 3,
  res = 300,
  width  = 6,
  height = 6,
  units = "in",
  nframes = 25
)

colnames(vcf) <- gsub('-','.',colnames(vcf))

gganimate::anim_save("hap_anim_cols.gif", hap_anim)












layout <- "DAE
FFF"

#base::message(paste0("Stitching plots"))
crosshap_stitched <-
  patchwork::wrap_plots(mid, left, right) +
  patchwork::guide_area() +
  patchwork::plot_layout(design = layout, guides = "collect")


