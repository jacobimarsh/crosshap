#' UMAP haplotype visualization helper
#'
#' prepare_hap_umap() builds a large composite ggplot2 object ready for faceting
#' and animation (see examples) for visualizing SNP alleles (coloured by marker
#' group) possessed by randomly sampled individuals within each haplotype. UMAP
#' coordinates for each SNP can be generated using umap::umap(), with the LD
#' matrix generated for run_haplotyping() as input. When fully rendered and
#' faceted, the resultant GIF intuitively visualizes the consistent alleles for
#' marker groups SNPs for each haplotype combination, despite variation across
#' non-marker group SNP alleles.
#'
#' @param umap_in Epsilon values passed through run_haplotyping().
#' @param vcf Input vcf
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param nsamples Number of times to sample each haplotype group, will directly
#' translate to the number of frames in animation. Should be the same as the
#' nframes passed to gganimate::animate().
#' @param hetmiss_as If hetmiss_as = "allele", heterozygous-missing SNPs './N'
#' are recoded as 'N/N', if hetmiss_as = "miss", the site is recoded as missing.
#' @param het_as If het_as = "alt", heterozygous SNPs are treated as 'ALT/ALT'
#' alleles, if het_as = "ref", if het_as = "ref", they are kept as 'REF/REF'.
#'
#' @export
#'
#' @return A large ggplot2 object.
#'
#' @examples hap_umap <- umap::umap(LD, min_dist = 2, spread = 2.5, n_neighbors = MGmin)
#'
#' hap_gg <- prepare_hap_umap(hap_umap, HapObject = Haplotypes_MGmin30_E0.6, vcf = vcf, nsamples = 25)
#'
#' hap_gganim <- hap_gg +
#'   ggplot2::facet_wrap(~hap) +
#'   gganimate::transition_states(Frame, transition_length = 0, state_length = 1)
#'
#' hap_anim <- gganimate::animate(hap_gganim, renderer = gganimate::gifski_renderer(), nframes = 25)
#'
#' gganimate::anim_save("hap_anim.gif", hap_anim)
#'


#umap_in <- umap(protLD, min_dist = 2, spread = 2.5, n_neighbors = MGmin)

prepare_hap_umap <- function(umap_in, hetmiss_as = 'allele', het_as = "alt", HapObject, vcf, nsamples = 25){

ID_bin_vcf <- dplyr::select(vcf, -c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
  dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,
                                             base::ifelse(x=='1|0'|x=='0|1',1,
                                                          base::ifelse(x=='1|1',2,
                                                                       switch(hetmiss_as, "allele" = base::ifelse(x=='1|.'|x=='.|1',1,
                                                                                                                  base::ifelse(x=='0|.'|x=='.|0',0,NA)),
                                                                              "miss" = NA))))}) %>%
  tibble::rownames_to_column('ID')

x_y_MGfile <-
  umap_in$layout %>% dplyr::as_tibble() %>%
  cbind(rownames(umap_in$data)) %>% dplyr::as_tibble() %>%
  dplyr::rename(UMAP1 = V1, UMAP2 = V2, ID = 'rownames(umap_in$data)') %>%
  dplyr::left_join(HapObject$MGfile %>%
                     dplyr::select(ID, MGs, cluster), by = "ID") %>%
  dplyr::mutate(cluster = as.character(cluster))

x_y_vcf <- dplyr::left_join(x_y_MGfile, ID_bin_vcf, by = "ID")

x_y_vcf[x_y_vcf == "0"] <- NA

x_y_vcf_long <-
  x_y_vcf %>% dplyr::select(UMAP1, UMAP2, MGs, cluster, HapObject$Indfile$Ind) %>%
  tidyr::gather(HapObject$Indfile$Ind, key = "Ind", value = "allele") %>%
  dplyr::left_join(HapObject$Indfile, hap, by = "Ind")

framenum_ID <- NULL
for (i in c(ifelse("0" %in% HapObject$Indfile$hap, 0, NULL),HapObject$Hapfile$hap)){
  framenum_ID <- rbind(framenum_ID, cbind(1:nsamples, dplyr::filter(HapObject$Indfile, hap == i)$Ind %>%
                                      sample(size = nsamples, replace = T)) %>% dplyr::as_tibble())
}
framenum_ID <- framenum_ID %>% dplyr::rename(Frame = V1, Ind = V2)

## xyvl_framenum <- x_y_vcf_long %>% dplyr::filter(Ind %in% framenum_ID$Ind) %>%
##   dplyr::left_join(framenum_ID, by = 'Ind') %>%
##   dplyr::mutate(hap = ifelse(hap==0,'Unassigned',paste0('Hap ', hap))) %>%
##   dplyr::mutate(MG_cols = ifelse(!is.na(allele), MGs, NA)) %>%
##   dplyr::mutate(MG_cols = ifelse(!is.na(allele) & is.na(MGs),0,MG_cols)) %>%
##   dplyr::mutate(MG_cols = ifelse(MG_cols == '0','Non-MG (0)',MG_cols))

xyvl_framenum <- x_y_vcf_long %>% dplyr::filter(Ind %in% framenum_ID$Ind) %>%
  dplyr::left_join(framenum_ID, by = 'Ind') %>%
  dplyr::mutate(hap = ifelse(hap==0,'Unassigned',paste0('Hap ', hap))) %>%
  dplyr::mutate(MG_cols = ifelse(allele == ifelse(het_as == "ref",2,c(1,2)), MGs, NA)) %>%
  dplyr::mutate(MG_cols = ifelse(allele == ifelse(het_as == "ref",2,c(1,2)) & is.na(MGs),0,MG_cols)) %>%
  dplyr::mutate(MG_cols = ifelse(MG_cols == '0','Non-MG (0)',MG_cols))


pre_anim_gg <- ggplot2::ggplot(xyvl_framenum %>% dplyr::arrange(is.na(MG_cols), dplyr::desc(MG_cols)) , ggplot2::aes(x = UMAP1, y = UMAP2)) +
  ggplot2::geom_point(alpha = 0.4, ggplot2::aes(colour = MG_cols), size = 0.25)  +
  ggplot2::geom_text(data = framenum_ID %>% dplyr::mutate(allele = 1:nrow(framenum_ID)) %>%
                       dplyr::left_join(HapObject$Indfile, by = "Ind") %>%
                       dplyr::mutate(hap = ifelse(hap==0,'Unassigned',paste0('Hap ', hap))),
                     mapping = ggplot2::aes(x = -Inf, y = Inf, label = Ind), hjust = 0, vjust = 1, size = 3) +
  ggplot2::theme_void() +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 10)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 5, alpha = 0.7), title = "Alt allele"))
return(pre_anim_gg)
}
