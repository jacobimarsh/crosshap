#' UMAP haplotype visualization helper
#'
#' prepare_hap_umap() builds a large composite ggplot2 object ready for faceting
#' and animation (see vignette) for visualizing SNP alleles (coloured by Marker
#' Group) possessed by individuals with each haplotype. UMAP
#' coordinates for each SNP can be generated using umap::umap(), with the LD
#' matrix generated for run_haplotyping() as input. When fully rendered and
#' faceted, the resultant GIF intuitively visualizes the shared loci within each
#' Marker Group that are constant within each haplotype combination.
#'
#' @param umap_in UMAP results produced for a haplotype object at a given epsilon.
#' @param vcf Input vcf.
#' @param HapObject Haplotype object created by run_haplotyping().
#' @param epsilon Epsilon matching the haplotype object used for umap_in.
#' @param nsamples Number of times to sample each haplotype group, will directly
#' translate to the number of frames in animation. Should be the same as the
#' nframes passed to gganimate::animate().
#' @param hetmiss_as If hetmiss_as = "allele", heterozygous-missing SNPs './N'
#' are recoded as 'N/N', if hetmiss_as = "miss", the site is recoded as missing.
#'
#' @importFrom rlang ".data"
#'
#' @export
#'
#' @return A large ggplot2 object.

prepare_hap_umap <- function(umap_in, hetmiss_as = 'allele', HapObject, epsilon, vcf, nsamples = 25){

#Extract haplotype results for given epsilon
for (x in 1:length(HapObject)){
    if(HapObject[[x]]$epsilon == epsilon){
      HapObject_eps <- HapObject[[x]]
    }
  }

ID_bin_vcf <- dplyr::select(vcf, -c(1,2,4:9)) %>% tibble::column_to_rownames('ID') %>%
  dplyr::mutate_all(function(x){base::ifelse(x=='0|0',0,
                                             base::ifelse(x=='1|0'|x=='0|1',1,
                                                          base::ifelse(x=='1|1',2,
                                                                       switch(hetmiss_as, "allele" = base::ifelse(x=='1|.'|x=='.|1',1,
                                                                                                                  base::ifelse(x=='0|.'|x=='.|0',0,NA)),
                                                                              "miss" = NA))))}) %>%
  tibble::rownames_to_column('ID')

#Attach UMAP X and Y coordinates to Varfile (Marker Group assignments)
x_y_Varfile <-
  suppressWarnings(umap_in$layout %>% dplyr::as_tibble() %>%
  cbind(rownames(umap_in$data)) %>% dplyr::as_tibble(.name_repair = 'unique') %>%
  dplyr::rename('UMAP1' = 'V1', 'UMAP2' = 'V2', 'ID' = 'rownames(umap_in$data)') %>%
  dplyr::left_join(HapObject_eps$Varfile %>%
                     dplyr::select('ID', 'MGs', 'cluster'), by = "ID") %>%
  dplyr::mutate(cluster = as.character(.data$cluster)))

#Attach to VCF with allele states for each individual
x_y_vcf <- dplyr::left_join(x_y_Varfile, ID_bin_vcf, by = "ID")

x_y_vcf[x_y_vcf == "0"] <- NA

#Convert to long format (warning: can be an extremely large tibble)
x_y_vcf_long <-
  x_y_vcf %>% dplyr::select('UMAP1', 'UMAP2', 'MGs', 'cluster', HapObject_eps$Indfile$Ind) %>%
  tidyr::gather(HapObject_eps$Indfile$Ind, key = "Ind", value = "allele") %>%
  dplyr::left_join(HapObject_eps$Indfile, by = "Ind")

#
framenum_ID <- NULL
#For each haplotype, randomly sample individuals, and assign them each a frame number
for (i in c(ifelse("0" %in% HapObject_eps$Indfile$hap, 0, NULL),HapObject_eps$Hapfile$hap)){
  framenum_ID <- rbind(framenum_ID, cbind(1:nsamples, dplyr::filter(HapObject_eps$Indfile, .data$hap == i)$Ind %>%
                                      sample(size = nsamples, replace = T)) %>% dplyr::as_tibble())
}
framenum_ID <- framenum_ID %>% dplyr::rename('Frame' = 'V1', 'Ind' = 'V2')

#Attach UMAP X & Y with VCF information to the frame number individual information
xyvl_framenum <- suppressWarnings(
  x_y_vcf_long %>% dplyr::filter(.data$Ind %in% framenum_ID$Ind) %>%
  dplyr::left_join(framenum_ID, by = 'Ind') %>%
  dplyr::mutate(hap = ifelse(.data$hap==0,'Unassigned',paste0('Hap ', .data$hap))) %>%
  dplyr::mutate(MG_cols = ifelse(!is.na(.data$allele), .data$MGs, NA)) %>%
  dplyr::mutate(MG_cols = ifelse(!is.na(.data$allele) & is.na(.data$MGs),0,.data$MG_cols)) %>%
  dplyr::mutate(MG_cols = ifelse(.data$MG_cols == '0','Non-MG (0)',.data$MG_cols)))

#Add dummy duplicate of each frame as hap A with all alternate alleles for the SNPs
#placed in Marker Groups to have a static 'All MGs' frame in the top left.
MGrows <- xyvl_framenum %>%
  dplyr::filter(.data$hap == "Hap A") %>%
  dplyr::mutate(hap = "All MGs", Ind = "N/A", MG_cols = .data$MGs)

#Add it back to the composite tibble
xyvlf_MG <- xyvl_framenum %>% rbind(MGrows)

#Build a large ggplot object with all frames glued on top of one another
#From this you can facet by haplotype and animate.
pre_anim_gg <- ggplot2::ggplot(xyvlf_MG %>% dplyr::arrange(is.na(.data$MG_cols), dplyr::desc(.data$MG_cols)) , ggplot2::aes(x = .data$UMAP1, y = .data$UMAP2)) +
  ggplot2::geom_point(alpha = 0.4, ggplot2::aes(colour = .data$MG_cols), size = 0.25)  +
  ggplot2::geom_text(data = framenum_ID %>% dplyr::mutate(allele = 1:nrow(framenum_ID)) %>%
                       dplyr::left_join(HapObject_eps$Indfile, by = "Ind") %>%
                       dplyr::mutate(hap = ifelse(.data$hap==0,'Unassigned',paste0('Hap ', .data$hap))),
                     mapping = ggplot2::aes(x = 1.05*min(xyvl_framenum$UMAP2), y = 1.05*max(xyvl_framenum$UMAP1), label = .data$Ind), hjust = 0, vjust = 1, size = 2.5) +
  ggplot2::theme_void() +
  ggplot2::theme(panel.border = ggplot2::element_rect(colour = 'grey90', linewidth = 1, fill = NA)) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 10)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 5, alpha = 0.7), title = "Alt allele")) +
  ggplot2::scale_color_manual(labels = c(paste0("MG",1:(length(unique(xyvl_framenum$MG_cols))-2)), "Non-MG (0)", "REF"),
                                values = c(c("#F763E0","#7997FF","#00C0B8","#39B600","#BB9D00","#F37B59","#FF6C90","#BF80FF","#00B4EF","#00BF7D","#85AD00","#E08B00")[1:(length(unique(xyvl_framenum$MG_cols))-2)], '#696969'), na.value = "grey80")
return(pre_anim_gg)
}
