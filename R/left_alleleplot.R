#' Left SNP-info alleleplot
#'
#' Internal function that plots the mean frequencies of all alleles for loci
#' within each marker group. Makes use of $Varfile information. The first
#' two snippets of code calculates the mean frequencies of each allele type for
#' each marker group, before plotting. May be missing some axis to allow for
#' 'crosshap' stitching.
#'
#' @param HapObject Haplotype object created by crosshap::run_haplotyping
#'
#' @return
#' @export
#'
#' @example
#' build_left_alleleplot(Haplotypes_MP_E2)
#'

build_left_alleleplot <- function(HapObject) {
pre_leftplotdata <- HapObject$Varfile %>% dplyr::filter(cluster > 0) %>%
  select(-avPheno) %>%
  spread(key, nInd) %>%
  rename(ref = '0', het = '1', alt = '2')

leftplot_data <- aggregate(pre_leftplotdata$ref,
          list(pre_leftplotdata$cluster),
          mean) %>% rename(cluster = 'Group.1', ref = 'x') %>%
  left_join(aggregate(pre_leftplotdata$het,
                      list(pre_leftplotdata$cluster),
                      mean) %>% rename(cluster = 'Group.1', het = 'x'), by = "cluster") %>%
  left_join(aggregate(pre_leftplotdata$alt,
                      list(pre_leftplotdata$cluster),
                      mean) %>% rename(cluster = 'Group.1', alt = 'x'), by = "cluster")

left_alleleplot <- ggplot(leftplot_data %>% gather("Type", "nInd", 2:4) %>%
           mutate(Type = factor(Type, levels = c("ref", "het", "alt"))),
         aes(x = nInd,
             y = as.character(cluster),
             fill = Type,
             color = Type),
         ylim()) +
  geom_bar(aes(), position = 'stack', stat = "identity", width = 0.8, color = "black") +
  scale_x_reverse(breaks = scales::pretty_breaks(),
                  expand = c(0,0)) +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_text(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.key.size = unit(5,
                               "mm"),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0.1),
                           "cm"),
        plot.title = element_blank()) +
  scale_fill_manual(values = c("#FFFFFF", "#73D055FF", '#440154FF')) +
  xlab("Allele count") +
  scale_y_discrete(position = "right", limits = rev,
                   labels = c(paste0("MG",as.character(max(HapObject$Varfile$cluster):1))))

return(left_alleleplot)
}
