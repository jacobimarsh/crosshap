#' Hap/MG summary tables
#'
#' build_summary_tables() builds summary tables for each haplotype and Marker
#' Group with some of the information shown in the peripheral crosshap plots.
#' It is an internal function called by crosshap_viz(), though can be called
#' separately to build stand-along grob tables.
#'
#' @param HapObject Haplotype object created by run_haplotyping().
#'
#' @export
#'
#' @return A list containing two TableGrob objects.
#'

build_summary_tables <- function(HapObject){

#Filter out unassigned individuals
no0Varfile <- HapObject$Varfile %>% dplyr::filter(MGs != 0)

#Format MG data in clean tibble to be build as tablegrob
MGdata <- dplyr::left_join(
no0Varfile %>% dplyr::count(MGs) %>% dplyr::rename(nSNP = 'n'),
stats::aggregate(no0Varfile$phenodiff,
                 base::list(no0Varfile$MGs),
                 mean) %>% dplyr::rename(MGs = 'Group.1', phenodiff = 'x') %>%
  tibble::as_tibble(),
by = "MGs") %>%
  dplyr::left_join(

stats::aggregate(no0Varfile$meanr2,
                 base::list(no0Varfile$MGs),
                 mean) %>% dplyr::rename(MGs = 'Group.1', meanR2 = 'x') %>%
  tibble::as_tibble(),
by = "MGs") %>%
  dplyr::left_join(
    stats::aggregate(no0Varfile$AltAF,
                       base::list(no0Varfile$MGs),
                       mean) %>% dplyr::rename(MGs = 'Group.1', AltAF = 'x') %>%
      tibble::as_tibble(),
    by = "MGs") %>%
  dplyr::mutate_if(is.double, function(x){round(x, digits = 2)})

#Build basic tableGrob using MGdata with alternating shaded rows
basic_MGgrob <- gridExtra::tableGrob(MGdata %>% tibble::column_to_rownames('MGs'),
                          theme = ggpp::ttheme_gtstripes(
  colhead = list(bg_params = list(fill = "white"),
                 fg_params = list(fontface = 2L)),
  rowhead = list(bg_params = list(fill = "white"),
                 fg_params = list(fontface = 2L))
))

#Add a line at the bottom of the grob
MG_botline <- gtable::gtable_add_grob(basic_MGgrob,
                             grobs = grid::segmentsGrob(
                               x0 = grid::unit(0,"npc"),
                               y0 = grid::unit(0,"npc"),
                               x1 = grid::unit(1,"npc"),
                               y1 = grid::unit(0,"npc"),
                               gp = grid::gpar(lwd = 1)),
                             t = nrow(basic_MGgrob), b = nrow(basic_MGgrob), l = 2, r = ncol(basic_MGgrob))

#Add a line  under column names
MG_colnamesline <- gtable::gtable_add_grob(MG_botline,
                        grobs = grid::segmentsGrob(
                          x0 = grid::unit(0,"npc"),
                          y0 = grid::unit(1,"npc"),
                          x1 = grid::unit(1,"npc"),
                          y1 = grid::unit(1,"npc"),
                          gp = grid::gpar(lwd = 1)),
                        t = 2, b = 2, l = 2, r = ncol(MG_botline))

#Add a line at the top above column names
MG_final <- gtable::gtable_add_grob(MG_colnamesline,
                        grobs = grid::segmentsGrob( # line across the bottom
                          x0 = grid::unit(0,"npc"),
                          y0 = grid::unit(1,"npc"),
                          x1 = grid::unit(1,"npc"),
                          y1 = grid::unit(1,"npc"),
                          gp = grid::gpar(lwd = 1)),
                        t = 1, b = 1, l = 2, r = ncol(MG_colnamesline))

#MGtable <- ggplot2::ggplot() + MG_final + ggplot2::theme_minimal()

#The next few lines progressively organise and build the data for the hap table
#First, calculate phenotype averages for each haplotype
hap_pheno <- HapObject$Indfile %>%
  dplyr::filter(hap != 0) %>%
  dplyr::group_by(hap) %>%
  dplyr::summarise(phenav = mean_na.rm(Pheno)) %>%
  tidyr::spread(hap, phenav) %>%
  tibble::as_tibble() %>%
  dplyr::mutate_if(is.double, function(x){signif(x, 3)}) %>%
  dplyr::mutate_if(is.double, as.character) %>%
  dplyr::mutate(rname = "Pheno") %>%
  tibble::column_to_rownames("rname")

#Build a table summarising metadata frequency across haplotypes
temp_meta <- suppressMessages(HapObject$Indfile %>%
                                   dplyr::group_by(hap, Metadata) %>%
                                   dplyr::summarise(counts = length(Metadata)) %>%
                                   dplyr::filter(hap != 0) %>%
                              tidyr::spread('hap', 'counts'))
temp_meta$Metadata[is.na(temp_meta$Metadata)] <- "NA"
temp_meta[is.na(temp_meta)] <- 0
hap_meta <- tibble::column_to_rownames(temp_meta, "Metadata") %>% as.matrix()

#Extract total frequency of each haplotype
hap_total <- HapObject$Hapfile %>%
  dplyr::select(hap, n) %>%
  tidyr::spread(hap, n) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(rname = "nTotal") %>%
  tibble::column_to_rownames("rname")

#Glue together
hapdata <- rbind(hap_pheno, hap_meta, hap_total)

#Don't need hap_meta when metadata isn't present
nometa_data <- rbind(phen, n)

#Ensures table has proper formatting without metadata
basic_hapgrob <- gridExtra::tableGrob(if(nrow(hap_meta) == 1){nometa_data}else{hapdata},
                          theme = ggpp::ttheme_gtstripes(
                            colhead = list(bg_params = list(fill = "white"),
                                           fg_params = list(fontface = 2L)),
                            rowhead = list(bg_params = list(fill = "white"),
                                           fg_params = list(fontface = 2L))))


hap_botline <- gtable::gtable_add_grob(basic_hapgrob,
                                grobs = grid::segmentsGrob(
                                  x0 = grid::unit(0,"npc"),
                                  y0 = grid::unit(0,"npc"),
                                  x1 = grid::unit(1,"npc"),
                                  y1 = grid::unit(0,"npc"),
                                  gp = grid::gpar(lwd = 1)),
                                t = nrow(basic_hapgrob), b = nrow(basic_hapgrob), l = 2, r = ncol(basic_hapgrob))

hap_colnamesline <- gtable::gtable_add_grob(hap_botline,
                                grobs = grid::segmentsGrob(
                                  x0 = grid::unit(0,"npc"),
                                  y0 = grid::unit(1,"npc"),
                                  x1 = grid::unit(1,"npc"),
                                  y1 = grid::unit(1,"npc"),
                                  gp = grid::gpar(lwd = 1)),
                                t = 2, b = 2, l = 2, r = ncol(hap_botline))

hap_final <- gtable::gtable_add_grob(hap_colnamesline,
                                grobs = grid::segmentsGrob(
                                  x0 = grid::unit(0,"npc"),
                                  y0 = grid::unit(1,"npc"),
                                  x1 = grid::unit(1,"npc"),
                                  y1 = grid::unit(1,"npc"),
                                  gp = grid::gpar(lwd = 1)),
                                t = 1, b = 1, l = 2, r = ncol(hap_colnamesline))
return(list(MG_final, hap_final))
}
