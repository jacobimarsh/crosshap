test_epsilon <- c(0.61)

set.seed(153)

crosshap::run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = test_epsilon)

test_that("test crosshap viz", {
  haplotype_viz2 <- crosshap_viz(Haplotypes_MGmin30_E0.61)
  vdiffr::expect_doppelganger("haplotype_viz2", haplotype_viz2)
})

test_that("test alt crosshap viz", {
alt_viz2 <- crosshap_viz(Haplotypes_MGmin30_E0.61, plot_left = 'pos', plot_right = 'cluster')
vdiffr::expect_doppelganger("haplotype_viz_alt2", alt_viz2)
})
