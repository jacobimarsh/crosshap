test_epsilon <- 0.61

set.seed(153)

testHapObject <- run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = test_epsilon)

test_that("test crosshap viz", {
  haplotype_viz2 <- crosshap_viz(HapObject = testHapObject, epsilon = 0.61)
  vdiffr::expect_doppelganger("haplotype_viz2", haplotype_viz2)
})

test_that("test alt crosshap viz", {
alt_viz2 <- crosshap_viz(HapObject = testHapObject, epsilon = test_epsilon, plot_left = 'pos', plot_right = 'cluster')
vdiffr::expect_doppelganger("haplotype_viz_alt2", alt_viz2)
})

test_that("test no labels crosshap viz", {
nolabs <- crosshap_viz(HapObject = testHapObject, epsilon = test_epsilon, hide_labels = TRUE)
vdiffr::expect_doppelganger("haplotype_viz_nolabs", nolabs)
})

test_that("test isolate_groups", {
isolate_wt <- crosshap_viz(HapObject = testHapObject, epsilon = test_epsilon, isolate_group = "wt")
vdiffr::expect_doppelganger("haplotype_viz_isolatewt", isolate_wt)
})
