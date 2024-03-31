test_epsilon <- 0.62

set.seed(155)

testHapObject <- run_haplotyping(vcf = crosshap::vcf,
                                 LD = crosshap::LD,
                                 pheno = crosshap::pheno,
                                 metadata = crosshap::metadata,
                                 epsilon = test_epsilon)

test_that("test crosshap viz", {
  haplotype_viz4 <- crosshap_viz(HapObject = testHapObject, epsilon = 0.62)
  vdiffr::expect_doppelganger("haplotype_viz4data", haplotype_viz4$data)
})

test_that("test alt crosshap viz", {
  alt_viz4 <- crosshap_viz(HapObject = testHapObject, epsilon = test_epsilon, plot_left = 'pos', plot_right = 'cluster')
  vdiffr::expect_doppelganger("haplotype_viz_alt4data", alt_viz4$data)
})

test_that("test no labels crosshap viz", {
  nolabs3 <- crosshap_viz(HapObject = testHapObject, epsilon = test_epsilon, hide_labels = TRUE)
  vdiffr::expect_doppelganger("haplotype_viz_nolabs3data", nolabs3$data)
})

test_that("test isolate_groups", {
  isolate_wt3 <- crosshap_viz(HapObject = testHapObject, epsilon = test_epsilon, isolate_group = "wt")
  vdiffr::expect_doppelganger("haplotype_viz_isolatewt3$data", isolate_wt3$data)
})
