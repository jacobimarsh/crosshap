set.seed(153)

umap_in <- umap::umap(crosshap::LD, min_dist = 2, spread = 2.5, n_neighbors = 30)
test_anim_gg2 <- crosshap::prepare_hap_umap(umap_in,
                                            HapObject = crosshap::Haplotypes_MGmin30_E0.6,
                                            vcf = crosshap::vcf,
                                            nsamples = 2)

test_that("test umap composite plot NAs", {
expect_equal(sum(is.na(test_anim_gg2$data)), 87427)
})
