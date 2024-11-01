set.seed(153)

umap_in <- umap::umap(crosshap::LD, min_dist = 2, spread = 2.5, n_neighbors = 30)
test_anim_gg2 <- prepare_hap_umap(umap_in,
                                  HapObject = crosshap::HapObject,
                                  epsilon = 0.6,
                                  vcf = crosshap::vcf,
                                  nsamples = 2)

test_that("test umap composite plot NAs", {
expect_equal(sum(is.na(test_anim_gg2$data)), 87563)
})
