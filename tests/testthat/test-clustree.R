test_epsilon <- c(0.39, 0.79, 0.99)

set.seed(154)

HapObject <- run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = test_epsilon)


test_that("test MG clustree", {
MGtree2 <- clustree_viz(HapObject)
vdiffr::expect_doppelganger("MGtreedata", MGtree2$data)
})

test_that("test hap clustree", {
haptree2 <- clustree_viz(HapObject = HapObject, type = 'hap')
vdiffr::expect_doppelganger("haptreedata", haptree2$data)
})
