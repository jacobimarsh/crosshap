test_epsilon <- c(0.4, 0.8, 1)

set.seed(153)

crosshap::run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = test_epsilon)


test_that("test MG clustree", {
MGtree2 <- clustree_viz(pheno = crosshap::pheno, test_epsilon)
vdiffr::expect_doppelganger("MGtreedata", MGtree2$data)
})

test_that("test hap clustree", {
haptree <- clustree_viz(pheno = crosshap::pheno, test_epsilon, type = 'hap')
vdiffr::expect_doppelganger("haptree", haptree)
})
