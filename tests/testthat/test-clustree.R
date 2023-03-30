test_epsilon <- c(0.4, 0.8, 1)

crosshap::run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = test_epsilon)


test_that("test MG clustree", {
MGtree <- clustree_viz(epsilon = test_epsilon, pheno = crosshap::pheno, type = 'MG')
vdiffr::expect_doppelganger("MGtree", MGtree)
})

test_that("test hap clustree", {
haptree <- clustree_viz(epsilon = test_epsilon, pheno = crosshap::pheno, type = 'hap')
vdiffr::expect_doppelganger("haptree", haptree)
})