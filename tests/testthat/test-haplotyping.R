crosshap::run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = 0.58)

test_that("multiplication works", {
  expect_equal(ncol(Haplotypes_MGmin30_E0.58$Hapfile), 11)
})
