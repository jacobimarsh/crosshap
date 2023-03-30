crosshap::run_haplotyping(vcf = crosshap::vcf,
                          LD = crosshap::LD,
                          pheno = crosshap::pheno,
                          metadata = crosshap::metadata,
                          epsilon = 0.59)

test_that("Marker Groups found in HapObject", {
  expect_equal(length(colnames(Haplotypes_MGmin30_E0.59$Hapfile)), 11)
})

test_that("number of NAs in each bucket", {
  expect_equal(sum(is.na(Haplotypes_MGmin30_E0.59$Varfile)), 252)
  expect_equal(sum(is.na(Haplotypes_MGmin30_E0.59$Indfile)), 54)
})

test_that("wrong vcf input doesn't return haplotype object", {
  expect_error(crosshap::run_haplotyping(vcf = crosshap::metadata,
                                         LD = crosshap::LD,
                                         pheno = crosshap::pheno,
                                         metadata = crosshap::metadata,
                                         epsilon = 0.79))
})
