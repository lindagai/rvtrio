context("TPED dimensions and columns are correct")
skip("already tested: 7/16/19")

#TODO: Fix these
# test_that("TPED's columns names (PIDs) are correct", {
#         expected.colnames<-colnames(vcf.geno)
#         expect_equal(colnames(tped)[1],expected.colnames)
# })
#
# test_that("TPED and PED have PIDs that are in the same order", {
#         tped.pids <-rownames(tped) #remove last 2 characters and remove duplicated
#         expect_equal(tped.pids, rownames(vcf.geno))
# })

test_that("Dimensions are correct", {
        #nrows = no. of SNPs
        #ncols = no. of alleles (2*no. of individuals)
        expect_equal(nrow(tped),nrow(vcf.geno))
        expect_equal(ncol(tped),ncol(vcf.geno)*2)
})

test_that("Entries are either 0 or 1", {
        expect_true(all.equal(unique(as.vector(as.matrix(tped))),c(0,1)))
})

test_that("TPED's row names (SNPs) are correct", {
        expect_equal(rownames(tped), rownames(vcf.geno))
})
