context("TPED dimensions and columns are correct")

#TODO: Fix these
test_that("TPED's column names (PIDs) are correct", {
        expected.snps <-colnames(vcf.geno)
        	#Remove the "_1" or "_2" (indicates which allele) from TPED PIDs
        tped.snps <- unique(substr(colnames(tped), 1, nchar(colnames(tped))-2))
        expect_equal(tped.snps, expected.snps)
})

test_that("TPED and PED have PIDs that are in the same order", {
        tped.pids <-rownames(tped)
        expect_true(identical(tped.pids, rownames(vcf.geno)))
})

test_that("Dimensions are correct", {
        #nrows = no. of SNPs
        #ncols = no. of alleles (2*no. of individuals)
        expect_equal(nrow(tped), nrow(vcf.geno))
        expect_equal(ncol(tped), ncol(vcf.geno)*2)
})

test_that("Entries are either 0 or 1", {
        expect_true(all.equal(unique(as.vector(as.matrix(tped))), c(0,1)))
})

test_that("TPED's row names (SNPs) are correct", {
		tped.snps <- rownames(tped)
        expect_true(identical(rownames(tped), rownames(vcf.geno)))
})
