
context("Inputs are correct")

test_that("Parameters are numeric or NA", {
        expect_equal(str_length("a"), 1)
        expect_equal(str_length("ab"), 2)
        expect_equal(str_length("abc"), 3)
})

########################################################

test_that("PID in PED are in the same order as in GENO/PED", {

})

########################################################

context("PED and VCF IDs match")

########################################################

context("PED dimensions and columns are correct")

########################################################

context("TPED dimensions and columns are correct")

test_that("Dimensions are correct", {
	#Rows: #number of SNPs
	#Cols: SNP name, then 2*number of individuals
		vcf.geno
		tped<-.getTPED(vcf.geno=)
        expect_equal(nrow(tped),nrow(vcf.geno))
        expect_equal(ncol(tped),3)
})

test_that("Entries are correct", {
		vcf.geno
		tped<-.getTPED(vcf.geno=)
        expect_equal(str_length("abc"), 3)
        expect_equal(ncol(tped),3)
})

########################################################

context("PIDs in PED file are in the same order as in GENO/TPED")

########################################################

context("MAP dimensions and columns are correct")

########################################################

context("Number of windows are correct")

########################################################

context("RV-TDT input files for windows do not contain overlapping SNPs")

########################################################

context("RV-TDT input files for windows have correct dimension and columns")

########################################################

context("All files used to run RV-TDT are deleted after run")

########################################################