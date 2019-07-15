########################################################

# Inputs

########################################################

context("Inputs are correct")

test_that("Parameters are numeric or NA", {

})

########################################################

context("PED and VCF IDs match")

########################################################

context("PED dimensions and columns are correct")

test_that("Dimensions are correct", {
	#of rows in PED = # of rows in VCF
	#colnames are famid, 

})

test_that("Column names are correct", {
	#Rows: #number of SNPs
	#Cols: SNP name, then 2*number of individuals
        expect_equal(nrow(tped),nrow(vcf.geno))
        expect_equal(ncol(tped),3)
})

test_that("Entries are correct", {
	#Rows: #number of SNPs
	#Cols: SNP name, then 2*number of individuals
	
})

########################################################

# Internal functions

########################################################

context("TPED dimensions and columns are correct")

attach(sm.vcf)
vcf.geno<-.geno(sm.vcf)$GT
tped<-.getTPED(vcf.geno=vcf.geno)

test_that("Dimensions are correct", {
	#Rows: #number of SNPs
	#Cols: SNP name, then 2*number of individuals
        expect_equal(nrow(tped),nrow(vcf.geno))
        expect_equal(ncol(tped),ncol(vcf.geno))
})

test_that("Entries are correct type", {
        expect_output(tped, "list")
        expect_output(is.numeric(tped[,-1]), TRUE) #Check whether TPED is only 0/1
           expect_output(is.numeric(tped[,1]), TRUE) #snp.names
             expect_output(anyNA(tped), FALSE) #Check for NAs
})

test_that("TPED/GENO and PED have PIDs that are in the same order", {
        expect_equal(tped[,1]), ped$pid)
})

########################################################

context("MAP dimensions and columns are correct")

attach(sm.vcf)
map<-.getMAP(tped)

test_that("Dimensions are correct", {
	#Cols: 1. Gene 2. VarID 3. MAF
	#Rows: Var under analysis
        expect_equal(nrow(map),nrow(vcf.geno))
        expect_equal(ncol(map),3)
})

test_that("Entries are correct type", {
        expect_output(map, "list")
        expect_output(map[,1])== "NA", TRUE) #Check whether gene names are NA
             expect_output(lengths(map[,2]), TRUE) #Check whether snp.names have > 3 characters
             expect_output(is.numeric(map[,3]), TRUE) #MAF is numeric
             expect_output(anyNA(map), FALSE) #Check for NAs
})

########################################################

context("RV-TDT input files for windows have correct dimension and columns")

test_that("Number of windows are correct", {
        #Check whether RV_TDT has correct dimensions for input files
        	#Check a few different window sizes (0, 10, 300, 1000)
})

test_that("Size of windows are correct", {
        	#Check whether the sizes of the windows for MAP are correct
        	#Check whether the size of the windows for TPED are correct
})

########################################################

context("All files used to run RV-TDT are deleted after run")

test_that("Check", {
        	#Check whether the scratch data directory has been deleted
})

########################################################