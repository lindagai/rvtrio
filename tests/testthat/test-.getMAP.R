context("MAP dimensions and columns are correct")

test_that("Dimensions are correct", {
        #Cols: 1. Gene 2. VarID 3. MAF
        #Rows: Var under analysis
        expect_equal(nrow(map),nrow(vcf.geno))
        expect_equal(ncol(map),3)
})

#TODO: should MAP elements be numeric or char? may need to remove as.vector
test_that("Entries are correct type", {
        expect_equal(unique(as.vector(map[,1])),"NA")
        expect_true(all.equal(as.vector(map[,2]),rownames(vcf.geno)))
        expect_true(is.numeric(as.numeric(map[,3]))) #MAF is numeric
        expect_false(anyNA(map)) #Check for NAs
})
