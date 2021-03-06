context("All files used to run RV-TDT are deleted after run")

filepath.to.RV_TDT <- "/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT"
RV_TDT.results <- RV_TDT(vcf, ped,
                         filepath.RV_TDT = filepath.to.RV_TDT)

test_that("Check whether the input_files directory has been deleted", {
        expect_false(file.exists(
                # system.file("data","input_files", package = "rvtrio")
                file.path(.libPaths(),"rvtrio", "data","input_files")
                                 )
        )
})

test_that("Check whether the results directory has been deleted", {
        expect_false(file.exists(file.path(.libPaths(),"rvtrio", "data","results")
                                 )
        )
})
