context("RV_TDT value")

library(testthat)
library(rvtrio)

test_check("rvtrio")

expect_equal(10, 10 + 1e-7)
expect_equal(10, 8)

test_that("RV_TDT returns a data frame", {
        expect_equal(str_length("a"), 1)
        expect_equal(str_length("ab"), 2)
        expect_equal(str_length("abc"), 3)
})

test_that("str_length is number of characters", {
        expect_equal(str_length("a"), 1)
        expect_equal(str_length("ab"), 2)
        expect_equal(str_length("abc"), 3)
})
