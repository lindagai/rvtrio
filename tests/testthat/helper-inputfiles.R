################################################################################

#helper-inputfiles.R

################################################################################

#Loads files to unit test helper functions

################################################################################

#1. Load VCF and PED
testthat::expect_true(file.exists(
        file.path(system.file("extdata", "hg38.vcf", package="rvtrio")
                )
        )
)

testthat::expect_true(file.exists(
        file.path(system.file("extdata", "hg38.ped.txt", package="rvtrio")
                )
        )
)


sm.vcf.fp <- file.path(system.file("extdata", "hg38.vcf", package="rvtrio"))
hg.assembly <- "hg38"
vcf <- VariantAnnotation::readVcf(sm.vcf.fp, hg.assembly)

fp.ped <- file.path(system.file("extdata", "hg38.ped.txt", package="rvtrio"))
ped <- read.table(fp.ped, header=TRUE, stringsAsFactors = FALSE)

################################################################################

#2. Load RV-TDT input functions (MAP, TPED, and PED)

vcf.geno <- VariantAnnotation::geno(vcf)$GT
new.ped <- rvtrio:::.getPED(ped, vcf.geno)
tped <- rvtrio:::.getTPED(vcf.geno)
map <- rvtrio:::.getMAP(tped)