################################################################################

#helper-inputfiles.R

################################################################################

#Loads files to unit test helper functions

################################################################################

#1. Load VCF and PED
testthat::expect_true(file.exists(file.path(system.file("data",
                                                        package="rvtrio"), "hg38.vcf")))

sm.vcf.fp <- file.path(system.file("data", package="rvtrio"), "hg38.vcf")
hg.assembly <- "hg38"
sm.vcf <- VariantAnnotation::readVcf(sm.vcf.fp, hg.assembly)

fp.ped <- file.path(system.file("data", package="rvtrio"), "hg38.ped.txt")
ped <- read.table(fp.ped,header=TRUE) 

################################################################################

#2. Load RV-TDT input functions (MAP, TPED, and PED)

#TODO: remove the triple colons?
vcf.geno <- VariantAnnotation::geno(sm.vcf)$GT
new.ped <- rvtrio:::.getPED(ped, vcf.geno)
tped <- rvtrio:::.getTPED(vcf.geno)
map <- rvtrio:::.getMAP(tped)