#Loads files to unit test helper functions

testthat::expect_true(file.exists(file.path(system.file("data", package="rvtrio"), "8q24.cleaned.10snps.vcf")))
sm.vcf.fp<-file.path(system.file("data", package="rvtrio"), "8q24.cleaned.10snps.vcf")
hg.assembly<-"hg19"
sm.vcf<-VariantAnnotation::readVcf(sm.vcf.fp, hg.assembly)

#TODO: remove the triple colons?
vcf.geno<-VariantAnnotation::geno(sm.vcf)$GT
tped<-rvtrio:::.getTPED(vcf.geno=vcf.geno)
map<-rvtrio:::.getMAP(tped)
map
