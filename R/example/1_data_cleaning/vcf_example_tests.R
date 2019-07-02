#Tests
#Are the duplicated sites removed as expected?
vcf.test<-vcf[-which(start(rowRanges(vcf)) %in% duplicate.sites),]

dup.test<-start(rowRanges(vcf.test))[duplicated(start(rowRanges(vcf.test)))]
dup.test

test<-vcf.test[which(start(rowRanges(vcf.test)) %in% duplicate.sites),]
test
