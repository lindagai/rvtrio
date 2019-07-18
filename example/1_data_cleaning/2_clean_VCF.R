#II. Clean VCF file

#Examine the genotypes of the VCF
table(geno(vcf)$GT)

#i. Set half-calls to missing; must include "/" separator so BEAGLE can read it
geno(vcf)$GT[geno(vcf)$GT == "."]<-"./."
geno(vcf)$GT[geno(vcf)$GT == "1/."]<-"./."
geno(vcf)$GT[geno(vcf)$GT == "0/."]<-"./."

table(geno(vcf)$GT)

#iii. Remove duplicated sites and multi-allelic SNPs
dups<-duplicated(start(rowRanges(vcf)))
duplicate.sites<-start(rowRanges(vcf))[duplicated(start(rowRanges(vcf)))]
vcf<-vcf[-which(start(rowRanges(vcf)) %in% duplicate.sites),]

#iv. Check for SNPs with duplicate names
names<-names(rowRanges(vcf))
names[1:5]
anyNA(names)
sum(duplicated(names(rowRanges(vcf))))
#0, so we are good

#v. Write out
filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.07_1_19.vcf"
writeVcf(vcf,filepath.filtered.vcf)