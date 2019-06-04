library(trio)
library(dplyr)
data(trio.data)

#Run example analysis

#TODO: example data in trio does not have unique identifier for pids/motid/fatid
#You need to include this info in the walkthrough

#TODO: Run through all your RV-TDT code with the example ped
# file to make sure it all works
#You can create a pedfile with MAF <0.05 to filter to afterwards.
#An analyst will probably have to filter their vcf using VariantAnnotation.
#They can use read.plink.

#For PLINK files

ped <- trio.ped1 %>%
        mutate(pid = paste0(famid,"-",pid)) %>%
        mutate(fatid = ifelse(fatid == 0, fatid,
                              paste0(famid,"-",fatid))) %>%
        mutate(motid = ifelse(motid == 0, motid,
                              paste0(famid,"-",motid)))

#This should be abstracted
geno<-ped2geno(ped)
ped <-trio.ped2[,1:6]

#TODO: .getMAP cannot handle the case where SNPs are missing.
#So the user must specift whether to include them

tped<-.getTPED(geno)
map<-.getMAP(tped)
map

testdf<-.runRV_TDT(ped, map, tped)
testdf



window.size=25
window.type="M"
i<-1
n.snps

RV_TDT(geno, ped, window.size=25)

#################### For VCF files ###################

ped<-filepath

#OR, can be extracted from vcf somehow
#figure this out later -- not used in your data

geno<-vcf2geno()


################### Actual scratch #####################################

mafs<-rowMeans(tped)
range(mafs)

table(ped$fatid)
table(ped$motid)

#TODO:
geno<-ped2geno(trio.ped1)
geno<-ped2geno(ped)

#matrix
#Columns: each SNP representing the genotypes of the respective SNP

# genotypes are coded by 0, 1, 2 (i.e. the number of minor alleles)
# matrix contains 3 ∗ t values for

#each SNP genotyped at the t trios, where each block of 3 values is
#composed of the genotypes of the
#father
#mother
#offspring
#of a specific trio.

geno[1:5,1:5]

#RV TDT tped columns:
# 1) SNP/variant id
# 2-n) genotype on every individual (for the rest of the columns)

#vcf.geno.bool<-apply(vcf.geno, 1, function(x) paste(gsub("\\/", "\t", x), collapse="\t"))
#toWrite<-paste(rownames(vcf.geno), vcf.geno.bool, sep=" ")
#writeLines(toWrite, filepath.tped)

tped<-t(geno)
#rownames are SNPs, colnames are pids

geno[1:5,1:5]
tped[1:5,1:5]

pid2<-mutate(pid2 = paste0(famid,"-",pid))

colMAFtrio(geno, changeMinor = FALSE)
