#TODO: Automate: BEAGLE phasing

################################################################################
#0. Set up
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

#Ensure X11 forwarding is set up so you can graph

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

devtools::install_github("lindagai/rvtrio")
#TODO: Fix importFrom in RV_TDT
library(VariantAnnotation)
library(trio)
library(rvtrio)
library(dplyr)
library(ggplot2)

################################################################################
#I. Load input files
#filepath.functions<-"filepath"
#source(filepath.functions)

filepath.vcf<-"/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.vcf, hg.assembly)

filepath.ped<-"/dcl01/beaty/data/gmkf/euro/peds/gmkf_euro_completetrios.csv"
ped <- read.csv(filepath.ped)

################################################################################
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

################################################################################
#III. Clean PED file for VCF
head(ped)

#i. Ensure column names are correct
ped<-ped %>%
  rename("famid"="Family.ID",
         "pid"="Individual.ID",
         "fatid"="Father.ID",
         "motid"="Mother.ID",
         "sex" = "Gender",
         "affected"="Clinical.Status")

head(ped)

#ii. Select only the necessary columns
ped<-ped %>%
        select( "famid", "pid", "fatid", "motid", "sex","affected")

head(ped)

#iii. Ensure sex, affected are coded as 0/1, not "male/female" or "unaffected/affected"
ped <- ped %>%
        mutate(sex = ifelse(sex == "male", 1, 0)) %>%
        mutate(affected = ifelse(affected == "Affected", 1, 0))

head(ped)

#iv. Ensure the PIDs in the VCF and PED files match
#All subjects in vcf must also appear (with the same ID) in ped.

#Examine PIDs in VCF and PED
vcf.pid<-colnames(geno(vcf)$GQ)
head(vcf.pid)

ped.pid<-ped$pid
head(ped.pid)

#Modify PIDs in PED to match VCF
ped <- ped %>%
        mutate(pid = paste0("H_TZ-",pid,"-",pid)) %>%
        mutate(fatid = ifelse(fatid == 0, "0",
                              paste0("H_TZ-",fatid,"-",fatid))) %>%
        mutate(motid = ifelse(motid == 0, "0",
                              paste0("H_TZ-",motid,"-",motid)))

#Identify any PIDs in VCF but not in PED and vice versa
ped.pid<-ped$pid
vcf.pid<-colnames(geno(vcf)$GT)
setdiff(vcf.pid,ped$pid)
setdiff(ped$pid,vcf.pid)

#These IDs in VCF contain a B, so we edit them in PED
pids.to.edit <-setdiff(ped.pid,vcf.pid)
pids.to.edit

#Remember to change the fatid and motid as well
head(ped)

ped <- ped %>%
        mutate(pid =  ifelse(pid %in% pids.to.edit,
               paste0(pid,"B"),pid)) %>%
        mutate(fatid =  ifelse(fatid %in% pids.to.edit,
                             paste0(fatid,"B"), fatid)) %>%
        mutate(motid =  ifelse(motid %in% pids.to.edit,
                             paste0(motid,"B"), motid))

head(ped)

#Checks
setdiff(vcf.pid,ped$pid)
setdiff(ped$pid,vcf.pid)
setdiff(ped$fatid,ped$pid)
setdiff(ped$motid,ped$pid)
head(ped)

#vii. Write out
filepath.ped.cleaned<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_07_01_2019.txt"
write.table(ped, filepath.ped.cleaned, sep=" ", col.names = TRUE, row.names = FALSE,quote = FALSE)

################################################################################
#IV. Haplotype phasing

#i. Install BEAGLE software to phase your data

#Un-comment this out if you don't already have BEAGLE4.0
#NOTE: Do not use BEAGLE4.1 or 5.0: these do not use pedigree information.

#Download vcf2beagle
# filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
# dl.beagle4<-paste0("wget -O ",filepath.beagle4," https://faculty.washington.edu/browning/beagle/beagle.r1399.jar")
# system(dl.beagle4)

#Phase vcf
# filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
# filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
# filepath.phased.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased"

filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.07_1_19.vcf"
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_07_01_2019.txt"
filepath.phased.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.07_1_19.phased"

phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)
phase.command

system(phase.command)

#Runtime: 47 minutes 45 seconds

################################################################################
#IV. Checking for Mendelian errors

#Phasing is done first because R's trio.check function does not accept missing observations
#TODO: does it matter whether you phase first?

filepath.phased.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.07_1_19.phased.vcf.gz"
hg.assembly<-"hg19"
vcf <- readVcf(filepath.phased.vcf, hg.assembly)
# test<-geno(vcf)$GT
# test[1:5,1:5]

table(geno(vcf)$GT)
# 0|0      0|1      1|0      1|1
# 13906325   339580   331028   200851

#For trio to read in the genotupes, we must replace all instances of 0|0 with 0/0, etc.
geno(vcf)$GT<-gsub("\\|", "\\/",geno(vcf)$GT)

table(geno(vcf)$GT)
# 0/0      0/1      1/0      1/1
# 14081894   360053   350895   205922

#v. Write out
filepath.formatted.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.07_1_19.phased.formatted.vcf"
writeVcf(vcf,filepath.formatted.vcf)

#The trio package can't read in "1/0", so we must rename it
#i. Set half-calls to missing; must include "/" separator so BEAGLE can read it
geno(vcf)$GT[geno(vcf)$GT == "1/0"]<-"0/1"

#This does not work bc the phased VCF does not include PIDs
trio.geno<-vcf2geno(vcf,ped)
trio.geno[1:5,1:5]

trio.geno

#test<-data.frame(cbind(rownames(trio.geno),trio.geno),c("pid",colnames(trio.geno)))
colnames(trio.geno)[1]<-"pid"
trio.geno[1:5,1:5]

#Check -- do the PIDs match as we expect? PUT THIS INTO TESTS
vcf.pid<-colnames(geno(vcf)$GQ)
setdiff(vcf.pid,ped$pid)
setdiff(ped$pid,vcf.pid)
setdiff(ped$fatid,ped$pid)
setdiff(ped$motid,ped$pid)

library(trio)
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_07_01_2019.txt"
ped <- read.table(filepath.ped,header=TRUE)
head(ped)
ped.test<-ped[,c("famid","pid")]
head(ped.test)

#Check
setdiff(trio.geno$pid,ped.test$pid)
setdiff(ped.test$pid,trio.geno$pid)
head(data.frame(trio.geno$pid,ped.test$pid))
typeof(ped.test$pid)
typeof(trio.geno$pid)

# geno.with.famid<-dplyr::left_join(data.frame(ped.test),data.frame(trio.geno))
# test<-geno.with.famid[1:5,1:5]
# geno.with.famid2<-mapply(geno.with.famid[,3:ncol(geno.with.famid)],as.character)
# geno.with.famid2[1:5]
# geno.with.famid3<-sapply(geno.with.famid2,as.numeric)
# geno.with.famid3[1:5]
# sapply(test2,as.numeric)
#
# typeof(geno.with.famid[3,3])
# dim(geno.with.famid)

#vii. Write out
filepath.geno.with.famid<-"/users/lgai/8q24_project/data/processed_data/geno.with.famid_07_02_2019.txt"
write.table(geno.with.famid, filepath.geno.with.famid, sep=" ", col.names = TRUE, row.names = FALSE,quote = FALSE)

filepath.geno.with.famid<-"/users/lgai/8q24_project/data/processed_data/geno.with.famid_07_02_2019.txt"
geno.with.famid <- read.table(filepath.geno.with.famid,header=TRUE)
geno.with.famid[1:5,1:5]

snp.names<-rownames(geno(vcf)$GT)
snp.names[1:5]

#TODO: Set the geno matrix's genotype entries to numeric, not integer, so you can run trio.check!
#You also need to fix the SNP names to match the VCF
filepath.geno.with.famid<-"/users/lgai/8q24_project/data/processed_data/geno.with.famid_07_02_2019.txt"
geno.with.famid<- read.table(filepath.geno.with.famid,header=TRUE)
geno.with.famid[,3:ncol(geno.with.famid)] <- lapply(geno.with.famid[,3:ncol(geno.with.famid)], as.numeric)
geno.with.famid[1:5,1:5]

#Check
# geno.with.famid[1:5,1:5]
# test<-colnames(geno.with.famid[,3:ncol(geno.with.famid)])
# head(test)
# tail(test)
# head(snp.names)
# tail(snp.names)
#snp.names don't match but it doesn't really matter

#Make sure it's numeric
#geno.with.famid[2:13] <- lapply(geno.with.famid[2:13], as.numeric)

# sm.geno.with.famid<-geno.with.famid[,1:20]
# trio.tmp <- trio::trio.check(dat=sm.geno.with.famid,is.linkage=FALSE)
# str(trio.tmp, max=1)
# test<-trio.tmp$errors
# library(dplyr)
# test %>% filter(famid=="GMKF0097")
# sort(table(trio.tmp$errors$famid),decreasing = TRUE)[1:10]

################################################################################

#Run this on the whole dataset after you're sure it works on the small test set!
trio.tmp <- trio::trio.check(dat=geno.with.famid,is.linkage=FALSE)
#takes awhile to run!

filepath.trio.tmp<-"/users/lgai/8q24_project/data/processed_data/trio.check.output_07_03_2019.txt"
saveRDS(trio.tmp,filepath.trio.tmp)

filepath.trio.tmp.errors<-"/users/lgai/8q24_project/data/processed_data/trio.mend.errors_07_03_2019.txt"
write.table(trio.tmp$errors, filepath.trio.tmp.errors, sep=" ", col.names = TRUE, row.names = FALSE,quote = FALSE)

################################################################################

#Graph the number of Mendelian errors
filepath.trio.tmp.errors<-"/users/lgai/8q24_project/data/processed_data/trio.mend.errors_07_03_2019.txt"
trio.tmp.errors <- read.table(filepath.trio.tmp.errors,header=TRUE)
head(trio.tmp.errors)

mend.err.sorted<-as.data.frame(sort(table(trio.tmp.errors$famid),decreasing = TRUE))
colnames(mend.err.sorted)<-c("famid","mend.errors")
mend.err.sorted[1:20,]

#Save this graph
mend.err.plot.fp<-"/users/lgai/8q24_project/data/processed_data/8q24.mendelian.error.plot.pdf"
pdf(file=mend.err.plot.fp)

ggplot(mend.err.sorted, aes(x=famid,y=mend.errors)) +
        geom_bar(stat="identity") +
        xlab("Family ID") + ylab("Mendelian error count") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank()
              )
dev.off()

################################################################################
#Remove 15 families with large number of

#Remove families with large number of Mendelian errors from VCF and PED

#Remove families with large number of Mendelian errors from VCF

#Remove families with large number of Mendelian errors from PED

################################################################################

#Now you are ready to start analyzing your data! Filtering by annotation is recommended for rare variants analysis.
