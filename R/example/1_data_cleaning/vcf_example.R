#TODO: Automate: BEAGLE phasing

################################################################################
#0. Set up
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

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

#i. Set half-calls to missing
geno(vcf)$GT[geno(vcf)$GT == "."]<-"."
geno(vcf)$GT[geno(vcf)$GT == "1/."]<-"."
geno(vcf)$GT[geno(vcf)$GT == "0/."]<-"."

#ii. Genotype formatting
#Ensure that the VCF's genotypes are in 0/0 format so BEAGLE can read it

#iii. Remove multi-allelic SNPs/duplicate sites
duplicate.sites<-start(rowRanges(vcf))[duplicated(start(rowRanges(vcf)))]
vcf<-vcf[-which(start(rowRanges(vcf)) %in% duplicate.sites),]

#iV. Remove unnamed SNPs
names<-names(rowRanges(vcf))
names
sum(duplicated(names(rowRanges(vcf))))
duplicate.snps<-names(rowRanges(vcf))[duplicated(names(rowRanges(vcf)))]
duplicate.snps
vcf<-vcf[-which(start(rowRanges(vcf)) %in% duplicate.sites),]

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

#Identify any PIDs in VCF but not in PED
ped.pid<-ped$pid

#Identify any PIDs in VCF but not in PED
setdiff(vcf.pid,ped.pid)
#Identify any PIDs in PED but not in VCF
setdiff(ped.pid,vcf.pid)

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
filepath.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.txt"
filepath.phased.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased"

phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)
phase.command

system(phase.command)

################################################################################

#Now you are ready to start analyzing your data! Filtering by annotation is recommended for rare variants analysis.

#test
filepath.vcf.test<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
vcf <- readVcf(filepath.vcf.test, hg.assembly)
table(geno(vcf)$GT)
#       .      0/.      0/0      0/1      1/.      1/1

filepath.vcf.phased<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.vcf"
vcf <- readVcf(filepath.vcf.phased, hg.assembly)
table(geno(vcf)$GT)
#     0|0      0|1      1|0      1|1
# 13906325   339580   331028   200851


