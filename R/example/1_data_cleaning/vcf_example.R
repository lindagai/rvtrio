
#TODO: Automate: BEAGLE phasing

################################################################################

devtools::install_github("lindagai/rvtrio")
#TODO: Fix importFrom in RV_TDT
library(VariantAnnotation)
library(trio)
library(rvtrio)
library(tidyverse)

################################################################################
#0. Set up

#i. Set working directory
setwd("/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio_example")

################################################################################
#I. Load input files
filepath.functions<-"filepath"
source(filepath.functions)

filepath.vcf<-"/dcl01/beaty/data/gmkf/euro/vcfs/filtered/8q24.recode.vcf"
filepath.ped<-"/dcl01/beaty/data/gmkf/euro/peds/gmkf_euro_completetrios.csv"
hg.assembly<-"hg19"

vcf <- readVcf(filepath_vcf, hg.assembly)
ped <- read.csv(filepath_ped)

################################################################################
#II. Clean VCF file

#Examine the genotypes of the VCF
table(geno(vcf)$GT)

#i. Remove missing and half-calls
geno(vcf)$GT[geno(vcf)$GT == "."]<-"./."
geno(vcf)$GT[geno(vcf)$GT == "./1"]<-"./."
geno(vcf)$GT[geno(vcf)$GT == "./0"]<-"./."
geno(vcf)$GT[geno(vcf)$GT == "1/."]<-"./."
geno(vcf)$GT[geno(vcf)$GT == "0/."]<-"./."

#ii. Genotype formatting
#Ensure that the VCF's genotypes are in 0/0 format so BEAGLE can read it

#iii. Remove multi-allelic SNPS/duplicate sites
duplicate.sites<-start(rowRanges(vcf))[duplicated(start(rowRanges(vcf)))])),
vcf<-vcf[-(which(start(rowRanges(vcf)) %in% duplicate.sites,]

#iii. Remove individuals with Mendelian errors
###################################### TODO: THIS SHOULD BE ITS OWN FILE
mendelian.errors<-.getMendelianErrors() #Include this in rvtrio?
sort(mendelian.errors)

#v. Write out
filepath.filtered.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.vcf"
writeVcf(vcf,filepath.filtered.vcf)

################################################################################
#III. Clean PED file for VCF

#i. Ensure column names are correct and select only the necessary columns
ped<-ped %>%
  rename("famid"="Family.ID",
         "pid"="Individual.ID",
         "fatid"="Father.ID",
         "motid"="Mother.ID",
         "sex" = "Gender",
         "affected"="Clinical.Status")
  select( "famid", "pid", "fatid", "motid", "sex","affected")

head(ped)

#ii. Ensure the PIDs in the VCF and PED files match

#Examine PIDs in VCF and PED
vcf.pid<-colnames(geno(vcf)$GQ)
ped.pid<-ped$pid
head(cbind(vcf.pid,ped.pid))

#Modify PIDs in PED to match VCF
ped <- ped %>% 
        mutate(pid = paste0(famid,"-",pid)) %>%
        mutate(fatid = ifelse(fatid == 0, fatid,
                              paste0(famid,"-",fatid))) %>%
        mutate(motid = ifelse(motid == 0, motid,
                              paste0(famid,"-",motid)))
                              
#Examine PIDs in VCF and PED

#iii. Ensure sex, affected are coded as 0/1, not "male/female" or "unaffected/affected"
ped <- ped %>%
        mutate(sex = ifelse(affected == 2, 1, 0)) %>%
        mutate(affected = ifelse(affected == 2, 1, 0))
	
#iv. Remove individuals with Mendelian errors
###################################### TODO: THIS SHOULD BE ITS OWN FILE
checks<-trio::trio.check(dat, is.linkage=)
trio.check$errors

#vii. Write out
filepath.ped<-"/users/lgai/8q24_project/data/processed_data/gmkf_euro_completetrios_ids_match_vcf.txt"
write.table(chr8.ped, filepath_ped, sep=" ", col.names = TRUE, row.names = FALSE,quote = FALSE)

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