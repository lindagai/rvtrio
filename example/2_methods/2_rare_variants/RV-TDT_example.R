
################################################################################

devtools::install_github("lindagai/rvtrio")
#TODO: Fix importFrom in RV_TDT
library(VariantAnnotation)
library(rvtrio)
library(tidyverse)

################################################################################
#0. Set up
setwd("/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio_example")
filepath.to.RV_TDT<-"/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT"

################################################################################
#I. Load input files

filepath.vcf<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"
filepath.vcf.ped<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"

ped<-read.table(filepath.vcf.ped,header=FALSE)
colnames(ped)<-c("pid","famid","fatid","motid","sex","affected")
vcf<-VariantAnnotation::readVcf(filepath.vcf, "hg19")

################################################################################
#II. Run RV-TDT with various inputs

#Run on the entire VCF (i.e. 1 window)
RV_TDT.results<-rvtrio::RV_TDT(vcf=vcf, vcf.ped = ped, rv.tdt.dir = filepath.to.RV_TDT)
RV_TDT.results

#Run RV-TDT on 25-marker windows

RV_TDT.results.360M.windows<-rvtrio::RV_TDT(vcf=vcf, vcf.ped = ped, rv.tdt.dir = filepath.to.RV_TDT,window.size=360,upper_cutoff=0.1)
head(RV_TDT.results.360M.windows)

#ii. Write out RV-TDT results in one .txt file
filepath.results<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio_example/example-vcf/RV_TDT.results.360M.windows.txt"
write.table(RV_TDT.results.360M.windows,filepath.results, sep="\t", row.names=FALSE,quote=FALSE)

################################################################################
#III. Graph RV-TDT results with ggplot2 (IF TESTING MULTIPLE WINDOWS)

filepath.results<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio_example/example-vcf/RV_TDT.results.360M.windows.txt"
RV_TDT.results<-read.table(filepath.results,sep="\t",header=TRUE, quote ="")
head(RV_TDT.results)

n.windows<-nrow(RV_TDT.results)

#Convert to long format
RV_TDT.results.long <- RV_TDT.results %>%
        gather(key = test, value = pval,
        CMC.Analytical,BRV.Haplo,CMC.Haplo,VT.BRV.Haplo,VT.CMC.Haplo,WSS.Haplo)

RV_TDT.results.long %>% head

#Bonferroni corrected
bonferroni.sig.level<-0.05/n.windows

#Plot
ggplot() +
  geom_line(data = RV_TDT.results.long, aes(group=test, color = test,
                               x = mid.window.pos, y = pval))+
    geom_hline(yintercept=bonferroni.sig.level, linetype=2, color = "red", size=2) +
    labs(title='RV-TDT results for window size = 360 SNPs, 0 overlap',
         x ='Position (hg19)', y = 'p-value at center of window')+
    guides(color=guide_legend("RV-TDT test type")) +
  scale_linetype_manual(name = "Bonferroni-corrected significance", values = 2,
                        guide = guide_legend(override.aes = list(color = c("red"))))
