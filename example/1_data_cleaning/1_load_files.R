
################################################################################
#0. Set up
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu
ssh -X lgai@jhpce01.jhsph.edu -o ForwardX11Timeout=336h

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