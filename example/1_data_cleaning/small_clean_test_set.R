#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu
ssh -X lgai@jhpce01.jhsph.edu -o ForwardX11Timeout=336h

#Ensure X11 forwarding is set up so you can graph

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#Create small test set

filepath.mend.err.rm.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.formatted.mend.err.rm.07_15_19.vcf"
hg.assembly<-"hg19"
vcf<-VariantAnnotation::readVcf(filepath.mend.err.rm.vcf, hg.assembly)

#Create small test set
sm.vcf<-vcf[1:10,]

filepath.sm.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.formatted.mend.err.rm.07_15_19.10snps.vcf"
VariantAnnotation::writeVcf(sm.vcf,filepath.sm.vcf)

filepath.cluster<-"/users/lgai/8q24_project/data/processed_data/vcfs/8q24.cleaned.phased.formatted.mend.err.rm.07_15_19.10snps.vcf"
filepath.destination<-" '/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/8q24.cleaned.phased.formatted.mend.err.rm.07_15_19.10snps.vcf'"

#rename this 8q24.10snps.vcf

scp.command<-paste0("scp lgai@jhpce01.jhsph.edu:", filepath.cluster, " ", filepath.destination)
scp.command
system(scp.command)

