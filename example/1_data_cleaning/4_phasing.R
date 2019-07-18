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