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
