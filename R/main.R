######################

#rare variants analysis

######################

#User-specified information

#filter by annotation information
annotation.information<-

filepath.vcf<-
window.size<-
window.type<-
param<-

############################################

#NOTE: These should give you DFs

vcfs.filepaths<-c()
get.rv_tdt.results(vcfs.filepaths, param)

#TODO:  this could be abstracted
# Avoid for loops :/

for (i in vcfs){
filepath.vcfs<-vcfs[I]
filepath.rv_tdt.results<-
rv_tdt.results<-get.rv_tdt.results(filepath.vcf)
write.table(rv_tdt.results, filepath.rv_tdt.results)
}

############################################

rvTDT.results<-get.rvTDT.results()
filepath.rvTDT.results<-
write.table(rv_tdt.results, filepath.rv_tdt.results)

############################################

scan_trio.results<-get.scan_trio.results()
filepath.rvTDT.results<-
write.table(rv_tdt.results, filepath.rv_tdt.results)
