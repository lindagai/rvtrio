#' Calculating RV-TDT statistic for user-defined windows
#' @param filepath.vcf absolute filepath of vcf
#' @param filepath.ped absolute filepath of vcf
#' @param window.type type of window, either number of markers ("M") or width of kilobase interval ("K")
#' @param window.size window size size of window
#' @param param parameters for RV-TDT

#' @return results data frame containing results from RV-TDT

#' @export
#'
#' @examples RV_TDT(filepath.vcf, filepath.ped)
#' @examples RV_TDT(filepath.vcf, filepath.ped, window.size=10, window.type="M")
#'
#' @export
#'

########################################################

RV_TDT<-function(plink.ped, window.size=25, window.type = "M", param){

	#Check to see if RV_TDT has been installed
	#TODO: Is there any way to install other software automatically?
	#Install RV_TDT if it has not already been installed
	# if (!.checkRV_TDT()){
	# 	.installRV_TDT()
	# }
		ped <-plink.ped[,1:6]
        tped<-.getTPED(plink.ped)
        map<-.getMAP(tped)
        results<-.runRV_TDT(window.size, window.type, ped, map, tped)

}

########################################################

# Helper Functions - RV_TDT

########################################################

.checkRV_TDT()<-function(){
	#TODO: Fix this
	wd<-getwd()
	filepath.RV_TDT<-file.path(wd,"RV-TDT")
	f<-file.path(wd,"RV-TDT")
	file.exists(f)
}

########################################################

.installRV_TDT<-function(){

    # cd "/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/rv-tdt-master"
    # remove rvTDT zip file

	#TODO: Fix this
	wd<-getwd()
	filepath.RV_TDT<-file.path(wd,"rv-tdt-master")

	dl.RV_TDT<-paste0("wget -O ",filepath.RV_TDT,)
	system(dl.RV_TDT)

	make.rvTDT<-"make rvTDT"
	system(make.rvTDT)


}

########################################################

.getTPED<-function(plink.ped){

        tped<-as.data.frame(t(plink.ped[,7:ncol(plink.ped)]))
        return(tped)

}

############################

#1. gene
#2. variant id. The variant id must matches with the variant id in tped file
#3. maf

.getMAP<-function(tped){

        #TODO: fix the gene.id variable and rename it after the window name
        gene.id<-"NA"
        n.snps<-nrow(tped)
        gene.id.vec<-rep(gene.id,n.snps)
        snps<-rownames(tped)
        mafs<-rowMeans(tped)
        map<-cbind(gene.id.vec,snps,mafs)
        return(map)

}

############################

.runRV_TDT<-function(ped, map, tped, window.size=0,window.type  = "M"){

        #return rv_tdt results as df
        if (window.size==0){
                n.windows<-1
                window.size<-nrow(map)
        } else {
                n.windows<- max(n.snps - window.size + 1,1)
        }

        results.df<-as.data.frame(matrix(data=NA,nrow=n.windows, ncol=8))
        colnames(results.df)<-c(
                "X.gene",
                "CMC.Analytical","BRV.Haplo","CMC.Haplo","VT.BRV.Haplo","VT.CMC.Haplo","WSS.Haplo",
                "mid.window.pos"
        )

        #TODO: Get rid of for loop and use a list instead
        #Make this into a function, and use apply instead?
        #Indexing is very slow for data frames

        for (i in (1:n.windows)){
                input.filepaths<-.getInputFilesForWindow(ped, map, tped, window.type, window.size,i)
                curr.window.result<-.calculateRV_TDTOnWindow(input.filepaths)
                results[i,]<-c(curr.window.result, mid.position)
        }

        return(results.df)

}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow<-function(ped, map, tped, window.type, window.size, i){

        n.snps<- nrow(map)

        start.index<-i
        end.index<-min(i+window.size-1,n.snps)
        mid.position<- floor(start.index + (end.index - start.index)/2)

        #TODO: fix this to get the directory functions.R is in
        data.dir <-getwd()

        file.param<-paste0("/R/window",i,".",
                           start.index,"-",end.index,window.type)
        file.param

        filepath.map.new<-paste0(data.dir,file.param,".map")
        filepath.tped.new<-paste0(data.dir,file.param,".tped")
        filepath.ped<-paste0(data.dir,"/R/pedfile.ped")
        filepath.ped

        sm.map<-map[start.index:end.index,]
        sm.tped<-tped[start.index:end.index,]

        write.table(sm.map,filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = TRUE,quote = FALSE)

        if (i==1){
                write.table(ped,filepath.ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        }

        filepaths<-c(filepath.tped.new, filepath.ped, filepath.map.new)

}

############################

.runRV_TDTOnWindow<-function(input.filepaths, param){

        .calculateRV_TDTOnWindow(input.filepaths, param)
        results<-.extract.results(input.filepaths)
        .clean.up.rv_tdt(input.filepaths)
        return(results)

}

############################

.calculateRV_TDTOnWindow<-function(input.filepaths, param=NULL){

        filepath.tped<-paste0("'",input.filepaths[1],"'")
        filepath.phen<-paste0("'",input.filepaths[2],"'")
        filepath.map<-paste0("'",input.filepaths[3],"'")

        #TODO: Add checks for the param variable
        if (is.null(param)){
                adapt<-100
                alpha<-0.00001
                permut<-1000
                u<-0.5 #TODO: Fix this to be 0.01
        } else {
                adapt<-param[1]
                alpha<-param[2]
                permut<-param[3]
                u<-param[4]
        }

        #TODO: Fix this to match the filepath in the install RV_TDT directories
        rv.tdt.dir<-"'/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/rv-tdt-master/rvTDT'"
        gene.name<-"NA"

        #TODO: Add directories for each window?
        rv.tdt.results.dir<-paste0(" ./",gene.name)
        rv.tdt.results.dir

        command<-paste0(rv.tdt.dir, rv.tdt.results.dir,
                        " -G ", filepath.tped,
                        " -P ", filepath.phen,
                        " -M ", filepath.map,
                        " --adapt ", adapt,
                        " --alpha ", alpha,
                        " --permut ", permut,
                        " -u ", u
        )

        print(command)
        system(command)

}

########################################################

# Helper Functions - .calculateRV_TDTOnWindow

########################################################

.extract.results<-function(file.param){

	filepath.results<-paste0(" ./",gene.name,"_pval/",gene.name,".pval")
	pval.df<-read.table(filepath.results,comment.char = "",header=TRUE)
	return(pval.df[1,])

} 

############################

.clean.up.rv_tdt<-function(input.filepaths){
        
        #delete all input files
        filepath.tped<-paste0("'",input.filepaths[1],"'")
        filepath.map<-paste0("'",input.filepaths[3],"'")
        
        delete.input.files<-paste("rm",
                                  filepath.tped,
                                  filepath.map
        )
        
        delete.input.files
        system(delete.input.files)
        
        #delete all output files
        filepath.rv.tdt.results <-paste0()
        command<-paste("rm -rf", filepath.rv.tdt.results)
        command
        system(delete.input.files)
        
}