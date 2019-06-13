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

#TODO: Make sure all filepaths would work on both Mac OS X and Windows

########################################################

RV_TDT<-function(plink.ped=NULL, vcf = NULL, vcf.ped = NULL, window.size=0, window.type = "M", param){

        #TODO: Add some error-checking

        #Check to see if RV_TDT has been installed
        #TODO: Is there any way to install other software automatically?
        #Install RV_TDT if it has not already been installed
        # if (!.checkRV_TDT()){
        # 	.installRV_TDT()
        # }

        #Obtain ped/tped and df of snp.names/positions based on file input type
        if(!is.null(plink.ped)){
                ped <-plink.ped[,1:6]
                tped<-.getTPED(plink = plink.ped)

        } else if (!is.null(vcf) & !is.null(vcf.ped)){
                geno<-as.data.frame(geno(vcf)$GT)
                snps<-as.data.frame(names(geno))
                names(snps)<-"pids"
                ped<-left_join(ped,snps, by=c("pid" ="pids"))
                tped<-.getTPED(vcf.geno = geno)

        } else {
                print("Check your input files.")
                return(NULL)
        }

        map<-.getMAP(tped)
        results<-.runRV_TDT(ped, map, tped)
        return(results)
     
}

########################################################

# Helper Functions - RV_TDT

########################################################

.checkRV_TDT<-function(){
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

.getTPED<-function(plink = NULL, vcf.geno = NULL){

        if(!is.null(plink)){
                tped<-as.data.frame(t(plink[,7:ncol(plink)]))
        } else if (!is.null(vcf.geno)) {
                tped <- as.data.frame(unlist
                                      (lapply(vcf.geno,data.table::tstrsplit, "/"),
                                              recursive = FALSE)
                )
                tped<-sapply(tped, unfactor)
                rownames(tped)<-rownames(vcf.geno)
        }
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
        mafs<-as.vector(rowMeans(tped))
        map<-as.data.frame(cbind(gene.id.vec,snps,mafs))
        return(map)

}

############################

.runRV_TDT<-function(ped, map, tped, window.size=0,window.type  = "M"){
	
	    n.snps<- nrow(map)

        #return rv_tdt results as df
        if (window.size==0){
                n.windows<-1
                window.size<-n.snps
        } else {
                n.windows<- max(n.snps - window.size + 1,1)
        }

        results.df<-as.data.frame(matrix(data=NA,nrow=n.windows, ncol=10))
        colnames(results.df)<-c(
        		"gene.name",
                "CMC.Analytical","BRV.Haplo","CMC.Haplo","VT.BRV.Haplo","VT.CMC.Haplo","WSS.Haplo",
                "start.pos","mid.window.pos","end.pos"
        )

        #TODO: Get rid of for loop and use a list instead
        #Make this into a function, and use apply instead?
        #Indexing is very slow for data frames

        for (i in (1:n.windows)){
        	    start.index<-i
        		end.index<-min(i+window.size-1,n.snps)
        		mid.index<- floor(start.index + (end.index - start.index)/2)
                input.filepaths<-.getInputFilesForWindow(ped,map,tped,window.type,start.index,end.index,i)
                curr.window.result<-.runRV_TDTOnWindow(input.filepaths,adapt=100,alpha=0.00001,permut=1000,u=0.2)
                results.df[i,1:7]<-curr.window.result
                #pos.info<-.getWindowPos(i,start.index,mid.index,end.index)
                results.df[i,1:7]<-curr.window.result
          
        }

        return(results.df)

}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow<-function(ped, map, tped, window.type,start.index,end.index,i){
        
        #TODO: fix this to get the directory functions.R is in
        data.dir<-"./R/"
        
        file.param<-paste0("window",i,".",
                           start.index,"-",end.index,window.type)
        file.param
        
        filepath.map.new<-paste0(data.dir,file.param,".map")
        filepath.tped.new<-paste0(data.dir,file.param,".tped")
        filepath.tped.new
        filepath.ped<-paste0(data.dir,"pedfile.ped")
        filepath.ped
        
        sm.map<-map[start.index:end.index,]
        sm.tped<-tped[start.index:end.index,]
        
        write.table(sm.map,filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = TRUE,quote = FALSE)
        
        if (i==1){
        		#Make sure ped file sample IDs are the same order as tped
        		head(ped)
                write.table(ped,filepath.ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        }
        
        filepaths<-c(filepath.tped.new, filepath.ped, filepath.map.new)
}

############################

.runRV_TDTOnWindow<-function(input.filepaths,adapt=100,alpha=0.00001,permut=1000,u=0.2){

		#TODO: Fix u to be 0.01
		#TODO: Add in all the parameters that RV-TDT includes
		#TODO: Modify command accordingly

        .calculateRV_TDTOnWindow(input.filepaths)
        results<-.extract.results()
        .clean.up.rv_tdt(input.filepaths)
        return(results)

}

############################

.calculateRV_TDTOnWindow<-function(input.filepaths,adapt=100,alpha=0.00001,permut=1000,u=0.2){
	
		#TODO: Fix u to be 0.01
		#TODO: Add in all the parameters that RV-TDT includes
		#TODO: Modify command accordingly

        #filepath.tped<-paste0("'",input.filepaths[1],"'")
        #filepath.phen<-paste0("'",input.filepaths[2],"'")
        #filepath.map<-paste0("'",input.filepaths[3],"'")
        
        filepath.tped<-paste0(input.filepaths[1])
        filepath.phen<-paste0(input.filepaths[2])
        filepath.map<-paste0(input.filepaths[3])

        #TODO: Fix this to match the filepath in the install RV_TDT directories
        rv.tdt.dir<-"'/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/rv-tdt-master/rvTDT'"
        gene.name<-"NA"

        #TODO: Add directories for each window?
        rv.tdt.results.dir<-paste0(" ./R/",gene.name)
        rv.tdt.results.dir

        command<-paste0(rv.tdt.dir, " ", rv.tdt.results.dir,
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

############################

.getWindowPos<-function(i,start.index,mid.index,end.index){
	
	# Get start, end, middle position of each window 
	# Add to results df for easier graphing
	
}

########################################################

# Helper Functions - .calculateRV_TDTOnWindow

########################################################

.extract.results<-function(){
        gene.name<-"NA"
        filepath.results<-paste0("./R/",gene.name,"_pval/",gene.name,".pval")
        filepath.results
        pval.df<-read.table(filepath.results,comment.char = "",header=TRUE)
        pval.df
        return(pval.df)

}

############################

.clean.up.rv_tdt<-function(input.filepaths){
	
	    #TODO: Fix getwd() to be the R package directory
	    # Use filepath and not paste0 so it works on Windows

        #delete all input files
        filepath.tped<-paste0("'",input.filepaths[1],"'")
        filepath.map<-paste0("'",input.filepaths[3],"'")

        delete.input.files<-paste("rm",
                                  filepath.tped,
                                  filepath.map
        )

        print(delete.input.files)
        #system(delete.input.files)

        #delete all output files
        currwd<-getwd()
        filepath.pval<-paste0("'",currwd,"/R/NA_pval/NA.pval","'")
        filepath.results<-paste0("'",currwd,"/R/NA_rvTDT/NA.rvTDT","'")        
        delete.output.files<-paste("rm", filepath.pval, filepath.results)
        print(delete.output.files)
        #system(delete.output.files)
        
        #delete all output directories
        dir.pval<-paste0("'",currwd,"/R/NA_pval","'")
        dir.results<-paste0("'",currwd,"/R/NA_rvTDT","'")        
        delete.output.dir<-paste("rmdir", dir.pval, dir.results)
        print(delete.output.dir)
        #system(delete.output.dir)

}