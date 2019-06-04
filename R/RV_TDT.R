#' Calculating RV-TDT statistic for user-defined windows, either 
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

RV_TDT<-function(geno, ped, window.size=25, window.type = "M", param){

        tped<-.getTPED(geno)
        map<-.getMAP(tped)
        results<-.runRV_TDT(window.size, window.type, ped, map, tped)

}

########################################################

# Helper Functions - RV_TDT

########################################################

.getTPED<-function(geno){
       
        tped<-t(geno)
        return(tped)

}

############################

#1. gene 
#2. variant id. The variant id must matches with the variant id in tped file
#3. maf

.getMAP<-function(tped){
	
		#TODO: fix this
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
                print(paste0(i, "-",input.filepaths))
                #curr.window.result<-.runRV_TDTOnWindow(input.filepaths,param)
                #results[i,]<-c(curr.window.result, mid.position)
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
        write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

        if (i==1){
                write.table(ped,filepath.ped, sep="\t",col.names=TRUE,row.names = FALSE,quote = FALSE)
        }

        filepaths<-c(filepath.tped.new, filepath.ped, filepath.map.new)

}

############################

.runRV_TDTOnWindow<-function(input.filepaths, param){

        runRV_TDTOnWindow(input.filepaths, param)
        results<-extract.results()
        clean.up.rv_tdt()
        return(results)

}

############################

.runRV_TDTOnWindow<-function(input.filepaths, param=NULL){

        filepath.tped<-input.filepaths[1]
        filepath.phen<-input.filepaths[2]
        filepath.map<-input.filepaths[3]
        
        adapt<-100
        alpha<-0.00001
        permut<-1000
        u<-0.01
        
        #TODO: Fix this
        
        rv.tdt.dir<-"x"
        
        rv.tdt.results.dir<-paste0(data.dir,file.param)
        rv.tdt.results.dir

        command<-paste0(rv.tdt.dir, rv.tdt.results.dir,
                        " -G ", filepath.tped,
                        " -P ", filepath.phen,
                        " -M ", filepath.map,
                        " --adapt 100 --alpha 0.00001 --permut 1000",
                        " -u 0.01"
        )

        command

        system(work.dir)
        system(command)

}

############################

.extract.results<-function(file.param)){

	filepath<-paste0(rv.tdt.results.dir,file.param)
	result<-read.table(filepath)
	return(result)
	
}

############################

.clean.up.rv_tdt<-function(input.filepaths, param){
	
	    filepath.tped<-input.filepaths[1]
        filepath.map<-input.filepaths[3]
        filepath.rv.tdt.results
        

	#delete all extra files
	command<-paste("rm",
	filepath.tped,
	filepath.map,
	)
	
	command
	
	system(command)
	
	command<-paste("rm -rf", filepath.rv.tdt.results)
	command

}
