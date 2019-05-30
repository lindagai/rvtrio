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

RV_TDT<-function(filepath.vcf, filepath.ped, window.size=25, window.type = "M", param){

        vcf<-readVcf(filepath.vcf)
        geno<-VariantAnnotation::geno(vcf)
        rm(vcf)

        ped<-read.table(filepath.ped)
        tped<-.getTPED(geno)
        map<-.getMAP(tped)

        results<-.runRV_TDT(window.size, window.type, ped, map, tped)

}

########################################################

# Helper Functions - RV_TDT

########################################################

.getTPED<-function(geno){
	#TODO: Fix this
	
	    #RV TDT tped columns:
        # 1) SNP/variant id
        # 2-n) genotype on every individual (for the rest of the columns)
        
        vcf.geno.bool<-apply(vcf.geno, 1, function(x) paste(gsub("\\/", "\t", x), collapse="\t"))
        toWrite<-paste(rownames(vcf.geno), vcf.geno.bool, sep=" ")
        writeLines(toWrite, filepath.tped)
        return(tped)

}

############################

.getMAP<-function(tped){
	
		#TODO: fix this
		gene.id<-"NA"
        n.snps<-nrow(tped)
        gene.id.vec<-rep(gene.id,n.snps)
        snps<-tped[,1]
        mafs<-rowMeans(tped[,-1])
        map<-cbind(gene.id.vec,snps,mafs)
        return(map)
        
}

############################

.runRV_TDT<-function(window.size,window.type, ped, map, tped){

        #return rv_tdt results as df

        if (window.size = NULL){
                n.windows<-1
        } else {
                n.snps<- nrow(map)
                n.windows<- n.snps - window.size + 1
        }

        results.df<-as.data.frame(matrix(nrow=n.windows, n.col = 8))
        colnames(results)<-c(
                "X.gene",
                "CMC.Analytical","BRV.Haplo","CMC.Haplo","VT.BRV.Haplo","VT.CMC.Haplo","WSS.Haplo",
                "mid.window.pos"
        )

        #TODO: Get rid of for loop

        for i in (1:n.windows){
 
                start.index<-(i-1)*window.size-(i-1)*overlap+1
                end.index<-min(i*window.size-(i-1)*overlap,n.snps)
                mid.position<- floor(start.index + (end.index - start.index)/2)
        
                input.filepaths<-.getInputFilesForWindow(map,ped,tped)
                curr.window.result<-.runRV_TDTOnWindow(input.filepaths,param)
                results[i,]<-c(curr.window.result, mid.position)
                
        }

        return(results.df)

}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow<-function(window.size, ped, map, tped, i){
	
		#TODO: fix this to get the directory functions.R is in
		data.dir <-getwd()

        file.param<-paste0("window",i,".",
                           start.index,"-",end.index,window.type)

        filepath.map.new<-paste0(data.dir,file.param,".map")
        filepath.tped.new<-paste0(data.dir,file.param,".tped")
        filepath.tped.new

        sm.map<-map[start.index:end.index,]
        sm.tped<-tped[start.index:end.index,]

        write.table(sm.map,filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

        filepaths<-c(filepath.tped.new, ped, filepath.map.new)

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
