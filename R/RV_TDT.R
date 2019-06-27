#' Calculating RV-TDT statistic for user-defined windows

#' `RV_TDT()` returns a data frame containing the RV_TDT statistic for a VCF or PLINK ped file. Note that RV_TDT only works for Linux and Mac OS X.

#' @param plink.ped PLINK ped file
#' @param vcf vcf file
#' @param vcf.ped data frame containing pedigree information for the VCF
#' @param window.type type of window, either number of markers ("M") or width of kilobase interval ("K") (doesn't work yet)
#' @param window.size size of window in number of markers
#' @param param parameters for RV-TDT

#' @param adapt To reduce computational time, adaptive permutation is used in /rvTDT/. Every /$adapt/ permutations (default: 500 permutations), the program will check if we should keep doing permutation (which means this gene looks promising to reach the desired alpha level), or we should give up on this gene (which means this gene will not reach the desired alpha level based on the permutations we have done so far, or we have done enough permutations)

#' @param alpha alpha level in adaptive permutation

#' @param permut The maximum number of permutations

#' @param lower_cutoff The cutoffs to determine which variants we should include in the analysis. In this example, the third column of map file is the number of minor allele counts, and here we only include the variants who have minor allele counts less than 100

#' @param upper_cutoff See lower_cutoff

#' @param minVariants The minimum number of variant sites for a gene. Genes with variant site number less than minVariants will be excluded from analysis (after check missing);

#' @param maxMissRatio The max missing ratio allowed for a variant. The variants with missing ratio greater than maxMissRatio will be excluded from analysis. In this example, we generated the genetic data file without any missing genotypes, so /--maxMissRatio 1/ is used here.

#' @return results data frame containing results from RV-TDT

#' @importFrom VariantAnnotation start geno rowRanges

#' @export
#'
#' @examples RV_TDT(filepath.vcf, filepath.ped)
#' @examples RV_TDT(filepath.vcf, filepath.ped, window.size=10, window.type="M")
#'
#' @export
#'

########################################################

RV_TDT<-function(plink.ped=NULL, vcf = NULL, vcf.ped = NULL, rv.tdt.dir, window.size=0, window.type = "M", adapt = 500, alpha = 0.00001, permut = 2000, lower_cutoff = 0, upper_cutoff = 100, minVariants = 3, maxMissRatio = 1){

	    #Extract parameters
	    parameters<-c(adapt, alpha, permut,lower_cutoff, upper_cutoff, minVariants,maxMissRatio)

        #Obtain ped/tped and df of snp.names/positions based on file input type
        if(!is.null(plink.ped)){
                ped <-plink.ped[,1:6]
                tped<-.getTPED(plink = plink.ped)
                
                #Ensure pids in PED file are in the same order as in GENO/TPED
        } else if (!is.null(vcf) & !is.null(vcf.ped)){
        	    #Get positions and SNP names
                snp.pos.df<-as.data.frame(cbind(names(vcf),start(rowRanges(vcf))))
                colnames(snp.pos.df)<-c("snp.name","pos")
                
                #Get genotypes and SNP names
                geno<-as.data.frame(geno(vcf)$GT)
                snps<-as.data.frame(names(geno))
                
                rm(vcf)

                #Ensure pids in PED file are in the same order as in GENO/TPED
                names(snps)<-"pids"
                ped<-dplyr::left_join(ped,snps, by=c("pid" ="pids"))
                tped<-.getTPED(vcf.geno = geno)
        } else {
                print("Check your input files.")
                return(NULL)
        }

        #Get input files
        map<-.getMAP(tped)

        results<-.runRV_TDT(ped, map, tped, rv.tdt.dir, window.size, snp.pos.df, param = parameters)
        return(results)

}

########################################################

# Helper Functions - RV_TDT

########################################################

.getTPED<-function(plink = NULL, vcf.geno = NULL){

        if(!is.null(plink)){
                tped<-as.data.frame(t(plink[,7:ncol(plink)]))
        } else if (!is.null(vcf.geno)) {
                tped <- as.data.frame(unlist
                                      (lapply(vcf.geno,data.table::tstrsplit, "/"),
                                              recursive = FALSE)
                )
                tped<-sapply(tped, varhandle::unfactor)
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

.runRV_TDT<-function(ped, map, tped, rv.tdt.dir, window.size=0, snp.pos.df, param, window.type  = "M"){

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
                curr.window.result<-.runRV_TDTOnWindow(input.filepaths,rv.tdt.dir,param)
                pos.info<-.getWindowPos(start.index,mid.index,end.index,snp.pos.df)
                results.df[i,1:7]<-curr.window.result
                results.df[i,8:10]<-pos.info

        }

        return(results.df)

}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow<-function(ped, map, tped, window.type,start.index,end.index,i){

        #TODO: fix this to get the directory functions.R is in
        data.dir<-"./"

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

.runRV_TDTOnWindow<-function(input.filepaths,rv.tdt.dir,param){

        .calculateRV_TDTOnWindow(input.filepaths,rv.tdt.dir, param)
        results<-.extract.results()
        .clean.up.rv_tdt(input.filepaths)
        return(results)

}

############################

.calculateRV_TDTOnWindow<-function(input.filepaths,rv.tdt.dir,param){

        #TODO: Fix u to be 0.01
        #TODO: Add in all the parameters that RV-TDT includes
        #TODO: Modify command accordingly

        adapt<-param[1]

        alpha<-param[2]

        permut<-param[3]

        lower_cutoff<-param[4]

        upper_cutoff<-param[5]

        minVariants<-param[6]

        maxMissRatio <-param[7]

        filepath.tped<-paste0("'",input.filepaths[1],"'")
        filepath.phen<-paste0("'",input.filepaths[2],"'")
        filepath.map<-paste0("'",input.filepaths[3],"'")
		rv.tdt.dir<-paste0("'",rv.tdt.dir,"'")

        #filepath.tped<-paste0(input.filepaths[1])
        #filepath.phen<-paste0(input.filepaths[2])
        #filepath.map<-paste0(input.filepaths[3])

        gene.name<-"NA"

        #TODO: Add directories for each window?
        rv.tdt.results.dir<-paste0("./",gene.name)
        rv.tdt.results.dir

        command<-paste0(rv.tdt.dir, " ", rv.tdt.results.dir,
                        " -G ", filepath.tped,
                        " -P ", filepath.phen,
                        " -M ", filepath.map,
                        " --adapt ", adapt,
                        " --alpha ", alpha,
                        " --permut ", permut,
                        " --lower_cutoff ", lower_cutoff,
                                " --upper_cutoff ", upper_cutoff,
                                    " --minVariants ", minVariants,
                                    " --maxMissRatio ", maxMissRatio
        )

        #print(command)
        system(command)

}

############################

.getWindowPos<-function(start.index,mid.index,end.index,snp.pos.df){

	# Get start, end, middle position of each window
	start.pos<-as.character(snp.pos.df$pos[start.index])
	mid.pos<-as.character(snp.pos.df$pos[start.index])
	end.pos<-as.character(snp.pos.df$pos[end.index])
	pos.info<-c(start.pos,mid.pos,end.pos)

	# Add to results df for easier graphing
	return(pos.info)

}

########################################################

# Helper Functions - .calculateRV_TDTOnWindow

########################################################

.extract.results<-function(){
        gene.name<-"NA"
        filepath.results<-paste0("./",gene.name,"_pval/",gene.name,".pval")
        filepath.results
        pval.df<-read.table(filepath.results,comment.char = "",header=TRUE)
        pval.df
        return(pval.df)

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

        #print(delete.input.files)
        system(delete.input.files)

        #delete all output files
        currwd<-getwd()
        filepath.pval<-paste0("'",currwd,"/NA_pval/NA.pval","'")
        filepath.results<-paste0("'",currwd,"/NA_rvTDT/NA.rvTDT","'")
        delete.output.files<-paste("rm", filepath.pval, filepath.results)
        #print(delete.output.files)
        system(delete.output.files)

        #delete all output directories
        dir.pval<-paste0("'",currwd,"/NA_pval","'")
        dir.results<-paste0("'",currwd,"/NA_rvTDT","'")
        delete.output.dir<-paste("rmdir", dir.pval, dir.results)
        #print(delete.output.dir)
        system(delete.output.dir)

}
