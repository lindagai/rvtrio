#' Calculating RV-TDT statistic for user-defined windows

#' `RV_TDT()` returns a data frame containing the RV_TDT statistic for a VCF or PLINK ped file. Note that RV_TDT only works for Linux and Mac OS X.

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
#' @import splitstackshape

#' @export
#'
#' @examples RV_TDT(filepath.vcf, filepath.ped)
#' @examples RV_TDT(filepath.vcf, filepath.ped, window.size=10, window.type="M")
#'
#' @export
#'

########################################################

RV_TDT<-function(vcf, vcf.ped, rv.tdt.dir, window.size=0, window.type = "M", adapt = 500, alpha = 0.00001, permut = 2000, lower_cutoff = 0, upper_cutoff = 100, minVariants = 3, maxMissRatio = 1){
        
        #.checkInputs(vcf, vcf.ped, rv.tdt.dir)
        
        parameters <- c(adapt, alpha, permut,lower_cutoff,
                        upper_cutoff, minVariants,maxMissRatio)

        snp.pos.df <- .get.snp.pos.df(vcf) 
        
		vcf <- geno(vcf)$GT
		
        tped <- .getTPED(vcf)
        ped <- .getPED(ped, tped)
        map <- .getMAP(tped)
        
        results <- .runRV_TDT(ped, map, tped, rv.tdt.dir,
                              window.size, snp.pos.df, param = parameters)
                              
      	#Remove references to the gene.region.id, which is not necessary
        results.df <- results.df[,-1]
        
        #Delete the new pedfile -- put this in a function if it works
        data.dir <-paste0(.libPaths(),"/rvtrio/data/input_files/")
        filepath.ped <- paste0(data.dir, "pedfile.ped")
        filepath.ped
        command <- paste0("rm ", filepath.ped)
        system(command)
        
        return(results)
}

########################################################

# Helper Functions - RV_TDT

########################################################

.checkInputs <- function(vcf, vcf.ped, rv.tdt.dir){

	#Check vcf
	if(!is(vcf, "CollapsedVCF")) {
			stop("VCF must be an object of class collapsedVCF")
	} else if(is(vcf, "CollapsedVCF")){
		vcf <- VariantAnnotation::geno(vcf)$GT
		if(is.null(vcf))
			stop("VCF does not seem to contain the genotype data.")
	}
		
	allowed.geno.values <- c("0|0", "0|1", "1|0", "1|1")
	
    if (sum(!(as.vector(geno(vcf)$GT) %in% allowed.geno.values)) > 0){
                stop("Genotype entries must be of the format `0|0`, `0|1`, `1|0`, or `1|1`. Did you remember to phase the VCF?")
        }
	
	#Check vcf's ped
	if(length(colnames(vcf.ped)) !=6 ) {
		stop("The PED file must have 6 columns named .")
		}
	#Check rv.tdt.dir
	if(!file.exists(rv.tdt.dir)) {
		stop("Did you give the correct path to RV-TDT?")
		}
}

############################

.get.snp.pos.df<- function(vcf) {
        
        snp.pos.df <- as.data.frame(cbind(names(vcf), 
                                          start(rowRanges(vcf))))
        colnames(snp.pos.df) <- c("snp.name","pos")
        return(snp.pos.df)
        
}

############################

.getPED <- function(ped, tped){

        pids <- as.data.frame(colnames(tped))
        colnames(pids) <- "pids"
        ped <- dplyr::left_join(pids, ped, by=c("pids" ="pid")) 
        return(ped)
}

############################

.getTPED <- function(vcf.geno){
        
        tped <- as.matrix(splitstackshape::cSplit(vcf.geno, 
                                                  colnames(vcf.geno), c("|")))
        rownames(tped) <-rownames(vcf.geno)
        return(tped)
}

############################

.getMAP <- function(tped){
        
        gene.id <- "NA"
        n.snps <- nrow(tped)
        gene.id.vec <- rep(gene.id, n.snps)
        snps <- rownames(tped)
        mafs <- as.vector(rowMeans(tped))
        map <- as.data.frame(cbind(gene.id.vec, snps, mafs),
                             stringsAsFactors=FALSE)
        return(map)
        
}

############################

.runRV_TDT <- function(ped, map, tped, rv.tdt.dir, window.size=0, snp.pos.df, param, window.type  = "M"){
        
        n.snps<- nrow(map)
        
        #return rv_tdt results as df
        if (window.size == 0){
                n.windows <- 1
                window.size <- n.snps
        } else {
                n.windows <- max(n.snps - window.size + 1,1)
        }
        
        results.df <- as.data.frame(matrix(data=NA, nrow=n.windows, ncol=10))
        colnames(results.df) <- c(
                "gene.name",
                "CMC.Analytical","BRV.Haplo","CMC.Haplo",
                "VT.BRV.Haplo","VT.CMC.Haplo","WSS.Haplo",
                "start.pos","mid.window.pos","end.pos"
        )
        
        for (i in (1:n.windows)){
                start.index <- i
                end.index <- min(i + window.size-1, n.snps)
                mid.index <- floor(start.index + (end.index - start.index)/2)
                input.filepaths <- .getInputFilesForWindow(ped, map, tped,
                                                           window.type,
                                                           start.index, end.index,
                                                           i)
                curr.window.result <- .runRV_TDTOnWindow(input.filepaths,
                                                         rv.tdt.dir, param)
                pos.info <- .getWindowPos(start.index, mid.index,
                                          end.index, snp.pos.df)
                results.df[i,1:7] <- curr.window.result
                results.df[i,8:10] <- pos.info
                
        }
        
        return(results.df)
        
}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow <- function(ped, map, tped, window.type, start.index, end.index, i){
        
        #TODO: fix this to get the directory functions.R is in
        #data.dir <- "./"
        data.dir <-paste0(.libPaths(),"/rvtrio/data/input_files/")
        
        file.param <- paste0("window",i,".",
                             start.index, "-", end.index, window.type)
        file.param
        
        filepath.map.new <- paste0(data.dir, file.param, ".map")
        filepath.tped.new <- paste0(data.dir, file.param, ".tped")
        filepath.tped.new
        filepath.ped <- paste0(data.dir, "pedfile.ped")
        filepath.ped
        
        sm.map <- map[start.index:end.index,]
        sm.tped <- tped[start.index:end.index,]
        
        write.table(sm.map,filepath.map.new, sep="\t", 
                            col.names=FALSE, row.names = FALSE, quote = FALSE)
        write.table(sm.tped,filepath.tped.new, sep="\t", 
                            col.names=FALSE, row.names = FALSE, quote = FALSE)
        
        if (i==1){
                head(ped)
                write.table(ped, filepath.ped, sep="\t", 
                            col.names=FALSE, row.names = FALSE, quote = FALSE)
        }
        
        filepaths <- c(filepath.tped.new, filepath.ped, filepath.map.new)
        
}

############################

.runRV_TDTOnWindow <- function(input.filepaths, rv.tdt.dir,param){

        .calculateRV_TDTOnWindow(input.filepaths, rv.tdt.dir, param)
        results <- .extract.results()
        .clean.up.rv_tdt(input.filepaths)
        return(results)

}

############################

.calculateRV_TDTOnWindow <- function(input.filepaths, rv.tdt.dir, param){

        adapt <- param[1]

        alpha <- param[2]

        permut <- param[3]

        lower_cutoff <- param[4]

        upper_cutoff <- param[5]

        minVariants <- param[6]

        maxMissRatio <- param[7]

        filepath.tped <- paste0("'",input.filepaths[1],"'")
        filepath.phen <- paste0("'",input.filepaths[2],"'")
        filepath.map <- paste0("'",input.filepaths[3],"'")
		rv.tdt.dir <- paste0("'",rv.tdt.dir,"'")

        #filepath.tped<-paste0(input.filepaths[1])
        #filepath.phen<-paste0(input.filepaths[2])
        #filepath.map<-paste0(input.filepaths[3])

        #gene.name <- "NA"
        #rv.tdt.results.dir<-paste0("./",gene.name)
        data.dir <-paste0(.libPaths(),"/rvtrio/data/results/", gene.name)
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

        print(command)
        system(command)

}

############################

.getWindowPos<-function(start.index,mid.index,end.index,snp.pos.df){

	start.pos<-as.character(snp.pos.df$pos[start.index])
	mid.pos<-as.character(snp.pos.df$pos[start.index])
	end.pos<-as.character(snp.pos.df$pos[end.index])
	pos.info<-c(start.pos,mid.pos,end.pos)
	return(pos.info)

}

########################################################

# Helper Functions - .calculateRV_TDTOnWindow

########################################################

.extract.results <- function(){
	
        gene.name<-"NA"
        filepath.results <- paste0("./",gene.name,"_pval/",gene.name,".pval")
        filepath.results
        pval.df <- read.table(filepath.results,comment.char = "",header=TRUE)
        pval.df
        return(pval.df)

}

############################

.clean.up.rv_tdt <- function(input.filepaths){

        #delete all input files
        filepath.tped <- paste0("'", input.filepaths[1], "'")
        filepath.map <- paste0("'", input.filepaths[3], "'")

        delete.input.files <- paste("rm",
                                  filepath.tped,
                                  filepath.map
        )

        #print(delete.input.files)
        system(delete.input.files)

        #delete all output files
        currwd <- getwd()
        filepath.pval <- paste0("'", currwd, "/NA_pval/NA.pval","'")
        filepath.results <- paste0("'", currwd, "/NA_rvTDT/NA.rvTDT","'")
        delete.output.files <- paste("rm", filepath.pval, filepath.results)
        #print(delete.output.files)
        system(delete.output.files)

        #delete all output directories
        dir.pval <- paste0("'", currwd, "/NA_pval", "'")
        dir.results <- paste0("'",currwd, "/NA_rvTDT", "'")
        delete.output.dir <- paste("rmdir", dir.pval, dir.results)
        #print(delete.output.dir)
        system(delete.output.dir)

}
