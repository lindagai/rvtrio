#' Calculating RV-TDT statistic for user-defined windows
#'
#' `RV_TDT` returns a data frame containing the RV_TDT statistic for a VCF file. Note that `RV_TDT` only works for Linux and Mac OS X.
#'
#' @param vcf vcf file in 'CollapsedVCF' format

#' @param ped data frame containing pedigree information for the VCF

#' @param filepath.RV_TDT filepath to RV-TDT program. Can be downloaded here: \href{https://github.com/statgenetics/rv-tdt}{https://github.com/statgenetics/rv-tdt}.

#' @param window.type type of window, either number of markers ("M") or width of kilobase interval ("K") (doesn't work yet)

#' @param window.size size of window in number of markers

#' @param adapt To reduce computational time, adaptive permutation is used in /rvTDT/. Every /$adapt/ permutations (default: 500 permutations), the program will check if we should keep doing permutation (which means this gene looks promising to reach the desired alpha level), or we should give up on this gene (which means this gene will not reach the desired alpha level based on the permutations we have done so far, or we have done enough permutations)

#' @param alpha alpha level in adaptive permutation

#' @param permut The maximum number of permutations

#' @param lower_cutoff The cutoffs to determine which variants we should include in the analysis. In this example, the third column of map file is the number of minor allele counts, and here we only include the variants who have minor allele counts less than 100

#' @param upper_cutoff See lower_cutoff

#' @param minVariants The minimum number of variant sites for a gene. Genes with variant site number less than minVariants will be excluded from analysis (after check missing);

#' @param maxMissRatio The max missing ratio allowed for a variant. The variants with missing ratio greater than maxMissRatio will be excluded from analysis. In this example, we generated the genetic data file without any missing genotypes, so /--maxMissRatio 1/ is used here.

#' @return results data frame containing results from RV-TDT

#' @import dplyr VariantAnnotation
#' @importFrom splitstackshape cSplit
#' @importFrom methods is
#' @importFrom utils head read.table write.table file_test
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom BiocGenerics start
#'
#' @export
#'

#' @examples
#' \donttest{fp.ped <- system.file("inst", "extdata", "hg38.ped.txt", package = "rvtrio")
#' ped <- read.table(fp.ped,header=TRUE)
#' head(ped)
#' 
#' fp.vcf <- system.file("inst", "extdata", "hg38.vcf", package = "rvtrio")
#' hg.assembly <- "hg38"
#' vcf <- VariantAnnotation::readVcf(fp.vcf, hg.assembly)
#' 
#' RV_TDT.results <- rvtrio::RV_TDT(vcf, ped, filepath.RV_TDT = filepath.to.RV_TDT)
#' RV_TDT.results <- rvtrio::RV_TDT(vcf, ped, filepath.RV_TDT, window.size=0, window.type = "M")
#' }

########################################################

RV_TDT <- function(vcf, ped, filepath.RV_TDT, window.size=0, window.type = "M", adapt = 500, alpha = 0.00001, permut = 2000, lower_cutoff = 0, upper_cutoff = 100, minVariants = 3, maxMissRatio = 1){

        curr.wd <- getwd()

        .intializeEnv(vcf, ped, filepath.RV_TDT)

        parameters <- c(adapt, alpha, permut,lower_cutoff,
                        upper_cutoff, minVariants, maxMissRatio)

        snp.pos.df <- .getSnpPosDF(vcf)

		vcf <- geno(vcf)$GT

        tped <- .getTPED(vcf)
        ped <- .getPED(ped, vcf)
        map <- .getMAP(tped)

        results <- .runRV_TDT(ped, map, tped, filepath.RV_TDT,
                              window.size, snp.pos.df, param = parameters)

       .restoreEnv(curr.wd)

        return(results)
}

########################################################

# Helper Functions - RV_TDT

########################################################


.intializeEnv <- function(vcf, ped, filepath.RV_TDT) {

	    .checkInputs(vcf, ped, filepath.RV_TDT)
        rvtrio.dir <- file.path(.libPaths(),"rvtrio")
        setwd(rvtrio.dir)
        data.dir <-"./data"
        dir.create(file.path(data.dir), showWarnings = FALSE)

}

########################################################

# Helper Helper Functions - .intializeEnv

########################################################

.checkInputs <- function(vcf, ped, filepath.RV_TDT){

        .checkVCF(vcf)
        .checkPED(ped)
        .checkVCFandPED(vcf, ped)
        .checkFilepathToRV_TDT(filepath.RV_TDT)

}

##################################

.checkVCF <- function(vcf){
        if(!is(vcf, "CollapsedVCF")) {
                stop("VCF must be an object of class collapsedVCF")
        }
        else if(is(vcf, "CollapsedVCF")){
                 vcf <- geno(vcf)$GT
                 if(is.null(vcf))
                         stop("VCF does not seem to contain the genotype data.")
        }

        allowed.geno.values <- c("0|0", "0|1", "1|0", "1|1")
        genotype.entries <- unique(as.vector(vcf))

        if (sum(!(genotype.entries %in% allowed.geno.values)) > 0){
                stop("Genotype entries must be of the format `0|0`, `0|1`, `1|0`, or `1|1`. Did you remember to phase the VCF?")
        }
}

##################################

.checkPED <- function(ped){
	
        expected.cols.ped <- c("pid", "famid", "fatid", "motid", "sex", "affected")

        if(!setequal(expected.cols.ped, colnames(ped))){
                stop("PED must contain 6 columns named `pid`, `famid`, `fatid`, `motid`, `sex``, and `affected`. \n" ,
                     "These correspond to the personal ID, family ID, father ID, mother ID, \n", "sex (1 for male, 0 for female), and case(1)/control(0).")
        }

        ids.kid1 <- ped$fatid != 0
        ids.kid2 <- ped$motid != 0

        if(any(ids.kid1 != ids.kid2)){
                stop("fatid and motid must both be either zero or non-zero.")
        }

}

##################################

.checkVCFandPED <- function(vcf, ped){
        n.samples.vcf <- ncol(geno(vcf)$GT)
        n.samples.ped <- nrow(ped)
        vcf.ids <- colnames(geno(vcf)$GT)

        if(n.samples.vcf != n.samples.ped) {
                stop("The VCF and the PED have a different number of samples.")
        } else {
                if(!any(vcf.ids %in% ped$pid)) {
                        stop("The VCF has samples that are not present in the PED.")
                }
                if(!any(ped$pid %in% vcf.ids)) {
                        stop("The PED has samples that are not present in the VCF.")
                }
        }

}

##################################

.checkFilepathToRV_TDT <- function(filepath.RV_TDT){
        if(!file_test("-f", filepath.RV_TDT)) {
                stop("Filepath given to RV-TDT does not exist.  \n Did you give the correct path to RV-TDT?")
        }

}

########################################################

.getSnpPosDF<- function(vcf) {

        snp.pos.df <- data.frame(cbind(names(vcf),
                                          BiocGenerics::start(
                                                  SummarizedExperiment::rowRanges(vcf)
                                                  )
                                          )
                                    )
        colnames(snp.pos.df) <- c("snp.name","pos")
        return(snp.pos.df)

}

############################

.getTPED <- function(vcf.geno){

        tped <- as.matrix(splitstackshape::cSplit(vcf.geno,
                                                  colnames(vcf.geno), c("|")))
        rownames(tped) <-rownames(vcf.geno)
        return(tped)

}

############################

.getPED <- function(ped, vcf){

	    ped <- ped %>%
	    		select("pid", "famid", "fatid", "motid", "sex", "affected")

        pids <- data.frame(colnames(vcf), stringsAsFactors = FALSE)
        colnames(pids) <- "pids"
        ped <- dplyr::left_join(pids, ped, by=c("pids" ="pid"))
        return(ped)

}

############################

.getMAP <- function(tped){

        gene.id <- "NA"
        n.snps <- nrow(tped)
        gene.id.vec <- rep(gene.id, n.snps)
        snps <- rownames(tped)
        mafs <- as.vector(rowMeans(tped))
        map <- data.frame(cbind(gene.id.vec, snps, mafs),
                             stringsAsFactors=FALSE)
        return(map)

}

############################

.restoreEnv <- function(curr.wd){

       .deletePED()
       .deleteInputDir()
       .deleteResultsDir()
        setwd(curr.wd)


}

########################################################

# Helper Helper Function - .restoreEnv

########################################################

.deletePED <- function(){

        data.dir <- file.path(.libPaths(),"rvtrio", "data","input_files")
        filepath.ped <- file.path(data.dir, "pedfile.ped")
        command <- paste0("rm ", filepath.ped)
        system(command)

}

############################

.deleteResultsDir <- function(){

        results.dir <- file.path(.libPaths(),"rvtrio","data","results")
        delete.results.dir <- paste("rmdir", results.dir)
        system(delete.results.dir)

}

############################

.deleteInputDir <- function(){

        input.dir <- file.path(.libPaths(),"rvtrio","data","input_files")
        delete.input.dir <- paste("rmdir", input.dir)
        system(delete.input.dir)

}

############################

.runRV_TDT <- function(ped, map, tped, filepath.RV_TDT, window.size=0, snp.pos.df, param, window.type  = "M"){

        n.snps<- nrow(map)

        #return rv_tdt results as df
        if (window.size == 0){
                n.windows <- 1
                window.size <- n.snps
        } else {
                n.windows <- max(n.snps - window.size + 1,1)
        }

        results.df <- data.frame(matrix(data=NA, nrow=n.windows, ncol=10))
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
                                                         filepath.RV_TDT, param)
                pos.info <- .getWindowPos(start.index, mid.index,
                                          end.index, snp.pos.df)
                results.df[i,1:7] <- curr.window.result
                results.df[i,8:10] <- pos.info

        }

        #Remove gene.id column, which is not applicable here
        results.df <- results.df[, -1]

        return(results.df)

}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow <- function(ped, map, tped, window.type, start.index, end.index, i){

        data.dir <-"./data/input_files/"

        dir.create(file.path(data.dir), showWarnings = FALSE)

        file.param <- paste0("window",i,".",
                             start.index, "-", end.index, window.type)

        filepath.map.new <- paste0(data.dir, file.param, ".map")
        filepath.tped.new <- paste0(data.dir, file.param, ".tped")
        filepath.ped <- paste0(data.dir, "pedfile.ped")

        sm.map <- map[start.index:end.index,]
        sm.tped <- tped[start.index:end.index,]

        write.table(sm.map,filepath.map.new, sep="\t",
                    col.names=FALSE, row.names = FALSE, quote = FALSE)
        write.table(sm.tped,filepath.tped.new, sep="\t",
                    col.names=FALSE, row.names = TRUE, quote = FALSE)

        if (i==1){
                head(ped)
                write.table(ped, filepath.ped, sep="\t",
                            col.names=FALSE, row.names = FALSE, quote = FALSE)
        }

        filepaths <- c(filepath.tped.new, filepath.ped, filepath.map.new)

}

############################

.runRV_TDTOnWindow <- function(input.filepaths, filepath.RV_TDT,param){

        .calculateRV_TDTOnWindow(input.filepaths, filepath.RV_TDT, param)
        results <- .extractResults()
        .cleanUpRV_TDT(input.filepaths)
        return(results)

}

############################

.calculateRV_TDTOnWindow <- function(input.filepaths, filepath.RV_TDT, param){

        adapt <- param[1]

        alpha <- param[2]

        permut <- param[3]

        lower_cutoff <- param[4]

        upper_cutoff <- param[5]

        minVariants <- param[6]

        maxMissRatio <- param[7]

        #NOTE: Do *NOT* attempt to replacte paste0 with file.path
		#file.path will break if there are spaces in the filepaths
        filepath.tped <- paste0("'", input.filepaths[1], "'")
        filepath.phen <- paste0("'", input.filepaths[2], "'")
        filepath.map <- paste0("'", input.filepaths[3], "'")
        filepath.RV_TDT <- paste0("'", filepath.RV_TDT, "'")

        results.dir <- "./data/results/"
        dir.create(file.path(results.dir), showWarnings = FALSE)
        gene.name <- "NA"
        rv.tdt.results.dir <- paste0(results.dir, gene.name)

        command<-paste0(filepath.RV_TDT, " ", rv.tdt.results.dir,
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

        system(command)

}

############################

.getWindowPos <- function(start.index,mid.index,end.index,snp.pos.df){

        start.pos <- as.numeric(as.character(snp.pos.df$pos[start.index]))
        mid.pos <- as.numeric(as.character(snp.pos.df$pos[start.index]))
        end.pos <- as.numeric(as.character(snp.pos.df$pos[end.index]))
        pos.info <- c(start.pos, mid.pos, end.pos)
        return(pos.info)

}

########################################################

# Helper Functions - .calculateRV_TDTOnWindow

########################################################

.extractResults <- function(){

        gene.name<-"NA"
        results.dir <- paste0("./data/results/")
        filepath.results <- paste0(results.dir, gene.name, "_pval/", gene.name, ".pval")
        filepath.results
        pval.df <- read.table(filepath.results, comment.char = "", header=TRUE)
        pval.df
        return(pval.df)

}

############################

.cleanUpRV_TDT <- function(input.filepaths){

        ### DELETE INPUT FILES
        	#NOTE: Do *NOT* attempt to replacte paste0 with file.path
		#file.path will break if there are spaces in the filepaths
        filepath.tped <- paste0("'", input.filepaths[1], "'")
        filepath.map <- paste0("'", input.filepaths[3], "'")

        delete.input.files <- paste("rm",
                                    filepath.tped,
                                    filepath.map
        )

        system(delete.input.files)

        ### DELETE OUTPUT FILES
        #results.dir <- paste0("./data/results")
        results.dir <- file.path(.libPaths(),"rvtrio","data","results")

        filepath.pval <- file.path(results.dir,
                                   "NA_pval", "NA.pval")
        filepath.results <- file.path(results.dir,
                                      "NA_rvTDT", "NA.rvTDT")
        delete.output.files <- paste("rm", filepath.pval, filepath.results)
        system(delete.output.files)

        ### DELETE OUTPUT DIRECTORIES
        dir.pval <- file.path(results.dir, "NA_pval")
        dir.results <- file.path(results.dir, "NA_rvTDT")
        delete.output.dir <- paste("rmdir", dir.pval, dir.results)
        system(delete.output.dir)

}
