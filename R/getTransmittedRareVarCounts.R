#' Get counts of the number of times each rare variants is transmitted from parents to offspring
#'
#' `getTransmittedRareVarCounts()` returns a data frame containing the count of rare variants transmitted from the parents to the affected offspring.
#'
#' @param vcf vcf file
#'
#' @param ped data frame containing pedigree information for the VCF
#'
#' @param cutoff The MAF cutoff to determine which variants are "rare" (default is 0.01), for which transmission counts will calculated.  
#'
#' @return results data frame containing the names of the SNVs that are considered rare, the number of times the rare variant was transmitted from parents to the affected offspring in the dataset, and the genomic position of the rare SNPs.
#'
#' @import dplyr VariantAnnotation
#' @importFrom methods is
#' @importFrom utils head read.table write.table type.convert 
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom BiocGenerics start
#'
#'
#' @export
#'
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
#' #Get no. of rare variants that are transmitted from parents to offspring
#' transmitted.rare.var.ct <- rvtrio::getTransmittedRareVarCounts(vcf, ped)
#' }

################################################################################

getTransmittedRareVarCounts <- function(vcf, ped, cutoff = 0.01){
	
	.checkInputFiles(vcf, ped)
	snp.pos.df <- .getPositionsOfSnps(vcf)
	vcf <- geno(vcf)$GT	
	mafs <- .getMAFs(vcf)
	rare.var.geno <- .getRareVarGeno(vcf, mafs, cutoff)
	transmitted.rare.var.ct <- .getTransmittedRareVar(rare.var.geno, ped, snp.pos.df)
	return(transmitted.rare.var.ct)
	
}

################################################################################

.checkInputFiles <- function(vcf, ped){

        .checkVCF(vcf)
        .checkPED(ped)
        .checkVCFandPED(vcf, ped)

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

################################################################################

.getMAFs <- function(vcf){
	#Put this in its own file
        geno.matrix <- .getGenotypeMatrix(vcf)
        
        #Get MAF of each SNP
        #Divide by 2 because there are 2 alleles per individual
        mafs <- rowMeans(geno.matrix)/2
        
        return(mafs)
        
}

################################################################################

.getGenotypeMatrix <- function(vcf){
	
	    #Calculate MAFs
        geno.matrix <- vcf
        
        #Replace genotypes with minor allele counts
        geno.matrix[geno.matrix == "0|0"]<- 0
        geno.matrix[geno.matrix == "1|0"]<- 1
        geno.matrix[geno.matrix == "0|1"]<- 1
        geno.matrix[geno.matrix == "1|1"]<- 2
        
        #Convert from string to numeric vector
        geno.matrix <- type.convert(geno.matrix)
        
        return(geno.matrix)
	
}

################################################################################

.getRareVarGeno <- function(vcf, mafs, cutoff = 0.01){
        
        rare.var.snps <- names(mafs[which(mafs < cutoff & mafs > 0)])
                                    snps.in.geno <- rownames(vcf)
                                    rare.var.geno <- vcf[snps.in.geno %in% rare.var.snps,]
                                    return(rare.var.geno)
                                    
}

################################################################################

.getTransmittedRareVar <- function(rare.var.geno, ped, snp.pos.df){
        
        #Convert rare.var.geno to allele counts - OK
        rare.var.geno <- .getGenotypeMatrix(rare.var.geno)
        
        #Get children with rare variants - OK
        child.with.RV.geno <- .getRareVarInChildren(rare.var.geno, ped)
        child.rv.ct <- rowSums(child.with.RV.geno)
        
        #Allocate table  - OK
        rv.in.children <- rownames(child.with.RV.geno)
        rare.snp.table <- data.frame(matrix(ncol=2, nrow=length(rv.in.children)),
        stringsAsFactors=FALSE)
        
        colnames(rare.snp.table) <- c("rare.snp", "trans.ct")
        rare.snp.table$rare.snp <- rv.in.children
        
        for (i in 1:length(rv.in.children)){
                #Check if the as.character is necessary
                curr.rv <- as.character(rv.in.children[i])
                
                children.with.this.rv <- child.with.RV.geno[curr.rv,]
                children.with.this.rv <- names(children.with.this.rv[children.with.this.rv >=1]
                )
                n.trans.SNPs <- 0
                
                for (j in 1:length(children.with.this.rv)){
                        curr.child <- children.with.this.rv[j]
                        child.rv.ct <- child.with.RV.geno[curr.rv, curr.child]
                        child.rv.ct
                                             
                        parents.of.children.with.this.rv <- ped %>%
                                filter_(~pid %in% children.with.this.rv) %>%
                                select_(.dots = c('fatid', 'motid'))
                        
                        fatid <- as.character(parents.of.children.with.this.rv$fatid)
                        motid <- as.character(parents.of.children.with.this.rv$motid)
                        
                        parent1.rv.ct <- rare.var.geno[curr.rv, fatid]
                        parent2.rv.ct <- rare.var.geno[curr.rv, motid]
                        
                        parents.rv.ct <- parent1.rv.ct + parent2.rv.ct
                        
                        n.trans.SNPs.trio <- min(child.rv.ct, parents.rv.ct)
                        n.trans.SNPs <- n.trans.SNPs + n.trans.SNPs.trio
                }
                
                n.trans.SNPs
                rare.snp.table[i,2] <- n.trans.SNPs      
                
        }
        
        rare.snp.table <- left_join(rare.snp.table, snp.pos.df, by = c("rare.snp" = "snp"))
        return(rare.snp.table)
}

       
################################################################################

.getRareVarInChildren <- function(rare.var.geno, ped){
       
        child.pids <- ped %>%
                filter_(~fatid!=0) %>%
                filter_(~motid!=0) %>%
                select_(.dots = c("pid")) %>%
                sapply(as.character) %>%
                as.vector
        
        study.participants <- colnames(rare.var.geno)
        child.geno <- rare.var.geno[, study.participants %in% child.pids]
        mafs.in.children <- .getMAFs(child.geno)
        children.with.RV.geno <- child.geno[mafs.in.children > 0, ]
        return(children.with.RV.geno)
        
}

################################################################################

.getPositionsOfSnps <- function(vcf){
        
        #Add SNP name and position to DF of results
        snp <- names(vcf)
        pos <- BiocGenerics::start(rowRanges(vcf))
        snp.pos.df <- data.frame(snp, pos, stringsAsFactors = FALSE)
        return(snp.pos.df)
        
}