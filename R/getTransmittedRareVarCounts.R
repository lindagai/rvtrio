################################################################################

# Get transmitted RV counts

################################################################################

getTransmittedRareVarCounts <- function(vcf, ped){
	
	vcf <- geno(vcf)$GT	
	mafs <- .getMAFs(vcf)
	rare.var.geno <- .getRareVarGeno(mafs)
	transmitted.rare.var.ct <- .getTransmittedRareVar(rare.var.geno, ped)
	return(transmitted.rare.var.ct)
	
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

.getTransmittedRareVar <- function(rare.var.geno, ped){
        
        #Convert rare.var.geno to allele counts - OK
        rare.var.geno <- .getGenotypeMatrix(rare.var.geno)
        
        #Get children with rare variants - OK
        child.with.RV.geno <- .getRareVarInChildren(rare.var.geno, ped)
        child.rv.ct <- rowSums(child.with.RV.geno)
        
        #Allocate table  - OK
        rv.in.children <- rownames(child.with.RV.geno)
        rare.snp.table <- as.data.frame(rep(NA,length(rv.in.children)
        )
        )
        
        rownames(rare.snp.table) <- rv.in.children
        colnames(rare.snp.table) <- "trans.ct"
        
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
                                filter(pid %in% children.with.this.rv) %>%
                                select(fatid, motid)
                        
                        fatid <- as.character(parents.of.children.with.this.rv$fatid)
                        motid <- as.character(parents.of.children.with.this.rv$motid)
                        
                        parent1.rv.ct <- rare.var.geno[curr.rv, fatid]
                        parent2.rv.ct <- rare.var.geno[curr.rv, motid]
                        
                        parents.rv.ct <- parent1.rv.ct + parent2.rv.ct
                        
                        n.trans.SNPs.trio <- min(child.rv.ct, parents.rv.ct)
                        n.trans.SNPs <- n.trans.SNPs + n.trans.SNPs.trio
                }
                
                n.trans.SNPs
                rare.snp.table[i,] <- n.trans.SNPs      
                
        }
        
        return(rare.snp.table)
}

       
################################################################################

.getRareVarInChildren <- function(rare.var.geno, ped){
        
        child.pids <- ped %>%
                filter(fatid!=0) %>%
                filter(motid!=0) %>%
                select(pid) %>%
                sapply(as.character) %>%
                as.vector
        
        study.participants <- colnames(rare.var.geno)
        child.geno <- rare.var.geno[, study.participants %in% child.pids]
        mafs.in.children <- .getMAFs(child.geno)
        children.with.RV.geno <- child.geno[mafs.in.children > 0, ]
        return(children.with.RV.geno)
        
}