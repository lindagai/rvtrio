% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RV_TDT.R
\name{RV_TDT}
\alias{RV_TDT}
\title{Calculating RV-TDT statistic for user-defined windows}
\usage{
RV_TDT(vcf, ped, filepath.RV_TDT, window.size = 0, window.type = "M",
  adapt = 500, alpha = 1e-05, permut = 2000, lower_cutoff = 0,
  upper_cutoff = 100, minVariants = 3, maxMissRatio = 1)
}
\arguments{
\item{vcf}{vcf file in 'CollapsedVCF' format}

\item{ped}{data frame containing pedigree information for the VCF}

\item{filepath.RV_TDT}{filepath to RV-TDT program. Can be downloaded here: \href{https://github.com/statgenetics/rv-tdt}{https://github.com/statgenetics/rv-tdt}.}

\item{window.size}{size of window in number of markers}

\item{window.type}{type of window, either number of markers ("M") or width of kilobase interval ("K") (doesn't work yet)}

\item{adapt}{To reduce computational time, adaptive permutation is used in /rvTDT/. Every /$adapt/ permutations (default: 500 permutations), the program will check if we should keep doing permutation (which means this gene looks promising to reach the desired alpha level), or we should give up on this gene (which means this gene will not reach the desired alpha level based on the permutations we have done so far, or we have done enough permutations)}

\item{alpha}{alpha level in adaptive permutation}

\item{permut}{The maximum number of permutations}

\item{lower_cutoff}{The cutoffs to determine which variants we should include in the analysis. In this example, the third column of map file is the number of minor allele counts, and here we only include the variants who have minor allele counts less than 100}

\item{upper_cutoff}{See lower_cutoff}

\item{minVariants}{The minimum number of variant sites for a gene. Genes with variant site number less than minVariants will be excluded from analysis (after check missing);}

\item{maxMissRatio}{The max missing ratio allowed for a variant. The variants with missing ratio greater than maxMissRatio will be excluded from analysis. In this example, we generated the genetic data file without any missing genotypes, so /--maxMissRatio 1/ is used here.}
}
\value{
results data frame containing results from RV-TDT
}
\description{
`RV_TDT` returns a data frame containing the RV_TDT statistic for a VCF file. Note that `RV_TDT` only works for Linux and Mac OS X.
}
\examples{
fp.ped <- system.file("extdata", "hg38.ped.txt", package = "rvtrio")
ped <- read.table(fp.ped,header=TRUE)
head(ped)

fp.vcf <- system.file("extdata", "hg38.vcf", package = "rvtrio")
hg.assembly <- "hg38"
vcf <- VariantAnnotation::readVcf(fp.vcf, hg.assembly)

RV_TDT.results <- rvtrio::RV_TDT(vcf, ped, filepath.RV_TDT = filepath.to.RV_TDT)
RV_TDT.results <- rvtrio::RV_TDT(vcf, ped, filepath.RV_TDT, window.size=0, window.type = "M")

}
