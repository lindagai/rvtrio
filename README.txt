######################################################

README.txt

######################################################

Description of folders and files in rvtrio folder.
Author: Linda Gai
Last update: 10/15/2019

######################################################

rvtrio

######################################################

rvtrio R package.
Contains wrapper functions for RV-TDT, calculating the number of transmitted rare variants from parents to child, and some unit tests.

=====

Contents:

=====

/data

(1) empty folder that RV_TDT.R uses to write the input and output files for RV-TDT, before erasing them

=====

/inst/extdata

Installed external datasets that rvtrio uses in unit tests and in roxygen2 documentation.

(1) hg38.ped.txt
Pedigree file from the cleaned Latino hg38 datasets.

(2) hg38.vcf
hg38 VCF from the cleaned Latino hg38 datasets, contains a small section from 8q24.

=====

/R

(1) getTransmittedRareVarCounts.R - calculates number of variants transmitted from parents to offspring 

(2) RV_TDT.R - wrapper function for RV-TDT

Each function is documented using roxygen2.

=====

/tests

(1) unit tests for RV_TDT using testthat

NOTE: getTransmittedRareVarCounts.R does not have unit tests.

=====

rvtrio.Rproj

Project file for package.

######################################################

Concerns and notes:

######################################################

(1) getTransmittedRareVarCounts.R does not have unit tests.

(2) In the roxygen2 documentation of RV_TDT and getTransmittedRareVarCounts, the examples are not checked (they are surrounded by /donttest{}).

For RV_TDT, this is because a file path to RV_TDT is required for the example to run. 

For getTransmittedRareVarCounts, this is because a file path to RV_TDT is required for the example to run. 