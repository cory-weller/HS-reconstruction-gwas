#!/usr/bin/env Rscript
#doSCRM.Rscript
library(scrm)
library(data.table)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)
referenceGroup <- args[1]
N0_ancestral <- as.numeric(args[2])
mu <- as.numeric(args[3])
chromosomeSize <- as.numeric(args[4])
nHaplotypes <- as.numeric(args[5])
chromosome_i <- as.numeric(args[6])


rr = 1.5 ### units = cM/Mb 
r <- rr * (1/100) * (1/chromosomeSize)

rho <- 4*N0_ancestral * r
theta <- 4*N0_ancestral*mu*chromosomeSize

options(scipen = 999)

# format command for SCRM
scrm.command.internal <- paste(nHaplotypes, "1", "-t", theta, "-r", rho, chromosomeSize, "-p 7 -n 1 1", "--transpose-segsites", sep=" ")

# run SCRM
scrm.list <- scrm(scrm.command.internal)$seg_sites

# format SCRM output as data.table
scrm.dt <- data.table(matrix(unlist(scrm.list), ncol=nHaplotypes, byrow=TRUE))

# format SCRM data.table as VCF

distances <- as.numeric(colnames(scrm.list[[1]]))
setnames(scrm.dt, as.character(1:nHaplotypes))
scrm.dt[,"chr" := paste("chr", chromosome_i, sep="")]
scrm.dt[,POS := round(chromosomeSize*distances)]
scrm.dt[!duplicated(POS)]
scrm.dt <- scrm.dt[!duplicated(POS)]
ref_alt <- sample(c("A_C", "A_G", "A_T", "C_A", "C_G", "C_T", "G_A", "G_C", "G_T", "T_A", "T_C", "T_G"), size=dim(scrm.dt)[1], replace=TRUE)
scrm.dt[, c("REF", "ALT") := tstrsplit(ref_alt, "_")]
scrm.dt[,"ID":= paste(chr, POS, "SNP", sep="_")]
scrm.dt[,"QUAL":= 999]
scrm.dt[,"FILTER":= "PASS"]
scrm.dt[,"INFO":= "REFCOUNT=100;ALTCOUNT=100"]
scrm.dt[,"FORMAT":= "GT"]
setnames(scrm.dt, c("chr"), c("#CHROM"))
setcolorder(scrm.dt, c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",1:nHaplotypes))
setnames(scrm.dt, c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", paste("line_", 1:nHaplotypes, sep="")))

# write output to file

write.table(scrm.dt, file="haplotypes.vcf", sep="\t", quote=FALSE, row.names=FALSE, col.names=ifelse(chromosome_i==1, TRUE, FALSE))
