#!/usr/bin/env Rscript
# args <- c("4CB-1-73_S73", "2L", "16", "5000", "dgrp2.recombination_map.bed")
args <- commandArgs(trailingOnly=TRUE)
    BAM_FILESTEM <- args[1]
    CHROMOSOME <- args[2]
    RABBIT_MAX_FOUNDERS <- as.numeric(args[3])
    RABBIT_MAX_SITES <- as.numeric(args[4])

library(data.table)
library(foreach)

bed_files <- list.files(pattern="*.bed")
if(length(bed_files) > 1) {
    cat("More than one bed file in this directory.\n
    Should only contain a single bed file for recombination rates.\n")
    exit
} else {
    RECOMBINATION_BED_FILE <- bed_files[1]
}

# Import and format genotypes for founders for this chromosome
vcf <- fread(file = paste(CHROMOSOME, ".RABBIT.vcf", sep=""),
              skip = "#CHROM",
              header = TRUE,
              showProgress = FALSE)
setnames(vcf, "#CHROM", "CHROM")
vcf[, "ID" := paste(CHROM, POS, "SNP", sep="_")]
setkey(vcf, "CHROM", "POS", "ID")

# Import and format read counts for this individual's chromosome
readcounts <- fread(file = paste(BAM_FILESTEM, CHROMOSOME, "readcounts", sep="."),
                    header = TRUE,
                    select = c("contig","position","variantID","refCount","altCount","otherBases","improperPairs"),
                    showProgress = FALSE)
setkey(readcounts, "contig", "position", "variantID")

# Import most likely ancestors (MLA) file
mla <- fread(file = paste(BAM_FILESTEM, CHROMOSOME, "mla", sep="."),
            header = TRUE,
            showProgress = FALSE)
mla[,rank := frank(-N, ties.method="random"), by = chromosome]
mla.chosen <- mla[rank <= RABBIT_MAX_FOUNDERS]


# Read in .bed file containing recombination rates
# Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
bed <- fread(file = RECOMBINATION_BED_FILE,
              header = FALSE,
              showProgress = FALSE,
              col.names = c("chr","start","stop","c")
            )
# subset chromosome
bed <- bed[chr == CHROMOSOME]
setkey(bed, chr, start)

# Generate column for cM (centiMorgans)
bed[, cM := c*((stop-start)/1e6)]

# Calculate cumulative cM at the END of each window
bed[, cumulative_cM := cumsum(cM), by=chr]

# Generate functions (cM ~ BP) to translate base pair position (in VCF) to units of cM
recombination_function <- new.env()

# Create recombination function based on .bed file
for(chr.i in unique(bed[,chr])) {
    recombination_function[[as.character(chr.i)]] <- approxfun(c(0, bed[chr==chr.i][,stop]), c(0,bed[chr==chr.i][,cumulative_cM]))
}

options(scipen=999)


chosen.founders <- mla.chosen[,lineID]

# If no most likely founders for this chromosome, then skip. This should not happen.
if(length(chosen.founders)==0) {
    cat("No most likely founders for chromosome ")
    cat(CHROMOSOME)
    cat(". this should not happen.")
    cat("\n")
    exit
}

# if only one founder was chosen for the entire chromosome (e.g., homozygous for entire chromosome for a single founder)
if(length(chosen.founders)==1) {
    writeLines(paste("Nonrecombinant Homozygous ", chosen.founders, sep=""), con=paste(BAM_FILESTEM, CHROMOSOME, ".RABBIT.out.csv", sep=""))
    exit
}

vcf <- vcf[,c("CHROM","POS","ID", chosen.founders), with=FALSE]

vcf[, nRef := apply(.SD, 1, function(x) sum(x == 0, na.rm=TRUE)), .SDcols=chosen.founders]
vcf[, nAlt := apply(.SD, 1, function(x) sum(x == 1, na.rm=TRUE)), .SDcols=chosen.founders]
vcf[, refFreq := nRef/(nRef+nAlt)]
vcf[, altFreq := nAlt/(nRef+nAlt)]
setkey(vcf, ID)

# Merge in sequenced allele freq
vcf.merge <- merge(vcf, readcounts, by.x="ID", by.y="variantID")[refFreq != 1 & altFreq != 1 & otherBases == 0 & improperPairs == 0]
vcf.merge[, c("otherBases","improperPairs") := NULL]
vcf.merge[, "imputed_ind" := ifelse(refCount > 0 & altCount > 0, "12",
                            ifelse(altCount==0, "1N",
                            ifelse(refCount==0, "2N", "NN")))]

# Select up to maximum number of SNPs
vcf.merge[,indx := 1:.N]

# Retain ALL heterozygous sites
vcf.merge[refCount>0 & altCount>0, marker := TRUE]
vcf.merge[, freq := min(refFreq, altFreq), by=indx]
nSites <- dim(vcf.merge[is.na(marker)])[1]
samplePercent <- 0.1
sampleNSites <- trunc(nSites * samplePercent)
retainmentPercent <- 0.002
retainNSites <- trunc(nSites * retainmentPercent)



if(dim(vcf.merge)[1] <= RABBIT_MAX_SITES) {
    vcf.merge[, marker := TRUE]
} else {
    while(dim(vcf.merge[marker == TRUE])[1] < RABBIT_MAX_SITES ) {
        indicesToMark <- vcf.merge[is.na(marker)][sample(.N, size=sampleNSites, replace=TRUE)][order(freq)][1:retainNSites][,indx]
        vcf.merge[indicesToMark, marker := TRUE]
    }
}

markers <- vcf.merge[marker==TRUE][sample(.N, size=min(RABBIT_MAX_SITES, .N), replace=FALSE)]
setkey(markers, "CHROM", "POS")
markers[, c("contig", "position", "refCount", "altCount", "freq", "indx", "marker", "refFreq", "altFreq", "nRef", "nAlt") := NULL]
setnames(markers, "CHROM", "Chromosome")
setnames(markers, "ID", "SNP")
markers[, cM := recombination_function[[chr.i]](POS)][]
markers[, POS := NULL]

markers[, (chosen.founders) := lapply(.SD, "+", 1), .SDcols=chosen.founders]
line.names=colnames(markers)[3:(length(colnames(markers))-1)]
nFoundersUsed <- length(colnames(markers))-4

setcolorder(markers, c("SNP","Chromosome","cM", line.names))
writeLines( paste("#founders,",nFoundersUsed, sep=""),
            con=paste(BAM_FILESTEM, ".", CHROMOSOME, ".RABBIT.in", sep="")
            )

write.table(t(as.matrix(markers)),
            file=paste(BAM_FILESTEM, ".", CHROMOSOME, ".RABBIT.in", sep=""),
            quote=FALSE,
            row.names=TRUE,
            col.names=FALSE,
            sep=",",
            na="N",
            append=TRUE)