#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
    population <- args[1]
    chromosome <- args[2]
    ind_n <- args[3]
    RABBIT_MAX_FOUNDERS <- as.numeric(args[4])
    RABBIT_MAX_SITES <- as.numeric(args[5])

library(data.table)
library(foreach)


# in R, combine information from three sources:
# 1. Readcounts, to find out which SNPs is the whole set we're considering
# 2. population VCF file, for all variants
# 3. alleleFreqs, for information content
# 4. recombination .bed file

# 3. Of remaining SNPs, merge in read depth and allele frequency
# 4. Subset by top N, according to criteria deemed important
# 5. Re-organize by chromosome + POS
# 6. Iterate over columns and print to file (effectively transposing without holding in memory)


# 1. read in population VCF, readcounts, and freqs for this individual


founders <- readLines(paste('../01_forward_simulator/', population, '.founders', sep=""))
vcf <- fread(paste(population, "/", chromosome, ".RABBIT.vcf", sep=""),
              skip="#CHROM",
              header=TRUE,
              showProgress=FALSE)

setnames(vcf, "#CHROM", "CHROM")
vcf[, "ID" := paste(CHROM, POS, "SNP", sep="_")]
setkey(vcf, "CHROM", "POS", "ID")

readcounts <- fread(paste(population, "/readcounts/", ind_n, ".", chromosome, ".readcounts", sep=""),
                    header=TRUE,
                    select=c("contig","position","variantID","refCount","altCount","otherBases","improperPairs"),
                    showProgress=FALSE)

setkey(readcounts, "contig", "position", "variantID")

mla <- fread(paste(population, "/MLA/", ind_n, ".", chromosome, ".mla", sep=""), header=TRUE, showProgress=FALSE)
mla[,rank := frank(-N, ties.method="random"), by=chromosome]
mla.chosen <- mla[rank <= RABBIT_MAX_FOUNDERS]


# Read in .bed file containing recombination rates
# Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
bed <- fread( file="../input_data/recombination_map.bed",
              header=FALSE,
              showProgress=FALSE,
              col.names=c("chr","start","stop","c")
            )

# subset chromosome
bed <- bed[chr == chromosome]

# Key .bed file by chromosome and window start
setkey(bed, chr, start)

# Generate column for cM (centiMorgans)
bed[, cM := c*((stop-start)/1e6)]

# Calculate cumulative cM at the END of each window
bed[, cumulative_cM := cumsum(cM), by=chr]

# Generate functions (cM ~ BP) to translate base pair position (in VCF) to units of cM
recombination_function <- new.env()

# Create recombination function based on .bed file
for(chromosome in unique(bed[,chr])) {
    recombination_function[[as.character(chromosome)]] <- approxfun(c(0, bed[chr==chromosome][,stop]), c(0,bed[chr==chromosome][,cumulative_cM]))
}

options(scipen=999)

for(chr.i in unique(bed[,chr])) {
    cat(chr.i)
    cat("\n")
    chosen.founders <- sort(mla.chosen[chromosome==chr.i][,lineID])

    # If no most likely founders for this chromosome, then skip. This should not happen.
    if(length(chosen.founders)==0) {
        cat("No most likely founders for chromosome ")
        cat(chr.i)
        cat("\n")
        next
    }

    # if only one founder was chosen for the entire chromosome (e.g., homozygous for entire chromosome for a single founder)
    if(length(chosen.founders)==1) {
        writeLines(paste("Nonrecombinant Homozygous ", chosen.founders, sep=""), con=paste(population, "/RABBIT_output/", ind_n, ".", chr.i, ".RABBIT.out.csv", sep=""))
        next
    }

    vcf.sub <- copy(vcf[CHROM==chr.i][,c("CHROM","POS","ID", chosen.founders), with=FALSE])

    vcf.sub[, nRef := apply(.SD, 1, function(x) sum(x == 0, na.rm=TRUE)), .SDcols=chosen.founders]
    vcf.sub[, nAlt := apply(.SD, 1, function(x) sum(x == 1, na.rm=TRUE)), .SDcols=chosen.founders]
    vcf.sub[, refFreq := nRef/(nRef+nAlt)]
    vcf.sub[, altFreq := nAlt/(nRef+nAlt)]
    setkey(vcf.sub, ID)

    # Merge in sequenced allele freq
    vcf.sub.merge <- merge(vcf.sub, readcounts, by.x="ID", by.y="variantID")[refFreq != 1 & altFreq != 1 & otherBases == 0 & improperPairs == 0]
    vcf.sub.merge[, c("otherBases","improperPairs") := NULL]
    vcf.sub.merge[, "imputed_ind" := ifelse(refCount > 0 & altCount > 0, "12",
                               ifelse(altCount==0, "1N",
                               ifelse(refCount==0, "2N", "NN")))]

    # Select up to maximum number of SNPs
    vcf.sub.merge[,indx := 1:.N]

    # Retain ALL heterozygous sites
    vcf.sub.merge[refCount>0 & altCount>0, marker := TRUE]
    vcf.sub.merge[, freq := min(refFreq, altFreq), by=indx]
    nSites <- dim(vcf.sub.merge[is.na(marker)])[1]
    samplePercent <- 0.1
    sampleNSites <- trunc(nSites * samplePercent)
    retainmentPercent <- 0.002
    retainNSites <- trunc(nSites * retainmentPercent)



    if(dim(vcf.sub.merge)[1] <= RABBIT_MAX_SITES) {
        vcf.sub.merge[, marker := TRUE]
    } else {
        while(dim(vcf.sub.merge[marker == TRUE])[1] < RABBIT_MAX_SITES ) {
            indicesToMark <- vcf.sub.merge[is.na(marker)][sample(.N, size=sampleNSites, replace=TRUE)][order(freq)][1:retainNSites][,indx]
            vcf.sub.merge[indicesToMark, marker := TRUE]
        }
    }

    markers <- vcf.sub.merge[marker==TRUE][sample(.N, size=min(RABBIT_MAX_SITES, .N), replace=FALSE)]
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
                con=paste(population, "/RABBIT_input/", ind_n, ".", chr.i, ".RABBIT.in", sep="")
              )

    write.table(t(as.matrix(markers)),
                file=paste(population, "/RABBIT_input/", ind_n, ".", chr.i, ".RABBIT.in", sep=""),
                quote=FALSE,
                row.names=TRUE,
                col.names=FALSE,
                sep=",",
                na="N",
                append=TRUE)
}