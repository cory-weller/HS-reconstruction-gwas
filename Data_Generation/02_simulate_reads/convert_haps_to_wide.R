#!/usr/bin/env Rscript

library(data.table)
library(foreach)

args <- commandArgs(trailingOnly=TRUE)

haplotypes_file <- args[1]
population <- args[2]
first_ind <- as.numeric(args[3])
last_ind <- as.numeric(args[4])
chr <- args[5]


fillVectorNAs <- function(x) {
  na.omit(x)[cumsum(!is.na(x))]
}

fillDataTableNAs <- function(DT, cols) {
  DT[, (cols) := lapply(.SD, function(x) fillVectorNAs(x)), .SDcols=cols]
}

convert_to_wide <- function(DT, chr) {
  # Collapse redundant breakpoints
  DT[, lineID_rleid := rleid(chromosome, haplotype, ind, lineID)]
  haps.collapsed <- DT[, list(chromosome,haplotype,start=min(start), stop=max(stop), ind, sex, lineID), by=lineID_rleid]
  haps.collapsed[,lineID_rleid := NULL]
  DT <- haps.collapsed[!duplicated(haps.collapsed)]
  

    unique_stops <- unique((sort(DT[,stop])))
    if(length(unique_stops) == 1) {
      unique_starts <- 1
    } else {
      unique_starts <- c(0,unique_stops[1:(length(unique_stops)-1)]) +1
    }
    
    DT.out <- data.table(start=unique_starts, stop=unique_stops)

    setkey(DT.out, start)
    setkey(DT, start)

    DT.out <- dcast(DT.out[DT], start+stop~haplotype, value.var="lineID")
    DT.out[, "chromosome" := chr]

    fillDataTableNAs(DT.out, c("1","2"))

    setcolorder(DT.out, c("chromosome","start","stop","1","2"))
    setnames(DT.out, c("1","2"), c("par1","par2"))
    return(DT.out[])
}


all_haplotypes <- fread(haplotypes_file)
all_haplotypes <- all_haplotypes[chromosome == chr]

for(ind_i in first_ind : last_ind) {
  haplotypes_file_out <- paste("./", population, "/", ind_i, ".haps.wide", sep="")
  ind_wide.dt <- convert_to_wide(all_haplotypes[ind==ind_i], chr=chr)
  write.table(ind_wide.dt, file=haplotypes_file_out, quote=F, row.names=F, col.names=T, sep="\t")
}
