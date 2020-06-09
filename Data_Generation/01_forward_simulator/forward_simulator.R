#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(digest)

args <- commandArgs(trailingOnly=TRUE)

# Debugging args:
#   args <- c("-bed","recombination.bed","-prefix","dgrp","-n0","500","-rate","1.0","-sex","dioecious","-nfounders","32","-ngenerations","2","-lineIDs","lines.txt","-chrx","X","-iter","1","-recombination","femaleOnly","-dmel","TRUE","-nthreads","2","-replacement","TRUE")

args.env <- new.env()
for(i in seq(1, length(args), 2)) {
    args.env[[args[i]]] <- args[i+1]
}

# Required arguments:
bed_file <- args.env[["-bed"]]
prefix <- args.env[["-prefix"]]
popSize <- as.numeric(args.env[["-n0"]])
growthRate <- as.numeric(args.env[["-rate"]])
sexModel <- args.env[["-sex"]]
n_founders <- as.numeric(args.env[["-nfounders"]])
n_generations <- as.numeric(args.env[["-ngenerations"]])
lineIDs_filename <- args.env[["-lineIDs"]]
x.chromosome <- args.env[["-chrx"]]
recombinationModel <- args.env[["-recombination"]]
n_threads <- as.numeric(args.env[["-nthreads"]])
iteration <- args.env[["-iter"]]
dmel <- args.env[["-dmel"]]
# if dmel is TRUE, arms 2L/2R and 3L/3R will be considered one chromosome during recombination

# Do not enter these two arguments unless you are generating recombinant inbred lines
n_RILs <- as.numeric(args.env[["-nRILs"]])
inbreed_generations <- as.numeric(args.env[["-inbreed_generations"]])

# Optional seed override
manual_seed <- args.env[["-seed"]]

# Set seed; use L'Ecuyer-CMRG seed for reproducible %dopar% results
if(! is.null(manual_seed)) {
  cat("manual seed: ")
  cat(manual_seed)
  cat("\n")
  set.seed(manual_seed,  kind = "L'Ecuyer-CMRG")
} else {
  # generate random seed by converting the combined input string to md5,
  # subsetting the first six characters,
  # and converting from hexadecimal to integer.
  input_string <- paste(bed_file, prefix, popSize, growthRate, sexModel, n_founders, n_generations, lineIDs_filename, x.chromosome, recombinationModel, n_threads, iteration, dmel, n_RILs, inbreed_generations, sep=" ")
  cat("input string: ")
  cat(input_string)
  starting_seed <-  strtoi(substr(digest(input_string, algo="md5"), 1, 7), base=16L)
  cat("\nstarting seed: ")
  cat(starting_seed)
  cat("\n")
  set.seed(starting_seed,  kind = "L'Ecuyer-CMRG")
}

filestem <- paste(prefix, "_", n_founders, "_F", n_generations, "_", iteration, sep="")


# If n_RILs argument exists, set flag
if(! length(n_RILs)==0) {make_RILs <- TRUE} else {make_RILs <- FALSE}


if(n_threads > 1) {
    library(doMC)
    registerDoMC(cores=n_threads)
}

initializePopulation <- function(bed, n_founders, sexModel) {
    DT <- foreach(chromosome=unique(bed$chr), .combine="rbind", .inorder=TRUE) %do% {
        data.table(chromosome, haplotype=1:2, start=1, stop=max(bed[chr==chromosome]$stop))
    }
    DT <- DT[rep(1:.N, each=n_founders)]
    DT[,founderID:=rep(1:n_founders, length(chromosomes)*2)]
    if(sexModel=="dioecious") {
            DT <- DT[rep(1:.N, each=2)]
            DT[, sex := rep(c("M","F"), n_founders*length(chromosomes)*2)]
    } else if(sexModel=="hermaphroditic") {
        DT[, sex := "H"]
    }
    setkey(DT, founderID, sex, chromosome, haplotype, start, stop)
    DT[, ind := rleid(sex,founderID)]

    # remove 2nd X chromosome from males
    return(DT[! (sex=="M" & haplotype==2 & chromosome=="X")][])
}


getGamete <- function(pop, ind.i, haplotype.i, parSex, recombinant) {
    if(recombinant==FALSE) {
        # simply return random selection of alleles
        foreach(chromosome.i=chromosomes, .combine="rbind", .inorder=TRUE) %do% {

            dt.out <- pop[.(ind.i, chromosome.i, sample(1:2, size=1))]
            dt.out[, haplotype := haplotype.i]
            return(dt.out[,c("chromosome","haplotype","start","stop", "founderID")])
                }

    } else if(recombinant==TRUE) {
        # pull out individual and recombine

        foreach(chromosome.i = chromosomes, .combine="rbind", .inorder=TRUE) %do% {
            if(chromosome.i=="X" & parSex=="M") {
                breakpoints <- breakpoints <- c(0, chromosome_sizes[[chromosome.i]])
            } else {
                breakpoints <- unique(trunc(sort(recombination_function[[chromosome.i]](runif(rpois(1, lambda=recombination_rates[[chromosome.i]]))))))
                breakpoints <- c(0, breakpoints, chromosome_sizes[[chromosome.i]])
            }

            N.recomb <- length(breakpoints) - 2
            startHaplotype <- sample(c(1,2), size=1)


            ranges <- data.table(
                "ind"=ind.i,
                "chromosome"=chromosome.i,
                "haplotype"=(((1 + startHaplotype) : (N.recomb + 1 + startHaplotype)) %% 2 + 1),
                "start"=1+breakpoints[1:(length(breakpoints)-1)],
                "stop"=breakpoints[2:length(breakpoints)])


            foreach(ind.i=ranges$ind, chromosome.i=ranges$chromosome, haplotype.j=ranges$haplotype,
                    start.i=ranges$start, stop.i=ranges$stop,.combine="rbind", .inorder=TRUE) %do% {
                    if(chromosome.i=="X" & parSex=="M") {
                        dt.out <- pop[.(ind.i, chromosome.i, haplotype.j)]
                        dt.out[,haplotype := haplotype.i]
                        return(dt.out[,c("chromosome","haplotype","start","stop", "founderID")])

                    } else {
                    dt.out <- pop[.(ind.i, chromosome.i, haplotype.j)][! (stop < start.i)  & ! (start > stop.i)]
                    dt.out[1,start := start.i]
                    dt.out[dim(dt.out)[1], stop := stop.i]
                    dt.out[,haplotype := haplotype.i]
                    return(dt.out[,c("chromosome","haplotype","start","stop", "founderID")])
                    }
            }
        }
    }
}

doSex <- function(pop, sexModel, n) {

    N.inds <- unique(pop$ind)
    if(sexModel=="hermaphroditic") {
        par1 <- sample(N.inds, size=1)
        par2 <- sample(N.inds, size=1)
    }else if(sexModel=="dioecious") {
        # pick first gamete
        par1 <- sample(N.inds, size=1)
        par1sex <- unique(pop[ind==par1]$sex)

        # pick second gamete from opposite sex
        if(length(N.inds)==2) {
            par2 <- unique(pop[ind!=par1]$ind)}
        else {
            par2 <- sample(unique(pop[sex!=par1sex]$ind), size=1)
        }
        par2sex <- unique(pop[ind==par2]$sex)
    }
    g1 <- getGamete(pop, par1, 1, par1sex, ifelse(recombinationModel=="femaleOnly" & par1sex=="M", FALSE, TRUE))
    g2 <- getGamete(pop, par2, 2, par2sex, ifelse(recombinationModel=="femaleOnly" & par2sex=="M", FALSE, TRUE))
    ind.out <- rbindlist(list(g1,g2))[!is.na(founderID)]
    end.out[,ind := n]
    if (length(unique(ind.out[chromosome=="X"]$haplotype))==1) {
        ind.out[, sex := "M"]
    } else {
        offspringSex <- ifelse(sexModel=="dioecious", "F", "H")
        ind.out[, sex := offspringSex]
    }
    return(ind.out[])
}

# Function for drosophila corrections, e.g. combining/splitting 2L/2R as a single chromosome
translateLandR <- function(DT, max2L=23100000L, max3L=24600000L) {
    splitting_chr2 <- DT[chromosome=="2" & max2L > start & max2L < stop]
    splitting_chr3 <- DT[chromosome=="3" & max3L > start & max3L < stop]

    # Restructure regions overlapping centromere on chromosome 2
    chr2L <- copy(splitting_chr2)
    chr2R <- copy(splitting_chr2)
    chr2L[, stop := max2L]
    chr2L[, chromosome := "2L"]
    chr2R[, start := 1]
    chr2R[, stop := stop - max2L]
    chr2R[, chromosome := "2R"]


    # Restructure regions overlapping centromere on chromosome 3
    chr3L <- copy(splitting_chr3)
    chr3R <- copy(splitting_chr3)
    chr3L[, stop := max3L]
    chr3L[, chromosome := "3L"]
    chr3R[, start := 1]
    chr3R[, stop := stop - max3L]
    chr3R[, chromosome := "3R"]

    # Restructure regions not overlapping centromere on chromosome 2
    dat.2Lor2R <- copy(DT[chromosome=="2"][stop < max2L | start > max2L])
    dat.2Lor2R[stop < max2L, chromosome := "2L"]
    dat.2Lor2R[start > max2L, chromosome := "2R"]
    dat.2Lor2R[chromosome=="2R", start := start - max2L]
    dat.2Lor2R[chromosome=="2R", stop := stop - max2L]

    # Restructure regions not overlapping centromere on chromosome 3
    dat.3Lor3R <- copy(DT[chromosome=="3"][stop < max3L | start > max3L])
    dat.3Lor3R[stop < max3L, chromosome := "3L"]
    dat.3Lor3R[start > max3L, chromosome := "3R"]
    dat.3Lor3R[chromosome=="3R", start := start - max3L]
    dat.3Lor3R[chromosome=="3R", stop := stop - max3L]

    # Keep chromosome X as is
    dat.X <- DT[chromosome=="X"]

    # Combine ranges
    dat.all <- rbindlist(list(
                dat.2Lor2R,
                dat.3Lor3R,
                chr2L,
                chr2R,
                chr3L,
                chr3R,
                dat.X
                ))

    # Reorder
    setkey(dat.all, ind, chromosome, haplotype, start)

    return(dat.all)
}

# Load bed file
bed <- fread(bed_file)
setnames(bed, c("chr", "start", "stop", "c"))
bed[chr==x.chromosome, chr := "X"]

# Correction for Drosophila
if(dmel==TRUE) {
    # stoare & add maximum value of 2L onto every start, stop for 2R
    # store & add maximum value of 3L onto every star,t stop for 3R
    # Reduce these later
    max2L <- as.integer(max(bed[chr=="2L"]$stop))
    max3L <- as.integer(max(bed[chr=="3L"]$stop))
    bed[chr=="2R", start := start + max2L]
    bed[chr=="2R", stop := stop + max2L]
    bed[chr=="3R", start := start + max3L]
    bed[chr=="3R", stop := stop + max3L]
    bed[chr %in% c("2L","2R"), chr := "2"]
    bed[chr %in% c("3L","3R"), chr := "3"]
}

# Get list of unique chromosome names within .bed file
chromosomes <- unique(bed$chr)

# Convert c (cM per Mb) to Morgans
bed[, M := c * ((stop-start)/1e8)]

# Create hash table with chr -> expected value for number of recombination events
#   e.g.,
#   > recombination_rates[["2L"]]
#   [1] 0.5533038
recombination_rates <- new.env()
for(chromosome in chromosomes) {
    recombination_rates[[chromosome]] <- sum(bed[chr==chromosome]$M)       # convert c (cM per Megabase) to Morgans
}

chromosome_sizes <- new.env()
for(chromosome in chromosomes) {
    chromosome_sizes[[chromosome]] <- max(bed[chr==chromosome]$stop)
}


# Create hash table with random value (0,1) -> recombination position, via linear interpolation of scaled cumulative sum of recombination rates
bed[, cumulative_M := cumsum(M), by=chr]
bed[, scaled := cumulative_M/max(cumulative_M), by=chr]

genomeSize <- sum(bed[, list(size=max(stop)), by=chr]$size)

recombination_function <- new.env()
for(chromosome in chromosomes) {
    recombination_function[[as.character(chromosome)]] <- approxfun(c(0, bed[chr==chromosome]$scaled), c(0,bed[chr==chromosome]$stop))
}

# Extract header column from (possibly zgipped) vcf file
lineIDs <- readLines(lineIDs_filename)

# Subset (size=N founders) founder IDs, in order that they appear in the header
used_IDs <- lineIDs[sort(sample(length(lineIDs), size=n_founders, replace=FALSE))]


# Initialize Population

pop <- initializePopulation(bed, n_founders, sexModel)
setkey(pop, ind, chromosome, haplotype, start)


# Iterate through generations

for(i in 1:n_generations) {
    print(i)
    pop2 <- foreach(n=1:trunc((popSize*growthRate**i)), .combine="rbind", .inorder=TRUE) %dopar% {
        doSex(pop, sexModel, n)
    }
    setkey(pop2, ind, chromosome, haplotype, start)
    #write.table(pop2, file=paste(stem, "_", i, ".txt", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    popSize <- trunc(popSize*growthRate)
    pop <- copy(pop2)
}

pop.outbred <- copy(pop2)

# Steps unique to generating RILs:
if(make_RILs == TRUE) {
        RILs <- foreach(nRIL=1:n_RILs, .combine="rbind", .errorhandling="remove", .inorder=TRUE) %do% {
        print(paste("RIL ", nRIL, sep=""))
        M.founder <- sample(unique(pop.outbred[sex=="M"][,ind]), size=1)
        F.founder <- sample(unique(pop.outbred[sex=="F"][,ind]), size=1)
        pop <- pop.outbred[ind %in% c(M.founder, F.founder)]

    # Do N many generations of inbreeding
        for(i in 1:inbreed_generations) {
        # Select a single male and female
        M.ind <- unique(pop[sex=="M"][,ind])
        F.ind <- unique(pop[sex=="F"][,ind])
        # Subset population to that single male and single female
        pop <- pop[ind %in% c(M.ind, F.ind)]

        # Generate first individual
        ind1 <- doSex(pop, sexModel, 1)
        ind1sex <- unique(ind1[,sex])

        # Generate second individual of different sex
        while(TRUE) {
            ind2 <- doSex(pop, sexModel, 2)
            if(unique(ind2[,sex]) != ind1sex) {
                break
            }
        }
        pop <- rbindlist(list(ind1, ind2))
        setkey(pop, ind, chromosome, haplotype, start)
    }

    # Return one single inbred individual (female)
    dt.return <- pop[sex=="F"]
    dt.return[, ind := nRIL][]
    return(dt.return)
    }

    # Fully homozygose RILs by forcing 2nd haplotype to be identical to the 1st
    RILs_2 <- copy(RILs[haplotype==1])
    RILs_2[, haplotype := 2]
    pop <- rbindlist(list(RILs[haplotype==1], RILs_2))

    # set n_founders to n_RILs for filename output
    n_founders <- n_RILs

    pop <- RILs

    # Translate founderID to lineID
    pop[, lineID := used_IDs[pop$founderID]]

    # Remove numeric founderID column
    pop[, founderID := NULL]

    # change X chromosome back to chrN
    pop[chromosome=="X", chromosome := x.chromosome]

    # Corrections for Drosophila (splits "2" into "2L" + "2R"; "3" into "3L" + "3R")
    if(dmel == TRUE) { pop <- translateLandR(pop) }
    setnames(pop, "ind", "RIL_ID")

    # Sort haplotype map
    setkey(pop, RIL_ID, chromosome, haplotype, start)

} else {
    # change X chromosome back to chrN
    pop[chromosome=="X", chromosome := x.chromosome]

    # Translate founderID to lineID
    pop[, lineID := used_IDs[pop$founderID]]

    # Corrections for Drosophila (splits "2" into "2L" + "2R"; "3" into "3L" + "3R")
    if(dmel == TRUE) { pop <- translateLandR(pop) }
    setkey(pop, ind, chromosome, haplotype, start)
}

# Write to file
    options(scipen=999)
    write.table(pop, file=paste(filestem, ".haps", sep=""), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    # write.table(used_IDs, file=paste(iteration, ".founders.txt", sep=""), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
