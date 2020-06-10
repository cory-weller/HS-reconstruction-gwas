#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)
library(foreach)
library(doMC)
registerDoMC(cores=5)

args <- commandArgs(trailingOnly=TRUE)
population <- args[1]
replicates <- as.numeric(args[2])
n_loci_subset <- as.numeric(args[3])
effect_size <- as.numeric(args[4])

# Change to TRUE if it is desired to sample hybrid swarm haplotypes with replacement.
# Ensure that at least N=5000 haplotypes are available for sampling, 
# or reduce number of individuals sampled for the GWAS mapping population,
# otherwise the simulation will fail.
sampleWithReplacement <- FALSE
mappingPopulationSize <- 5000


nFoundersSplit <- unlist(strsplit(unlist(strsplit(population, split="_F"))[1], split="_"))
nFounders <- as.numeric(nFoundersSplit[length(nFoundersSplit)])

print(replicates)
print(population)
print(effect_size)

if(grepl("F0$", population) == TRUE) {
  inbred <- TRUE
} else {
  inbred <- FALSE 
}

if(grepl("^RIL", population) == TRUE) {
  RIL <- TRUE
} else {
  RIL <- FALSE 
}

source("../liability.functions.R")

#all_sites <- make_sites_vcf_file(autosomes)
all_sites <- make_sites_vcf_file("2L")
n_loci <- 1000

# get total list of sites
sites <- all_sites[sample(.N, size=n_loci, replace=FALSE)][,c("CHROM","POS")]

# retrieve tabix info from sites
causative <- tabix_read(sites)
setnames(causative, "#CHROM", "CHROM")
  


# From here, get individual haplotype maps and convert to risk, then assign to case/control

if(inbred == FALSE) {
  haps <- load_haplotypes(population, RIL)
  origFounderIDs <- readLines(paste("../../01_forward_simulator/", population, ".founders", sep=""))
} else if(inbred == TRUE) {
  founding_lines <- readLines("../../01_forward_simulator/lines.txt")
  origFounderIDs <- sample(founding_lines, size=nFounders, replace=FALSE)
  haps <- data.table("start"= "1" ,
                      "stop"= "25000000" ,
                      "par1"= sample(origFounderIDs, size=10000, replace=TRUE)
                    )
  haps[, par2 := par1]
  haps[, chromosome := "2L"]
  haps[, ind := 1:10000]
}

if(RIL==TRUE) {
# replicate draws   
inds <- sample(unique(haps[,ind]), size=mappingPopulationSize, replace=TRUE)
}

vcf <- load_VCF("2L")


allFounderIDs <- colnames(causative)[colnames(causative) %like%  "RAL"]

replicate <- 1
while(replicate <= replicates) {
tryCatch({
#for(replicate in 1:replicates) {

permuted_haplotypes <- permuteHaplotypes(haps, origFounderIDs, allFounderIDs)
used_founders <- unique(c(permuted_haplotypes[,par1], permuted_haplotypes[,par2]))

causative_subset <- causative[, c("CHROM","POS",used_founders), with=FALSE]

# calculate allele freqs
causative_subset[, nRef := apply(.SD, 1, function(x) sum(x=="0/0", na.rm=TRUE)), .SDcols=colnames(causative_subset)[colnames(causative_subset) %like% "RAL"] ]
causative_subset[, nAlt := apply(.SD, 1, function(x) sum(x=="1/1", na.rm=TRUE)), .SDcols=colnames(causative_subset)[colnames(causative_subset) %like% "RAL"] ]
causative_subset[, nTot := nRef + nAlt]
causative_subset[, refFreq := nRef/nTot]
causative_subset[, altFreq := nAlt/nTot]
causative_subset[, effect := sapply(altFreq, function(x) sle(x, effect_size))]
causative_subset <- causative_subset[nRef != 0 & nAlt != 0]
setkey(causative_subset, CHROM, POS)

cat("calculating dosage at causative loci\n")
if(effect_size != 0) {
  dosage <- get_dosage_at_causative_loci(permuted_haplotypes, causative_subset) 
  dosage[, altDosage := apply(.SD, 1, function(x) sum(x=="1/1", na.rm=T)), .SDcols=c("hap1","hap2")]
  causative_subset <- causative_subset[POS %in% dosage[,POS]]

}
setkey(permuted_haplotypes, ind)

  cat(paste("rep ", replicate, "\n", sep=""))
  # Subset causaative loci
  if(effect_size != 0) {

  causative_subset <- causative_subset[sample(.N, size=n_loci_subset, replace=FALSE)]
  dosage_subset <- merge(causative_subset[,c("CHROM","POS","effect")], dosage, by=c("CHROM","POS"))
  dosage_subset <- dosage_subset[, list("liability"=sum(effect*altDosage)), by=list(ind_n)]
  dosage_subset[, risk := liability(liability)]
  if(RIL==TRUE) {
    chosen <- dosage_subset[sample(800, replace=sampleWithReplacement, size=mappingPopulationSize)]
    chosen[, rnd := runif(.N, min=0, max=1)]
    chosen[, group := ifelse(rnd <= risk, "case", "control")]
  } else {
  # if not RIL with real PVE  
    dosage_subset[, rnd := runif(.N, min=0, max=1)]
    dosage_subset[, group := ifelse(rnd <= risk, "case", "control")]
    chosen <- dosage_subset[ind_n %in% sample(10000, size=mappingPopulationSize, replace=FALSE)]
  }

} else {
  # when PVE == 0
    chosen <- data.table("ind_n"=1:10000)
    chosen[, "group" := sample(c("case","control"), size=10000, replace=TRUE)]
    chosen <- chosen[sample(.N, size=mappingPopulationSize, replace=FALSE)]
}
  cat("working on case population freqs\n")
  caseFreqs <- calcFreqs(permuted_haplotypes[.(chosen[group=="case",ind_n]), allow.cartesian=TRUE], vcf)
  cat("working on control population freqs\n")
  controlFreqs <- calcFreqs(permuted_haplotypes[.(chosen[group=="control",ind_n]), allow.cartesian=TRUE], vcf)
  cat("tabulating frequencies\n")
  setnames(caseFreqs, c("CHROM","POS","case.Ref", "case.Alt"))
  setnames(controlFreqs, c("CHROM","POS","control.Ref", "control.Alt"))
  # Merge case and control tables together
  setkey(caseFreqs, CHROM, POS)
  setkey(controlFreqs, CHROM, POS)
  freqs.merge <- merge(caseFreqs, controlFreqs)
  # Exclude rows that contain any zeros, which produce NaN in chisq test
  freqs.merge <- freqs.merge[case.Ref != 0 & case.Alt != 0 & control.Ref != 0 & control.Alt != 0]
  # Perform CHISQ test
  freqs.merge[,chisq := (((case.Ref)+(control.Ref)+((case.Alt))+((control.Alt))) * (((control.Ref)*(case.Alt)) - ((case.Ref)*(control.Alt)))**2) / (((case.Ref)+(control.Ref)) * (((case.Ref)+(case.Alt))) * ((control.Ref)+(control.Alt)) * (((case.Alt))+(control.Alt)))]
  if(inbred==TRUE | RIL == TRUE) {
    # perform correction for sample size of homozygote allele draws
    freqs.merge[, chisq := chisq/2]
  }
  freqs.merge[,P := pchisq(chisq, 1, lower.tail=FALSE)]
  # Convert to character and back to numeric for floating point precision rounding reasons
  freqs.merge[, P := as.numeric(as.character(P)) ]

  if(effect_size == 0) {
    N_case <- nrow(chosen[group=="case"])
    N_control <- nrow(chosen[group=="control"])
    min_P <- min(freqs.merge[,P])
    output <- data.table("population"=population,
                          "min_P" = min_P
                        )
    fwrite(output, file="permutations.log", quote=F, append=T, row.names=F, col.names=F, sep="\t")
    
  } else {
      freqs.fn <- paste(replicate, ".freqs.txt", sep="")
      causative.fn <- paste(replicate, ".causative.txt", sep="")
      fwrite(freqs.merge, file=freqs.fn, row.names=F, col.names=T, sep="\t", quote=FALSE)
      fwrite(causative_subset[,c("CHROM","POS","refFreq","altFreq","effect")], file=causative.fn, row.names=F, col.names=T, sep="\t", quote=FALSE)
  }
  replicate <- replicate + 1
}, error=function(e){})
}
