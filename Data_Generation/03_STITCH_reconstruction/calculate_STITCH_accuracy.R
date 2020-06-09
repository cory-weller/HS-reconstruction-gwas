#!/usr/bin/env Rscript

if(! require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if(! require(data.table)) {
  install.packages("data.table")
  library(data.table)
}

if(! require(foreach)) {
  install.packages("foreach")
  library(foreach)
}

if(! require(doMC)) {
  install.packages("doMC")
  library(doMC)
}
registerDoMC(cores=10)


if(! require(SNPRelate)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("SNPRelate")
  }
}

registerDoMC(cores=10)

args <- commandArgs(trailingOnly=TRUE)

population <- args[1]
population_size <- args[2]
n_generations <- unlist(strsplit(population, "_"))[4]
chromosome <- "2L"


# population <- "hybrid_swarm_32_F50"
# population_size <- "100"

# convert stitch VCF to GDS
vcf.filename <- paste(population, "_", population_size, "/", "stitch.2L.vcf.gz", sep="")
population.vcf.filename <- paste("../04_RABBIT/", population, "/", chromosome, ".RABBIT.vcf", sep="")
gds.filename <- paste(population, "_", population_size,"/", "stitch.2L.vcf.gds", sep="")
accuracy.filename <- paste(population, "_", population_size,"/", "accuracy.RDS", sep="")

if(! file.exists(gds.filename)) {
  cat("Building gds...\n")
  snpgdsVCF2GDS(vcf.fn = vcf.filename, out.fn = gds.filename)
} else {
  cat("gds already exists...\n")
}

# load fixed sites
population.vcf <- fread(population.vcf.filename)
chosen.founders <- colnames(population.vcf)[colnames(population.vcf) %like% "RAL"]
population.vcf[, nRef := apply(.SD, 1, function(x) sum(x == 0, na.rm=TRUE)), .SDcols=chosen.founders]
population.vcf[, nAlt := apply(.SD, 1, function(x) sum(x == 1, na.rm=TRUE)), .SDcols=chosen.founders]
population.vcf[, refFreq := nRef/(nRef+nAlt)]
population.vcf[, altFreq := nAlt/(nRef+nAlt)]
variable.sites <- population.vcf[refFreq != 0 & altFreq != 0][,POS]

if(! file.exists(accuracy.filename)) {
  cat("Opening gds\n")
  gds <- snpgdsOpen(gds.filename, allow.fork=TRUE)
  
  
  gds.dt <- data.table(id=read.gdsn(index.gdsn(gds, path="snp.id")),
                          chr=read.gdsn(index.gdsn(gds, path="snp.chromosome")),
                          pos=read.gdsn(index.gdsn(gds, path="snp.position")),
                          allele=read.gdsn(index.gdsn(gds, path="snp.allele"))
                      )
  
  N_sites <- dim(gds.dt)[[1]]
  N_individuals <- length(read.gdsn(index.gdsn(gds, path="sample.id")))
  
  if (N_individuals != population_size) {
    cat("gds sample count differs from population size! Exiting...\n")
    exit(1)
  }
  
  cat("Calculating STITCH genotype accuracy\n")
  o <- foreach(individual_N = 1:N_individuals, .combine="rbind") %dopar% {
    cat(paste(individual_N, "\n", sep=""))
    # read in true genotype vcf file
    truth.dt <- fread(cmd=paste("grep '^2L' ../02_simulate_reads/", population, "/", individual_N, ".geno", sep=""))
    setnames(truth.dt, c("chr", "pos", "A1", "A2"))
    setkey(truth.dt, "chr", "pos")
  # convert genotypes and calculate reference allele dosage
    # 3 = missing data
    # 2 = ref/ref
    # 1 = heterozygote
    # 0 = alt/alt
    truth.dt[, true_genotype := 3]
    truth.dt[A1=="0/0" & A2 == "0/0", true_genotype := 2]
    truth.dt[A1=="1/1" & A2 == "1/1", true_genotype := 0]
    truth.dt[A1=="0/0" & A2 == "1/1", true_genotype := 1]
    truth.dt[A1=="1/1" & A2 == "0/0", true_genotype := 1]
    

    gds.dt[, stitch_genotype := read.gdsn(index.gdsn(gds, "genotype"), start=c(individual_N,1), count=c(1, N_sites))]
    dat.merged <- merge(truth.dt, gds.dt, by=c("chr","pos"))[,c("chr","pos","true_genotype","stitch_genotype")]
    
    # remove fixed sites
    dat.merged <- dat.merged[pos %in% variable.sites]

    # get summary stats
    # numerator: # of sites with accurate stitch genotype (excluding unknowns)
      # i.e. [true_genotype != 3 & stitch_genotype != 3 & true_genotype == stitch_genotype]
    # denominator: # of sites where true and estimated genotypes have known values
      # i.e. [true_genotype != 3 & stitch_genotype != 3]

    accurate_sites <- nrow(dat.merged [true_genotype != 3 & stitch_genotype != 3 & true_genotype == stitch_genotype])
    total_sites <- nrow(dat.merged[true_genotype != 3 & stitch_genotype != 3])

    data.table( "population_size" = population_size,
                "individual_N" = individual_N,
                "accurate_sites" = accurate_sites,
                "total_sites" = total_sites)
  }
  o[, percent_correct := accurate_sites / total_sites]
  o[, "N" := population_size]
  o[, "n_generations" := n_generations]
  saveRDS(o, file=accuracy.filename)
} else {
  cat("accuracy already calculated...\n")
  cat("loading accuracy.RDS\n")
  o <- readRDS(accuracy.filename)
}
