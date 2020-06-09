#!/usr/bin/env/Rscript

library(pROC)
library(data.table)
library(ggplot2)
#library(cowplot)
library(foreach)
library(doMC)
library(ggthemes)
registerDoMC(cores=5)


#args <- commandArgs(trailingOnly=TRUE)
# zipFileName <- args[1]


## FUNCTIONS

getAuc <- function(count, population, n_loci, effect_size, replicate, distance=0) {
  filestem <- paste(population, "_", count, "_", n_loci, "_", effect_size, sep="")
  zip_fn <- paste(filestem, ".zip", sep="")
  
  freqs_fn <- paste(replicate, ".freqs.txt", sep="")
  causative_fn <- paste(replicate, ".causative.txt", sep="")
  
  freqs_cmd <- paste("unzip -p ", zip_fn, " ", freqs_fn, sep="")
  causative_cmd <- paste("unzip -p ", zip_fn, " ", causative_fn, sep="")
  
  freqs <- fread(cmd=freqs_cmd)
  causative <- fread(cmd=causative_cmd)
  
  setkey(freqs, CHROM, POS)
  setkey(causative, CHROM, POS)

  # To Do: add distance metric

  # causative / case SNP = 1, all others are controls = 0
  freqs[, causative_snp := 0]
  freqs[causative, causative_snp := 1]
  freqs[, minuslogp := -1*log10(P)]
  
  # weight SNPs by effect size?
  sumEffect <- sum(abs(causative$effect))
  causative[, weight := abs(effect)/sumEffect]
  causative[, N := trunc(weight/min(weight))]
  
  freqs[, idx := 1:.N]
  freqs[causative]
  idxs <- unlist(mapply(rep, x = freqs[causative][,idx], times = freqs[causative][,N]))
  freqs[idxs]
  freqs2 <- rbindlist(list(freqs, freqs[idxs]))
  res <- roc(freqs, response=causative_snp, predictor=minuslogp, direction="<")
  #res2 <- roc(freqs2, response=causative_snp, predictor=minuslogp)
  #g1 <- ggroc(res1) + geom_abline(slope=1, intercept=1, linetype="dashed", alpha=0.5) + labs(title=paste("Unscaled\n", capture.output(res1$auc), sep=""))
  #g2 <- ggroc(res2) + geom_abline(slope=1, intercept=1, linetype="dashed", alpha=0.5) + labs(title=paste("Y-contribution scaled to effect size\n", capture.output(res2$auc), sep=""))

  #g.all <- plot_grid(g1, g2, labels = c("A", "B"), align = "h", nrow=1)
  #plot(g.all)
  #print(freqs[causative][order(-minuslogp)])
  dat <- data.table(x1=abs(res$specificities - 1), y=res$sensitivities)
  dat.ag <- dat[dat[, .I[x1==min(x1)], by=y]$V1]
  dat.ag[, x2 := shift(x1, type="lag", n=1L, fill=1)]
  dat.ag[, "replicate" := replicate]
  dat.ag[, "count" := count]
  dat.ag[, "effect_size" := effect_size]
  dat.ag[, "n_loci" := n_loci]
  dat.ag[, "population" := population]
  
  return(dat.ag)
}

get_pAUC <- function(dt, pAUC_threshold) {
  dt1 <- dt[x1 < pAUC_threshold]
  dt1[which.max(x1), x2 := pAUC_threshold]
  sum(dt1[, list(y*(x2-x1))][, V1])
}

## RUN

for (zipFileName in list.files(pattern="*.zip")) {
  filestem = unlist(strsplit(zipFileName, split=".zip"))[1]
  RDSfilename <- paste(filestem, ".RDS", sep="")
  # if RDS not yet calculated, do so
  if(! file.exists(RDSfilename)) {
    print(filestem)
    zipFileNameSplit <- unlist(strsplit(zipFileName, split="_"))
    L <- length(zipFileNameSplit)
    effect_size <- unlist(strsplit(zipFileNameSplit[L], split=".zip"))[1]
    n_loci <- zipFileNameSplit[L-1]
    count <- zipFileNameSplit[L-2]
    population <- paste(zipFileNameSplit[1:(L-3)], collapse="_")
    dat <- foreach(replicate=1:50, .combine="rbind", .errorhandling="remove") %dopar% {
      getAuc(count, population, n_loci, effect_size, replicate)
    }
    saveRDS(dat, file=RDSfilename)
  }
}

# To convert AUC curve table to compressed form:

# 
not_run <- function() {
  dat <- rbindlist(lapply(list.files(pattern="*[0-9].RDS"), function(x) readRDS(x)))
  
  o <- foreach(x=seq(0.005,1,0.005), .combine="rbind") %do% {
    dat[x1 <= x & x2 >= x, list("x"=x, "y"=mean(y)), by=list(effect_size, n_loci, population)]
  }
  
  o[, c("n_founders","n_generations") := tstrsplit(population, "_")[3:4]]
  o[population %like% "outbred_", c("n_founders","n_generations") := tstrsplit(population, "_")[2:3]]
  o[population %like% "RILs_", "n_founders" := "8" ]
  o[population %like% "RILs_", "n_generations" := "F50" ]


  o[, n_loci := as.numeric(n_loci)]
  
  o[n_founders==32 & n_generations=="F5", label := "HS_32_F5"]
  o[n_founders==128 & n_generations=="F5", label := "HS_128_F5"]
  o[n_founders==128 & n_generations=="F1", label := "HS_128_F1"]
  o[n_founders==128 & n_generations=="F2", label := "HS_128_F2"]
  o[n_founders==32 & n_generations=="F1", label := "HS_32_F1"]
  o[n_founders==32 & n_generations=="F2", label := "HS_32_F2"]
  o[n_founders==8 & n_generations=="F50", label := "RILs (DSPR)"]
  o[n_founders==128 & n_generations=="F50", label := "Outbred_128_F50"]
  o[n_founders==32 & n_generations=="F50", label := "Outbred_32_F50"]
  o[n_founders==128 & n_generations=="F0", label := "ILs (DGRP)"]
  o[n_founders==128 & n_generations=="F0", label := "ILs (DGRP 128)"]
  o[n_founders==32 & n_generations=="F0", label := "ILs (DGRP (32))"]
  

  # Plot 3 types comparison
  ggplot(o[n_generations != "F1" & n_generations != "F2" ], aes(x=x, y=y, color=label)) + geom_line() + 
  geom_abline(slope=1, intercept=0) + xlim(0,1) + ylim(0,1) + facet_grid(n_loci ~ effect_size, labeller="label_both") +
  labs(x="1 - Specificity", y="Sensitivity") +
  theme_few(14)
  
    ggplot(o[n_founders==128 & n_generations %in% c("F0","F1","F2","F5","F50") ], aes(x=x, y=y, color=type, linetype=n_generations)) + geom_line() + 
    geom_abline(slope=1, intercept=0) + xlim(0,1) + ylim(0,1) + facet_grid(n_loci ~ effect_size, labeller="label_both") +
    labs(x="1 - Specificity", y="Sensitivity") +
    theme_few(14)
    
    ggplot(o[n_founders==128 & n_generations %in% c("F0","F1","F2","F5","F50") ], aes(x=x, y=y, color=type, linetype=n_generations)) + geom_line() + 
    geom_abline(slope=1, intercept=0) + xlim(0,1) + ylim(0,1) + facet_grid(n_loci ~ effect_size, labeller="label_both") +
    labs(x="1 - Specificity", y="Sensitivity") +
    theme_few(14)
    
    ggplot(o[!is.na(label)], aes(x=x, y=y, color=label)) + geom_line() + 
    geom_abline(slope=1, intercept=0) + xlim(0,1) + ylim(0,1) + facet_grid(n_loci ~ effect_size, labeller="label_both") +
    labs(x="1 - Specificity", y="Sensitivity") +
    theme_few(14)
  
  # Plot outbred comparison
  
  
  ggplot(o, aes(x=x, y=y, color=effect_size, linetype=n_loci)) + geom_line() + 
  geom_abline(slope=1, intercept=0) + xlim(0,1) + ylim(0,1) + facet_wrap(population~.) +
  labs(x="1 - Specificity", y="Sensitivity")
  
  
  dat[, grp := .GRP, by=list(replicate, count, effect_size, n_loci, population)]

  dat[, get_pAUC(.SD, 0.5), by=grp]

  pAUCs <- foreach(pAUC=1, .combine="rbind") %do% {
    dat[, .("pAUC_threshold"=pAUC, "pAUC"=get_pAUC(.SD, pAUC)), by=.(replicate, count, effect_size, n_loci, population)]
  }

  saveRDS(pAUCs, file="pAUCs.RDS")


  foreach(pAUC=c(0.2, 0.5), .combine="rbind") %do% {
    dat[grp<10][, .("pAUC_threshold"=pAUC, "pAUC"=get_pAUC(.SD, pAUC)), by=.(replicate, count, effect_size, n_loci, population)]
  }
  
  # Permutations of median pAUC difference?

  dt[sample(.N, size=.N/2, replace=FALSE)]

  o <- foreach(i=1:1000, .combine="rbind") %dopar% {
    g1 <- sample(1000, size=500, replace=FALSE)
    g2 <- setdiff(1:1000, g1)

    med1 <- median(dt[g1, pAUC])
    mean1 <- mean(dt[g1, pAUC])
    med2 <- median(dt[g2, pAUC])
    mean2 <- mean(dt[g2, pAUC])
    deltaMed <- max(abs(med2-med1))
    deltaMean <- max(abs(mean2-mean1))

    data.table("iteration"=i, deltaMed, deltaMean)
  }

  o[, deltamed := med2-med1]
  o[, deltamean := mean2-mean1]

  dt[, list(median(pAUC)), by=population]
}
# 




# 
# 
# 
#



# fwrite(dat, file="test_condensed.dat", quote=F, row.names=F, col.names=T, sep="\t")







# f0_5_0.2[, n_loci := 5]
# f0_5_0.2[, nGenerations := "F0"]
# f5_5_0.2[, n_loci := 5]
# f5_5_0.2[, nGenerations := "F5"]
# 
# dat <- rbindlist(list(f5_5_0.2, f0_5_0.2))
# dat[, effect_size := 0.2]
# 
# dat.ag <- dat[,list("mean_x"=mean(mean_x), "mean_y"=mean(mean_y), "sd_y"=sd(mean_y)), by=list(grp, n_loci, effect_size, nGenerations)]
# dat.ag[, bottom_y := mean_y - sd_y]
# dat.ag[, top_y := mean_y + sd_y]
# dat.ag[, type := paste(n_loci, effect_size, sep="_")]
# 
# ggplot(dat.ag, aes(x=mean_x, y=mean_y, color=nGenerations, group=nGenerations, ymin=bottom_y, ymax=top_y)) +
# geom_line() + scale_x_reverse() +
# labs(x="1 - Specificity", y="Sensitivity", title="100 replicates each from F5, 5 loci (low), 0.2 effect_size (low)") +
# ylim(0,1)
# 
# 
# 
# 
# 
# 
# 
# 
# f5_5_0.4 <- foreach(i=1:50, .combine="rbind") %dopar% {
#   foreach(j=1:2, .combine="rbind") %dopar% {
#     getAuc(j, "hybrid_swarm_32_F5", 5, 0.4, i)
#   }
# }
# 
# f5_10_0.2 <- foreach(i=1:50, .combine="rbind") %dopar% {
#   foreach(j=1:2, .combine="rbind") %dopar% {
#     getAuc(j, "hybrid_swarm_32_F5", 10, 0.2, i)
#   }
# }
# 
# f5_5_0.2 <- fread('f5_means.txt')
# f5_5_0.2[, n_loci := 5]
# f5_5_0.2[, effect_size := 0.2]
# 
# f5_10_0.2[, n_loci := 10]
# f5_10_0.2[, effect_size := 0.2]
# f5_5_0.4[, n_loci := 5]
# f5_5_0.4[, effect_size := 0.4]
# 
# dat <- rbindlist(list(f5_5_0.2, f5_5_0.4, f5_10_0.2))
# 
# dat.ag <- dat[,list("mean_x"=mean(mean_x), "mean_y"=mean(mean_y), "sd_y"=sd(mean_y)), by=list(grp, n_loci, effect_size)]
# dat.ag[, bottom_y := mean_y - sd_y]
# dat.ag[, top_y := mean_y + sd_y]
# dat.ag[, type := paste(n_loci, effect_size, sep="_")]
# ggplot(dat.ag, aes(x=mean_x, y=mean_y, fill=type, group=type, ymin=bottom_y, ymax=top_y)) + geom_line() + geom_ribbon(alpha=0.2) + scale_x_reverse()
# 
# 
# f5 <- fread('f5_means.txt')
# 
# ggplot(
#   o[,list(
#     "mean_x"=mean(mean_x), 
#     "mean_y"=mean(mean_y), 
#     "top"=quantile(mean_y, 0.75), 
#     "bottom"=quantile(mean_y, 0.25)
#   ), 
#   by=grp]
# 
# ggplot(o, aes(x=mean_x, y=mean_y, color=factor(replicate))) + geom_line() + scale_x_reverse()
# 
# 
# 
# dat <- data.table(x=res$specificities, y=res$sensitivities)
# dat[, grp := ceiling(500*1:(.N)/nrow(dat))]
# dat.ag <- dat[, list("mean_x"=mean(x), "mean_y"=mean(y)), by=grp]
# 
# ggplot(dat.ag, aes(x=V1, y=V2)) + scale_x_reverse() + geom_line()
# 
# 
# getAuc <- function(replicate, distance=0) {
#   freqs <- fread(paste(replicate, ".freqs.txt", sep=""))
#   causative <- fread(paste(replicate, ".causative.txt", sep=""))
#   setkey(freqs, CHROM, POS)
#   setkey(causative, CHROM, POS)
# 
#   # To Do: add distance metric
# 
#   # causative / case SNP = 1, all others are controls = 0
#   freqs[, causative_snp := 0]
#   freqs[causative, causative_snp := 1]
#   freqs[, minuslogp := -1*log10(P)]
#   res <- roc(freqs, response=causative_snp, predictor=minuslogp)
#   return(res)
# }
# 
# loadData <- function() {
# 
# 
# }
# 
# count_vals <- 1:10
# population_vals <- c("hybrid_swarm_32_F5")
# innerReplicate_vals <- 1:50
# n_loci <- 5
# effect_size <- 0.2
# tag_vals <- c("5_loci_0.2_effect_size")
# 
# o2 <- foreach(count=count_vals, .combine="rbind") %do% {
#   foreach(innerReplicate=innerReplicate_vals, .combine="rbind") %do% {
# 
#   }
# }
