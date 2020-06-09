#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(cores=12)

dat <- rbindlist(lapply(list.files(pattern="*[0-9].RDS"), function(x) readRDS(x)))

dat[, grp := .GRP, by=list(replicate, count, effect_size, n_loci, population)]

get_pAUC <- function(dt, pAUC_threshold) {
  dt1 <- dt[x1 < pAUC_threshold]
  dt1[which.max(x1), x2 := pAUC_threshold]
  sum(dt1[, list(y*(x2-x1))][, V1])
}

pAUCs_mini <- foreach(pAUC=c(0.5, 1), .combine="rbind") %dopar% {
  print(pAUC)
  dat[, .("pAUC_threshold"=pAUC, "pAUC"=get_pAUC(.SD, pAUC)), by=.(replicate, count, effect_size, n_loci, population)]
}

fwrite(pAUCs_mini, file="pAUCs.mini", quote=FALSE, row.names=F, col.names=TRUE, sep="\t")


pAUCs <- foreach(pAUC=c(0.000001, 0.00001, 0.0001, 0.001, seq(0.01, 1, 0.01)), .combine="rbind") %dopar% {
  print(pAUC)
  dat[, .("pAUC_threshold"=pAUC, "pAUC"=get_pAUC(.SD, pAUC)), by=.(replicate, count, effect_size, n_loci, population)]
}

saveRDS(pAUCs, file="pAUCs.RDS")


dat.ag <- dat[, .("pAUC"=mean(pAUC)), by=list(effect_size, n_loci, population, pAUC_threshold)]

ggplot(dat.ag, aes(x=pAUC_threshold, y=pAUC/pAUC_threshold, color=n_loci, linetype=effect_size)) + geom_line() + facet_grid(.~population)
