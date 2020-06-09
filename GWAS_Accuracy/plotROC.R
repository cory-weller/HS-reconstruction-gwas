#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)


get_pAUC <- function(dt, pAUC_threshold) {
  dt1 <- dt[x1 < pAUC_threshold]
  dt1[which.max(x1), x2 := pAUC_threshold]
  sum(dt1[, list(y*(x2-x1))][, V1])
}

if(file.exists("AUCs.tab")) {
  o <- fread("AUCs.tab")
} else {
  dat <- rbindlist(lapply(list.files(path="data/", pattern="*[0-9].RDS", full.names=TRUE), function(x) readRDS(x)))
  dat <- rbindlist(lapply(list.files(pattern="*[0-9].RDS"), function(x) readRDS(x)))

  o <- foreach(filename=list.files(path="data/", pattern="*[0-9].RDS", full.names=TRUE), i=1:600) %dopar% {
    print(i)
    dat <- readRDS(filename)
    dat.o <- foreach(j=1:50, .combine="rbind") %do% {
      AUC <- get_pAUC(dat[replicate==j], 1)
      return(cbind(dat[replicate==j], AUC))
    }
    return(dat.o)
  }

  o <- rbindlist(o)
  fwrite(o, file="AUCs.tab", row.names=F, col.names=T, quote=F, sep="\t")
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
o[n_founders==128 & n_generations=="F0", label := "ILs (DGRP 128)"]
o[n_founders==32 & n_generations=="F0", label := "ILs (DGRP 32)"]



ggplot(o[y==0], aes(x=effect_size, y=AUC)) + geom_boxplot() + facet_grid(n_loci ~ population)

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

o.ag <- foreach(x=seq(0.005,1,0.005), .combine="rbind") %do% {
  o[x1 <= x & x2 >= x, list("x"=x, "y"=mean(y)), by=list(effect_size, n_loci, population, label)]
}

o.ag[, n_loci := as.character(n_loci)]
o.ag[n_loci==1, n_loci := "Single-Locus"]
o.ag[n_loci==5, n_loci := "5 Loci"]
o.ag[n_loci==10, n_loci := "10 Loci"]
o.ag[, n_loci := factor(n_loci, levels=c("Single-Locus", "5 Loci", "10 Loci"))]

o.ag[, effect_size := as.character(effect_size)]
o.ag[effect_size==0.2, effect_size := "Smaller Effect"]
o.ag[effect_size==0.4, effect_size := "Larger Effect"]
o.ag[, effect_size := factor(effect_size, levels=c("Smaller Effect", "Larger Effect"))]


g1 <- ggplot(o.ag[label %in% c("ILs (DGRP 128)", "HS_32_F5", "HS_128_F5","RILs (DSPR)", "Outbred_128_F50")], aes(x=x, y=y, color=label)) +
geom_line() +
geom_abline(slope=1, intercept=0) +
xlim(0,1) + ylim(0,1) +
facet_grid(n_loci ~ effect_size) +
labs(x="1 - Specificity", y="Sensitivity") +
theme_few(10)

g2 <- ggplot(o[y==0 & label %in% c("ILs (DGRP 128)", "HS_32_F5", "HS_128_F5","RILs (DSPR)", "Outbred_128_F50")], aes(x=label, y=AUC, color=label)) +
facet_grid(n_loci ~ effect_size) +
theme_few(10) +
stat_summary(fun.data = quantiles_95, geom="boxplot", position=position_dodge(1))

ggsave(g1, file="ROC.svg", width=20, height=20, units="cm")
ggsave(g2, file="ROC.boxplots.svg", width=12, height=14, units="cm")



#### Supplement Plot

g3 <- ggplot(o.ag[label %in% c("ILs (DGRP 128)", "ILs (DGRP 32)", "HS_128_F1", "HS_128_F2", "HS_128_F5", "HS_32_F1", "HS_32_F2", "HS_32_F5")], aes(x=x, y=y, color=label)) +
geom_line() +
geom_abline(slope=1, intercept=0) +
xlim(0,1) + ylim(0,1) +
facet_grid(n_loci ~ effect_size) +
labs(x="1 - Specificity", y="Sensitivity") +
theme_few(10)

g4 <- ggplot(o[y==0 & label %in% c("ILs (DGRP 128)", "ILs (DGRP 32)", "HS_128_F1", "HS_128_F2", "HS_128_F5", "HS_32_F1", "HS_32_F2", "HS_32_F5")], aes(x=label, y=AUC, color=label)) +
facet_grid(n_loci ~ effect_size) +
theme_few(10) +
stat_summary(fun.data = quantiles_95, geom="boxplot", position=position_dodge(1))



ggsave(g3, file="ROC.supplement.svg", width=20, height=20, units="cm")
ggsave(g4, file="ROC.boxplots.supplement.svg", width=12, height=14, units="cm")




o.ag[, c("n_founders","n_generations") := tstrsplit(population, "_")[3:4]]
o.ag[population %like% "outbred_", c("n_founders","n_generations") := tstrsplit(population, "_")[2:3]]
o.ag[population %like% "RILs_", "n_founders" := "8" ]
o.ag[population %like% "RILs_", "n_generations" := "F50" ]




o[, list("AUC"=get_pAUC(


# Plot 3 types comparison
ggplot(o[n_generations != "F1" & n_generations != "F2" & ! (n_founders == 32 & n_generations=="F0")], aes(x=x, y=y, color=label)) + geom_line() +
geom_abline(slope=1, intercept=0) + xlim(0,1) + ylim(0,1) + facet_grid(n_loci ~ effect_size, labeller="label_both") +
labs(x="1 - Specificity", y="Sensitivity") +
theme_few(14)

ggplot(o[n_generations %in% c("F0","F1","F2","F5")], aes(x=x, y=y, color=label)) + geom_line() +
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
