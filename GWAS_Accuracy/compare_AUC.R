#!/usr/bin/env Rscript

library(data.table)
library(foreach)


dat <- fread('AUCs.tab')

dat[, c("x1","x2","y") := NULL]
dat <-  unique(dat)

setkey(dat, effect_size, n_loci, population)

effect_sizes <- c(0.2, 0.4)
locus_sizes <- c(1, 5, 10)
populations <- unique(dat$population)

o <- foreach(i = effect_sizes, .combine="rbind") %do% {
    foreach(j = locus_sizes, .combine="rbind") %do% {
        foreach(k = populations, .combine="rbind") %do% {
            foreach(l = populations, .combine="rbind") %do% {
                if(k!=l) {
                    dat.sub <- dat[.(i, j)][population %in% c(k,l)]
                    w = wilcox.test(AUC~population, data=dat.sub)
                    median1 <- median(dat.sub[population==k, AUC])
                    median2 <- median(dat.sub[population==l, AUC])
                    out <- data.table("effect_size" = i,
                                "locus_size" = j,
                                "population1" = k,
                                "population1_median_AUC" = median1,
                                "population2" = l,
                                "population2_median_AUC" = median2,
                                "W_stat" = w$statistic,
                                "W_p" = w$p.value)
                    return(out)
                }
            }
        }
    }
}

fwrite(o[duplicated(o[,W_p])], "table_S2.tab", quote=F, row.names=F, col.names=T, sep="\t")



# 270 comparisons: (10 choose 2 populations) * 3 locus sizes * 2 effect sizes
# Bonferroni threshold of 0.0001851852
