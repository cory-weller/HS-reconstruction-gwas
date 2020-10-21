#!/usr/bin/env Rscript

library(data.table)
library(foreach)
library(cowplot)
library(ggthemes)

setwd("/scratch/caw5cv/genome-reconstruction-revision/05_GWAS/data")

populations <- c("RILs_8_F50", "outbred_128_F50", "hybrid_swarm_128_F5", "hybrid_swarm_32_F5", "hybrid_swarm_128_F0")



getManhattan <- function(population, iteration, PVE, replicate, nLoci) {
    zipFilename <- paste(population, "_", iteration, "_", nLoci, "_", PVE, ".zip", sep="")
    causativeLociFilename <-  paste(replicate, ".causative.txt", sep="")
    freqsFilename <- paste(replicate, ".freqs.txt", sep="")
    freqs <- fread(cmd = paste("unzip -p ", zipFilename, " ", freqsFilename, sep=""))
    #freqs <- freqs[case.Ref + control.Ref > 100 & case.Alt + control.Alt > 100]
    #return(freqs)
    causativeLoci <- fread(cmd = paste("unzip -p ", zipFilename, " ", causativeLociFilename, sep=""))
    return(ggplot(freqs[P < 1], aes(x=POS, y = -1 * log10(P))) + geom_point(shape=21, alpha=0.5) + geom_vline(xintercept = causativeLoci$POS, color="red"))
}


if(! file.exists("gwas_causal_loci.tab")) {
    o <- foreach(population = populations, .combine = "rbind") %do% {
        foreach(nLoci = c("1", "10"), .combine = "rbind") %do% {
            foreach(PVE = c("0.2", "0.4"), .combine = "rbind") %do% {
                foreach(iteration = 1:10, .combine = "rbind") %do% {
                    foreach(replicate = 1:50, .combine = "rbind") %do% {
                        zipFilename <- paste(population, "_", iteration, "_", nLoci, "_", PVE, ".zip", sep="")
                        causativeLociFilename <-  paste(replicate, ".causative.txt", sep="")
                        freqsFilename <- paste(replicate, ".freqs.txt", sep="")
                        #freqs <- fread(cmd = paste("unzip -p ", zipFilename, " ", freqsFilename, sep=""))
                        causativeLoci <- fread(cmd = paste("unzip -p ", zipFilename, " ", causativeLociFilename, sep=""))
                        causativeLoci[, "iteration" := iteration]
                        causativeLoci[, "replicate" := replicate]
                        causativeLoci[, "population" := population]
                        causativeLoci[, "nLoci" := nLoci]
                        causativeLoci[, "PVE" := PVE]
                        return(causativeLoci)
                    }
                }
            }
        }
    }

fwrite(o, file = "gwas_causal_loci.tab", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
} else {
    o <- fread("gwas_causal_loci.tab")
}

singles <- o[nLoci == 1]
singles[population == "RILs_8_F50" & refFreq == 7/8, singleton := TRUE]
singles[population == "outbred_128_F50" & refFreq == 127/128, singleton := TRUE]
singles[population == "hybrid_swarm_128_F5" & refFreq == 127/128, singleton := TRUE]
singles[population == "hybrid_swarm_32_F5" & refFreq == 31/32, singleton := TRUE]
singles[population == "hybrid_swarm_128_F0" & refFreq == 127/128, singleton := TRUE]
singles[is.na(singleton), singleton := FALSE]
singles <- singles[singleton == TRUE & replicate == 1]


mytheme <- theme_few()

g1 <- getManhattan(population = "RILs_8_F50", iteration = 1, PVE=0.4, replicate = 6, nLoci = 1) + mytheme + labs(x="RILs")
g2 <- getManhattan(population = "hybrid_swarm_128_F0", iteration = 2, PVE=0.4, replicate = 7, nLoci = 1) + mytheme + labs(x="Inbred Lines")
g3 <- getManhattan(population = "hybrid_swarm_32_F5", iteration = 2, PVE=0.4, replicate = 7, nLoci = 1) + mytheme + labs(x="Hybrid Swarm (32)")
g4 <- getManhattan(population = "hybrid_swarm_128_F5", iteration = 1, PVE=0.4, replicate = 49, nLoci = 1) + mytheme + labs(x="Hybrid Swarm (128)")
g5 <- getManhattan(population = "outbred_128_F50", iteration = 2, PVE=0.4, replicate = 29, nLoci = 1) + mytheme + labs(x="Outbred F50")

realizations <- plot_grid(g1, g2, g3, g4, g5, ncol=1)

ggsave(realizations, file="realizations.svg", width=10, height=20, units="cm")


getManhattan(population = "RILs_8_F50", iteration = 4, PVE=0.2, replicate = 1, nLoci = 1)
getManhattan(population = "hybrid_swarm_128_F0", iteration = 4, PVE=0.2, replicate = 1, nLoci = 1)

population = "", iteration = 7, replicate = 1, nLoci = 1
population = "hybrid_swarm_128_F5", iteration = 4, replicate = 1, nLoci = 1
population = "", iteration = 10, replicate = 1, nLoci = 1
population = "", iteration = 2, replicate = 1, nLoci = 1


singles[, id := rleid(population)]
singles[, N := 1:.N]
singles[which.min(N), by=id]


o[refFreq >= 0.93][refFreq <= 0.97][nLoci == 1]

For 10 loci, no such determination (averaged over 10 loci, harder to set specific param set)




5% ref freq 95%
comparable effect size



freq of 0.1
effect of 0.2 each