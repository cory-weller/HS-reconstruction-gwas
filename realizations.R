#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(cowplot)
library(ggthemes)

setwd("/scratch/caw5cv/genome-reconstruction-revision/05_GWAS/data")

populations <- c("RILs_8_F50", "outbred_128_F50", "hybrid_swarm_128_F5", "hybrid_swarm_32_F5",
                  "hybrid_swarm_128_F0")


getManhattan <- function(population, iteration, PVE, replicate, nLoci) {
    zipFilename <- paste(population, "_", iteration, "_", nLoci, "_", PVE, ".zip", sep="")
    causativeLociFilename <-  paste(replicate, ".causative.txt", sep="")
    freqsFilename <- paste(replicate, ".freqs.txt", sep="")
    freqs <- fread(cmd = paste("unzip -p ", zipFilename, " ", freqsFilename, sep=""), colClasses=c(POS="numeric",case.Ref="numeric",case.Alt="numeric",control.Ref="numeric",control.Alt="numeric"))
    freqs[, chisq := NULL]
    freqs[, P := NULL]
    freqs[, chisq := (((case.Ref)+(control.Ref)+((case.Alt))+((control.Alt))) * (((control.Ref)*(case.Alt)) - ((case.Ref)*(control.Alt)))**2) / (((case.Ref)+(control.Ref)) * (((case.Ref)+(case.Alt))) * ((control.Ref)+(control.Alt)) * (((case.Alt))+(control.Alt)))]
    freqs[, ref := case.Ref + control.Ref]
    freqs[, alt := case.Alt + control.Alt]
    freqs[, tot := ref + alt ]
    freqs[, MAF := min(ref, alt)/tot, by=POS]
    freqs <- freqs[MAF >= 0.05]
    if(grepl("^RIL", population) == TRUE | grepl("F0$", population) == TRUE) {
        freqs[, chisq := chisq / 2]
        freqs[, halved := TRUE]
    }
    freqs[,P := pchisq(chisq, 1, lower.tail=FALSE)]

    causativeLoci <- fread(cmd = paste("unzip -p ", zipFilename, " ", causativeLociFilename, sep=""))
    return(ggplot(freqs[P < 1], aes(x=POS, y = -1 * log10(P))) + geom_point(shape=21, alpha=0.5) + geom_vline(xintercept = causativeLoci$POS, color="red"))
}


if(! file.exists("gwas_causal_loci.tab")) {
    o <- foreach(population = populations, .combine = "rbind") %do% {
        foreach(nLoci = c("1", "10"), .combine = "rbind") %do% {
            foreach(PVE = c("0.2", "0.4"), .combine = "rbind") %do% {
                foreach(iteration = 1:10, .combine = "rbind") %do% {
                    foreach(replicate = 1:50, .combine = "rbind") %do% {
                        zipFilename <- paste(   population, "_", iteration, "_", nLoci, "_", PVE, ".zip",
                                                sep="")
                        causativeLociFilename <-  paste(replicate, ".causative.txt", sep="")
                        freqsFilename <- paste(replicate, ".freqs.txt", sep="")
                        #freqs <- fread(cmd = paste("unzip -p ", zipFilename, " ", freqsFilename, sep=""))
                        causativeLoci <- fread( cmd = paste("unzip -p ", zipFilename, " ", 
                                                causativeLociFilename, sep="")
                                              )
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

g1 <- getManhattan( population = "RILs_8_F50", 
                    iteration = 1, 
                    PVE=0.4, 
                    replicate = 6, 
                    nLoci = 1) +
        mytheme +
        labs(x="RILs")


g2 <- getManhattan( population = "hybrid_swarm_128_F0", 
                    iteration = 2,
                    PVE=0.4,
                    replicate = 7,
                    nLoci = 1) +
        mytheme +
        labs(x="Inbred Lines")


g3 <- getManhattan( population = "hybrid_swarm_32_F5",
                    iteration = 2,
                    PVE=0.4,
                    replicate = 7,
                    nLoci = 1) +
        mytheme +
        labs(x="Hybrid Swarm (32)")


g4 <- getManhattan( population = "hybrid_swarm_128_F5",
                    iteration = 1,
                    PVE=0.4,
                    replicate = 49,
                    nLoci = 1) +
        mytheme +
        labs(x="Hybrid Swarm (128)")


g5 <- getManhattan( population = "outbred_128_F50",
                    iteration = 2,
                    PVE=0.4,
                    replicate = 29,
                    nLoci = 1) +
        mytheme +
        labs(x="Outbred F50")


realizations <- plot_grid(g1, g2, g3, g4, g5, ncol=1)

ggsave(realizations, file="realizations.png", width=10, height=20, units="cm")
