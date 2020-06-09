#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)
library(cowplot)

founder_selection <- fread('zcat reconstruction_founder_selection_stats.tab.gz', header=TRUE)

founder_selection[, founders_missed := N_true_founders_in_individual - N_true_founders_chosen ]
founder_selection.ag <- founder_selection[! chromosome %in% c("X","chr5"), list(
                    "meanFoundersMissed"=mean(founders_missed),
                    "meanFractionCovered"=mean(fraction_covered)
                ), by=list(population, Ne, mu, nFounders, coverage, maxCount, threshold)]

founder_selection.ag[Ne==1e+04, Ne_label := "N[e]: 10^4"]
founder_selection.ag[Ne==1e+05, Ne_label := "N[e]: 10^5"]
founder_selection.ag[Ne==1e+06, Ne_label := "N[e]: 10^6"]
founder_selection.ag[mu==1e-09, mu_label := "mu: 10^-9"]
founder_selection.ag[mu==5e-09, mu_label := "mu: 5 %*% 10^-9"]
founder_selection.ag[mu==1e-08, mu_label := "mu: 10^-8"]
founder_selection.ag[, Ne_label := factor(Ne_label, levels=c("N[e]: 10^4","N[e]: 10^5","N[e]: 10^6"))]
founder_selection.ag[, mu_label := factor(mu_label, levels=c("mu: 10^-9","mu: 5 %*% 10^-9","mu: 10^-8"))]
founder_selection.ag[, nFounders := factor(nFounders, levels=c(32,128))]
founder_selection.ag[, threshold := factor(threshold, levels=c(0.9, 0.95, 0.99, 0.999))]
founder_selection.ag[population=="scrm", population_label := "Coalescent-derived"]
founder_selection.ag[population=="dgrp", population_label := "DGRP-derived"]
founder_selection.ag[, population_label := factor(population_label, levels=c("DGRP-derived","Coalescent-derived"))]
founder_selection.ag[nFounders==32, founders_label := "32 Founders"]
founder_selection.ag[, founders_label := factor(founders_label, levels=c("32 Founders","128 Founders"))]
founder_selection.ag[nFounders==128, founders_label := "128 Founders"]
setnames(founder_selection.ag, "threshold", "HARP Threshold")


g.main <- ggplot(founder_selection.ag[(population=="scrm" & Ne == 1e+06 & mu == 5e-09) | (population=="dgrp" & coverage == "0.05X")], aes(
    x=maxCount,
    y=meanFractionCovered,
    color=`HARP Threshold`, group=`HARP Threshold`
)) +
geom_line() +
facet_grid(population_label~founders_label) +
ylim(0,1) +
theme_few(10) +
labs(x="Maximum MLA Set Size", y="Fraction of Genome Represented by MLA Set") 

# ggsave(g.main, "founder_selection_accuracy.svg", width=10, height=10, units="cm")

g.supplementA <- ggplot(founder_selection.ag[coverage == "0.05X" & population=="scrm"], aes(
    x=maxCount,
    y=meanFractionCovered,
    color=`HARP Threshold`, 
    linetype=founders_label
)) +
geom_line() +
facet_grid(mu_label~Ne_label, labeller =label_parsed) +
ylim(0,1) +
theme_few(10) +
labs(x="Maximum MLA Set Size", y="Fraction of Genome Represented by MLA Set") +
theme(legend.position = "none")
# guides(linetype=guide_legend(title="Haplotypes"))

g.supplementB <- ggplot(founder_selection.ag[population=="dgrp"], aes(
    x=maxCount,
    y=meanFractionCovered,
    color=`HARP Threshold`, 
    linetype=founders_label
)) +
geom_line() +
facet_grid(coverage~., labeller = label_both) +
theme_few(10) +
labs(x="Maximum MLA Set Size", y="Fraction of Genome Represented by MLA Set") +
guides(linetype=guide_legend(title="Haplotypes"))

g.supplement <- plot_grid(g.supplementA, g.supplementB, labels=c("A","B"), rel_widths=c(1, 0.9), axis="bt")



# ggsave(g.supplement, file="founder_selection_accuracy_supplement.svg", width=15, height=10, units="cm")