#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)
library(foreach)


loadData <- function() {
    # Load in data
    accuracy.DT <- fread('reconstruction_accuracy.tab', header=TRUE)
    recombination.DT <- fread('reconstruction_recombination_estimates.tab', header=TRUE)
    setnames(recombination.DT, "ind", "individual_iteration")
    # Exclude sex chromosomes, which do not have diploid genotype for all individuals
    recombination.DT <- recombination.DT[! chromosome %in% c("X","chr5")]
    accuracy.DT <- accuracy.DT[! chromosome %in% c("X","chr5")]

    # cast from long to wide, combining both haplotype recombination measurements into a single row
    recombination.wide.DT <- dcast(recombination.DT, individual_iteration+chromosome+population+coverage+nFounders+Ne+mu+population_iteration~haplotype, value.var=c("nRecombEst","nRecombTrue"))

    setkey(accuracy.DT, population, chromosome, individual_iteration, population_iteration, nFounders, coverage, Ne, mu)
    setkey(recombination.wide.DT, population, chromosome, individual_iteration, population_iteration, nFounders, coverage, Ne, mu)

    dat.merge <- merge(accuracy.DT, recombination.wide.DT)
    dat.merge[, recombination_Actual := (nRecombTrue_1 + nRecombTrue_2) ]
    dat.merge[, recombination_Estimate := (nRecombEst_1 + nRecombEst_2) ]
    dat.merge[, recombination_Estimate_delta := recombination_Estimate - recombination_Actual]
    dat.merge[, gid := .GRP, by=list(population, nFounders, coverage, Ne, mu)]
    dat.merge[population=="dgrp", "Hyper-recombinant" := ifelse(recombination_Estimate > 8, TRUE, FALSE)]
    dat.merge[population=="scrm", "Hyper-recombinant" := ifelse(recombination_Estimate > 9, TRUE, FALSE)]
    #dat.merge[recombination_Estimate > 10, "Recombinations" := "\u226511"]
    #dat.merge[recombination_Estimate <= 10, "Recombinations" := "\u226410"]
    dat.merge[, nFounders := factor(nFounders, levels=c(32, 128))]
    dat.merge[, coverage := factor(coverage, levels=c("0.005X","0.05X"))]

    return(dat.merge)
}
dat.merge <- loadData()


plot_subset.beeswarm <- function() {

    # manually jitter
    dat.jittered <- copy(dat.merge[ (population=="dgrp" & coverage=="0.05X") | (Ne==1e6 & mu==5e-9)])
    dat.jittered[percentMatches > 0.99999, percentMatches := sample(c(0.99999,0.999991), size=.N, replace=T)]

    dat.jittered[Ne==1e+04, Ne_label := "N[e]: 10^4"]
    dat.jittered[Ne==1e+05, Ne_label := "N[e]: 10^5"]
    dat.jittered[Ne==1e+06, Ne_label := "N[e]: 10^6"]
    dat.jittered[mu==1e-09, mu_label := "mu: 10^-9"]
    dat.jittered[mu==5e-09, mu_label := "mu: 5 %*% 10^-9"]
    dat.jittered[mu==1e-08, mu_label := "mu: 10^-8"]
    dat.jittered[, Ne_label := factor(Ne_label, levels=c("N[e]: 10^4","N[e]: 10^5","N[e]: 10^6"))]
    dat.jittered[, mu_label := factor(mu_label, levels=c("mu: 10^-9","mu: 5 %*% 10^-9","mu: 10^-8"))]
    dat.jittered[population=="dgrp", population_label := "DGRP-derived"]
    dat.jittered[population=="scrm", population_label := "Coalescent-derived"]



    ggplot(data=dat.jittered, mapping=aes(x=nFounders, y=percentMatches, shape=`Hyper-recombinant`)) +
    geom_beeswarm(data=dat.jittered, size=0.8, alpha=0.8, dodge.width=0.5) +
    scale_alpha_manual(values=c(0.5,1)) +
    scale_shape_manual(values=c(21,4)) +
    facet_grid(.~population_label) +
    theme_few(10) +
    scale_y_continuous(trans="logit", breaks=c(0.5, 0.90, 0.99, 0.999, 0.9999, 0.99999), labels=c("50%", "90%", "99%", "99.9%", "99.99%", "\u226599.999%")) +
    labs(x="Founding Lines", y="Percent Sites with Correct Genotype Estimate") +
        guides(shape = guide_legend(override.aes = list(size = 1)))
}

g <- plot_subset.beeswarm()

ggsave(g, file="Figure_04_reconstruction_accuracy.png", width=16, height=12, units="cm")















####

plot_dgrp.Recomb.beeswarm <- function() {
    dat.jittered <- copy(dat.merge[population=="dgrp"])

    dat.jittered[percentMatches > 0.99999, percentMatches := sample(c(0.99999,0.999991), size=.N, replace=T)]

    ggplot(data=dat.jittered, mapping=aes(x=nFounders, shape=`Hyper-recombinant`, y=percentMatches)) +
    geom_beeswarm(data=dat.jittered, size=0.8, alpha=0.8, dodge.width=0.75) +
    scale_shape_manual(values=c(21,4)) +
    facet_grid(~coverage, labeller="label_both") +
    theme_few(10) +
    scale_y_continuous(trans="logit", breaks=c(0.90, 0.99, 0.999, 0.9999, 0.99999), labels=c("90%", "99%", "99.9%", "99.99%", "\u226599.999%")) +
    labs(x="Founding Lines", y="Percent Sites with Correct Genotype Estimate", color="Estimated Recombination Count", shape="Hyper-recombinant") +
    guides(shape = guide_legend(override.aes = list(size = 3)))

}
#   plot_dgrp.Recomb.beeswarm()

plot_scrm.Recomb.beeswarm <- function() {
    dat.jittered <- copy(dat.merge[population=="scrm"])

    dat.jittered[percentMatches > 0.99999, percentMatches := sample(c(0.99999,0.999991), size=.N, replace=T)]


    dat.jittered[Ne==1e+04, Ne_label := "N[e]: 10^4"]
    dat.jittered[Ne==1e+05, Ne_label := "N[e]: 10^5"]
    dat.jittered[Ne==1e+06, Ne_label := "N[e]: 10^6"]
    dat.jittered[mu==1e-09, mu_label := "mu: 10^-9"]
    dat.jittered[mu==5e-09, mu_label := "mu: 5 %*% 10^-9"]
    dat.jittered[mu==1e-08, mu_label := "mu: 10^-8"]
    dat.jittered[, Ne_label := factor(Ne_label, levels=c("N[e]: 10^4","N[e]: 10^5","N[e]: 10^6"))]
    dat.jittered[, mu_label := factor(mu_label, levels=c("mu: 10^-9","mu: 5 %*% 10^-9","mu: 10^-8"))]


    ggplot(data=dat.jittered, mapping=aes(x=nFounders, shape=`Hyper-recombinant`, y=percentMatches)) +
    geom_beeswarm(data=dat.jittered, size=0.8, alpha=0.8, dodge.width=0.5) +
    scale_shape_manual(values=c(21,4)) +
    scale_alpha_manual(values=c(0.5,1)) +
    facet_grid(mu_label~Ne_label, labeller="label_parsed") +
    theme_few(10) +
    scale_y_continuous(trans="logit", breaks=c(0.5, 0.90, 0.99, 0.999, 0.9999, 0.99999), labels=c("50%", "90%", "99%", "99.9%", "99.99%", "\u226599.999%")) +
    labs(x="Founding Lines", y="Percent Sites with Correct Genotype Estimate", color="Estimated Recombination Count", shape="Estimated Recombination Count") +
    guides(shape = guide_legend(override.aes = list(size = 3)))
}

#   plot_scrm.Recomb.beeswarm()




dgrp_accuracy.g <- plot_dgrp.Recomb.beeswarm()
ggsave(dgrp_accuracy.g, file="Figure_S09_dgrp_imputation_accuracy.png", height=12, width=16, units="cm")

scrm_accuracy.g <- plot_scrm.Recomb.beeswarm()
ggsave(scrm_accuracy.g, file="Figure_S08_scrm_imputation_accuracy.png", height=16, width=20, units="cm")




####









get_dgrp.all <- function() {
    # Accuracy when recombination count estimate is off?
    dat.dgrp.all <- foreach("q"=seq(0.01,1,0.0025), .combine="rbind") %do% {
        dat.merge[population=="dgrp", list("quantile"=q, "accuracy"=quantile(percentMatches, q)), by=list(nFounders, coverage)]
    }
    return(dat.dgrp.all)
}


get_dgrp.CI95 <- function() {
    dgrp.CI95 <- dcast(dat.dgrp.all[quantile %in% c(0.025, 0.25, 0.5, 0.75, 0.975)], nFounders+coverage~quantile, value.var="accuracy")
    setnames(dgrp.CI95, c("0.025", "0.25", "0.5", "0.75", "0.975"), c("q0.025", "q0.25", "q0.5", "q0.75", "q0.975"))
    return(dgrp.CI95)
}

plot_dgrp.all <- function() {
    ggplot(data=dat.dgrp.all, mapping=aes(x=quantile, color=nFounders, shape=nFounders)) +
    scale_color_manual(values=c("blue","red")) +
    scale_shape_manual(values=c(1,3)) +
    geom_boxplot(data=dgrp.CI95, mapping=aes(width=0.2, x=0.014, ymin=q0.025, lower=q0.25, middle=q0.5, upper=q0.75, ymax=q0.975), stat="identity", show.legend=F) +
    geom_point(alpha=0.4, aes(y=accuracy)) +
    scale_x_log10(breaks=c(0.01, 0.1, 0.5, 1), labels=c("0.01", "0.1", "0.5", "1")) +
    scale_y_continuous(trans="logit", breaks=c(0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1), labels=c("86%","88%","90%","92%","94%","96%","98%","100%")) +
    facet_grid(.~coverage, labeller="label_both") +
    theme_few(18) +
    labs(x="Quantile", y="Imputation Accuracy", title="Accuracy for simulated DGRP imputations", color="Founding Lines", shape="Founding Lines") +
    geom_vline(xintercept=0.01, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.1, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.2, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.3, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.4, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.5, linetype="dashed", alpha=0.3) +
    geom_vline(xintercept=0.6, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.7, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.8, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.9, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=1, linetype="dashed", alpha=0.1) +
    theme(legend.direction = "horizontal") +
    theme(legend.position = "bottom")
}

plot_dgrp.beeswarm <- function() {
    dat.jittered <- copy(dat.merge[population=="dgrp"])

    dat.jittered[percentMatches > 0.99999, percentMatches := sample(c(0.99999,0.999991), size=.N, replace=T)]

    ggplot(data=dat.jittered, mapping=aes(x=nFounders, color=nFounders, y=percentMatches)) +
    geom_beeswarm(size=0.5) +
    facet_grid(~coverage) +
    theme_few(18) +
    scale_y_continuous(trans="logit", breaks=c(0.90, 0.99, 0.999, 0.9999, 0.99999), labels=c("90%", "99%", "99.9%", "99.99%", "99.999% +"))
}




plot_dgrp_quantiles <- function () {
    dgrp.byRecomb.CI95 <- dcast(dat.dgrp.byRecombEst[quantile %in% c(0.025, 0.25, 0.5, 0.75, 0.975)], nFounders+coverage+Recombinations~quantile, value.var="accuracy")
    setnames(dgrp.byRecomb.CI95, c("0.025", "0.25", "0.5", "0.75", "0.975"), c("q0.025", "q0.25", "q0.5", "q0.75", "q0.975"))

    g.dgrp.byRecomb <- ggplot(data=dat.dgrp.byRecombEst, mapping=aes(x=quantile, color=nFounders, shape=nFounders)) +
    scale_color_manual(values=c("blue","red")) +
    scale_shape_manual(values=c(1,3)) +
    geom_boxplot(data=dgrp.byRecomb.CI95, mapping=aes(width=0.2, x=0.014, ymin=q0.025, lower=q0.25, middle=q0.5, upper=q0.75, ymax=q0.975), stat="identity", show.legend=F) +
    geom_point(alpha=0.4, aes(y=accuracy)) +
    scale_x_log10(breaks=c(0.01, 0.1, 0.5, 1), labels=c("0.01", "0.1", "0.5", "1")) +
    scale_y_log10(breaks=c(0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1), labels=c("86%","88%","90%","92%","94%","96%","98%","100%")) +
    facet_grid(Recombinations~coverage, labeller="label_both") +
    theme_few(18) +
    labs(x="Quantile", y="Imputation Accuracy", title="Accuracy for simulated DGRP imputations\nby Estimated Recombination Counts", color="Founding Lines", shape="Founding Lines") +
    geom_vline(xintercept=0.01, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.1, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.2, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.3, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.4, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.5, linetype="dashed", alpha=0.3) +
    geom_vline(xintercept=0.6, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.7, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.8, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.9, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=1, linetype="dashed", alpha=0.1) +
    theme(legend.direction = "horizontal") +
    theme(legend.position = "bottom")
}

plot_scrm.quantiles <- function() {
    dat.scrm.q <- foreach("q"=seq(0.01,1,0.001), .combine="rbind") %do% {
        dat.merge[population=="scrm", list("quantile"=q, "accuracy"=quantile(percentMatches, q)), by=list(nFounders, Ne, mu)]
    }

    dat.scrm.q[Ne==1e+04, Ne_label := "N[e]: 10^4"]
    dat.scrm.q[Ne==1e+05, Ne_label := "N[e]: 10^5"]
    dat.scrm.q[Ne==1e+06, Ne_label := "N[e]: 10^6"]
    dat.scrm.q[mu==1e-09, mu_label := "mu: 10^-9"]
    dat.scrm.q[mu==5e-09, mu_label := "mu: 5 %*% 10^-9"]
    dat.scrm.q[mu==1e-08, mu_label := "mu: 10^-8"]

    dat.scrm.q[, Ne_label := factor(Ne_label, levels=c("N[e]: 10^4","N[e]: 10^5","N[e]: 10^6"))]
    dat.scrm.q[, mu_label := factor(mu_label, levels=c("mu: 10^-9","mu: 5 %*% 10^-9","mu: 10^-8"))]

    ggplot(data=dat.scrm.q, mapping=aes(x=quantile, y=accuracy, shape=factor(nFounders), color=factor(nFounders))) +
    scale_color_manual(values=c("blue","red")) +
    scale_shape_manual(values=c(1,3)) +
    geom_point(alpha=0.4) +
    scale_x_log10(breaks=c(0.01, 0.1, 0.5, 1), labels=c("0.01", "0.1", "0.5", "1")) +
    facet_grid(Ne_label~mu_label, labeller="label_parsed", scales="free") +
    theme_few(18) +
    labs(x="Quantile", y="Percent sites imputed correctly", title="Accuracy for seemingly highly recombinant imputations", shape="Founders", color="Founders") +
    geom_vline(xintercept=0.01, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.1, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.2, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.3, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.4, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.5, linetype="dashed", alpha=0.3) +
    geom_vline(xintercept=0.6, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.7, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.8, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=0.9, linetype="dashed", alpha=0.1) +
    geom_vline(xintercept=1, linetype="dashed", alpha=0.1) +
    theme(legend.direction = "horizontal") +
    theme(legend.position = c(0.85, -0.09))
}


# Accuracy by recombination count group
acc_by_recomb <- dat.merge[, list("min"=min(percentMatches), "max"=max(percentMatches), "median"=median(percentMatches), "mean"=mean(percentMatches), "percentAbove99"=sum(percentMatches >= 0.99)/.N, "percentAbove99.9"=sum(percentMatches >= 0.999)/.N, "sd"=sd(percentMatches), "N"=.N, "two_SE"=2*sd(percentMatches)/sqrt(.N)), by=list(coverage, Ne, mu, Recombinations, population, nFounders)]

# Accuracy overall, irrespective of recombination counts
acc_disregarding_recomb <- dat.merge[, list("min"=min(percentMatches), "max"=max(percentMatches), "median"=median(percentMatches), "mean"=mean(percentMatches), "percentAbove99"=sum(percentMatches >= 0.99)/.N, "percentAbove99.9"=sum(percentMatches >= 0.999)/.N, "sd"=sd(percentMatches), "N"=.N, "two_SE"=2*sd(percentMatches)/sqrt(.N)), by=list(coverage, Ne, mu, population, nFounders)]
