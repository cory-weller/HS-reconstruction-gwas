#!/usr/bin/env R

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


if(! require(ggthemes)) {
  install.packages("ggthemes")
  library(ggthemes)
}

if(! require(ggbeeswarm)) {
  install.packages("ggbeeswarm")
  library(ggbeeswarm)
}



# load in each accuracy table

mainstring <- "hybrid_swarm_32"
generations <- c("F1", "F2", "F5", "F50")



if(! file.exists("RABBIT_accuracy.tab")) {
  o <- foreach(generation=generations, .combine="rbind", .errorhandling="remove") %do% {
    accuracy_filename <- paste(mainstring, "_", generation, "/accuracy.txt", sep="")
    dat <- fread(accuracy_filename)
    dat[, n_generations := generation]
    return(dat)
  }
  fwrite(
    x=o,
    file = "RABBIT_accuracy.tab", 
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
} else {
  o <- fread("RABBIT_accuracy.tab")
}



o[, n_generations := factor(n_generations, levels=generations)]

# jitter accuracy values above 99.999%
o[percent_acc > 0.99999, percent_acc := runif(.N, min=0.99999, max=0.9999999)]


dat.ag <- o[, 
  list(
    qBottom=quantile(percent_acc, 0.025),
    qLower=quantile(percent_acc, 0.25),
    qMiddle=quantile(percent_acc, 0.5),
    qUpper=quantile(percent_acc, 0.75),
    qMax=quantile(percent_acc, 0.975)
  ),
  by=list(n_generations)
]



g <- ggplot(dat.ag, mapping=aes(
                            x=n_generations,
                            min=qBottom,
                            lower=qLower,
                            middle=qMiddle,
                            upper=qUpper,
                            max=qMax)
                            ) +
    geom_boxplot(stat="identity", color="black") +
    scale_y_continuous( trans="logit",
                        breaks=c(0.5, 0.90, 0.99, 0.999, 0.9999, 0.99999),
                        labels=c("50%", "90%", "99%", "99.9%", "99.99%", "\u226599.999%")
                      ) +
    theme_few(12) +
    labs( x="Generations",
          y="Genotype Accuracy (fraction correct sites)"
        )

ggsave(g, file="RABBIT_accuracy.png", width=10, height=10, units="cm")

# Jitter accuracy
o[percent_acc > 0.99999, percent_acc := runif(.N, min=0.99999, max=0.9999999)]

ggplot(data=o, mapping=aes(x=n_generations, y=percent_acc)) +
geom_beeswarm(size=0.8, alpha=0.8, dodge.width=0.5) +
scale_y_continuous(trans="logit",
  breaks=c(0.5, 0.90, 0.99, 0.999, 0.9999, 0.99999, 0.999999 ),
  labels=c("50%", "90%", "99%", "99.9%", "99.99%", "99.999%", "\u226599.9999%")
  )


scale_alpha_manual(values=c(0.5,1)) +
scale_shape_manual(values=c(21,4)) +
facet_grid(.~population_label) +
theme_few(10) +
scale_y_continuous(trans="logit", breaks=c(0.5, 0.90, 0.99, 0.999, 0.9999, 0.99999), labels=c("50%", "90%", "99%", "99.9%", "99.99%", "\u226599.999%")
labs(x="Founding Lines", y="Percent Sites with Correct Genotype Estimate") +
    guides(shape = guide_legend(override.aes = list(size = 1)))

plot_subset.beeswarm <- function(dat) {

    # manually jitter
    dat.jittered <- copy(dat)
    dat.jittered[percentMatches > 0.99999, percentMatches := sample(c(0.99999,0.999991), size=.N, replace=T)]

    ggplot(data=dat.jittered, mapping=aes(x=nFounders, y=percentMatches, shape=`Hyper-recombinant`)) +
    geom_beeswarm(data=dat.jittered, size=0.8, alpha=0.8, dodge.width=0.5) +
    scale_alpha_manual(values=c(0.5,1)) +
    scale_shape_manual(values=c(21,4)) +
    facet_grid(.~population_label) +
    theme_few(10) +
    scale_y_continuous(trans="logit", breaks=c(0.5, 0.90, 0.99, 0.999, 0.9999, 0.99999), labels=c("50%", "90%", "99%", "99.9%", "99.99%", "\u226599.999%")
    labs(x="Founding Lines", y="Percent Sites with Correct Genotype Estimate") +
        guides(shape = guide_legend(override.aes = list(size = 1)))
}

g <- plot_subset.beeswarm()

ggsave(g, file="Figure_04_reconstruction_accuracy.png", width=15, height=10, units="cm")





g <- ggplot(o, mapping=aes(x=N, y=percent_correct)) +
  geom_jitter(alpha=0.2, shape=21) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) +
  facet_grid(n_generations~.) +
  labs(x="Population Size", y="Genotype Accuracy (% sites)", title="STITCH accuracy\n(pseudoHaploid model)")

g <- ggplot(o, mapping=aes(x=N, y=percent_correct, color=n_generations)) +
  geom_point(alpha=0.5, shape=21, position = position_jitterdodge(0.15)) +
  geom_boxplot(outlier.shape=NA, alpha=0.5) +
  labs(x="Population Size", y="Genotype Accuracy (% sites)")

ggsave(g, file="RABBIT_accuracy.png", width=10, height=20, units="cm")



dat.ag <- o[, list(
    qBottom=quantile(percent_correct, 0.025),
    qLower=quantile(percent_correct, 0.25),
    qMiddle=quantile(percent_correct, 0.5),
    qUpper=quantile(percent_correct, 0.75),
    qMax=quantile(percent_correct, 0.975)),
    by=list(N,n_generations)]


g <- ggplot(dat.ag, mapping=aes(
                            x=N,
                            fill=n_generations,
                            min=qBottom,
                            lower=qLower,
                            middle=qMiddle,
                            upper=qUpper,
                            max=qMax)
                            ) +
    geom_boxplot(stat="identity", color="black") +
    labs(x="Population Size", y="Genotype Accuracy (fraction correct sites)") +
    theme_few(12)

ggsave(g, file="STITCH_accuracy.png", width=20, height=10, units="cm")
