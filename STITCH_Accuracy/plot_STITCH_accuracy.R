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

quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


# load in each accuracy table

mainstring <- "hybrid_swarm_32"
# generations <- c("F1", "F2", "F5", "F50")
generations <- c("F5")

#population_sizes <- seq(100,1000,100)
population_sizes <- c(100, 500, 1000, 2000, 3000, 4000, 5000)

# o <- foreach(filename = Sys.glob("./*/accuracy.RDS"), .combine="rbind") %do% {
#   readRDS(filename)
# }
# fwrite(o, file="STITCH_accuracies_summary.tab", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

o <- fread("STITCH_accuracies_summary.tab")
o[, N := as.numeric(as.character(N))]

# subset for plotting; exclude N=200-400 and 600-900
o <- o[N %in% c(100,500,1000,2000,3000,4000,5000)]


o[, N_label := factor(N, levels=c(0,100,500,1000,2000,3000,4000,5000))]

o.labels <- o[,list("medianAcc"=median(percent_correct), top95=quantile(percent_correct, 0.975)), by=list(N_label)]
o.labels[,val := paste(round(100*medianAcc,1), "%", sep="")]


g <- ggplot(o, mapping=aes(x=N_label, y=percent_correct)) +
  labs(x="Individuals Sequenced at 0.05X", y="Percent Sites with Correct Genotype Estimate") +
  theme_few(12) +
  stat_summary(fun.data = quantiles_95, geom="boxplot", width=0.4) +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(breaks=c(0.6, 0.7, 0.8, 0.9, 1.0), labels=c("60%", "70%", "80%", "90%", "100%"))
 # geom_text(data=o.labels, aes(x=N_label, y=medianAcc, label=val), size=3, position=position_nudge(x = -0.49, y = 0)) 


  
ggsave(g, file="STITCH_accuracy.svg", width=10,height=10, units="cm")