#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)
library(foreach)
library(ggthemes)

founders_128 <- scan("dgrp_128.founders", character())
o1 <- foreach(i = 1:10, .combine = "rbind") %do% {
    dt <- fread(paste("dgrp_128_", i, ".freqs", sep=""))
    dt[, "i" := i]
    setnames(dt, c("CHROM", "start", "stop", founders_128, "individual"))
}
o1[, population := "HS F5 128"]
o1[, Mb := stop/1e6]
haplotypes_128 <- fread("dgrp_128.haps")
o1.long <- melt(o1, measure.vars=founders_128, variable.name="Founder", value.name="HARP freq output")


founders_32 <- scan("dgrp_32.founders", character())
o2 <- foreach(i = 1:10, .combine = "rbind") %do% {
    dt <- fread(paste("dgrp_32_", i, ".freqs", sep=""))
    dt[, "i" := i]
    setnames(dt, c("CHROM", "start", "stop", founders_32, "individual"))
}
o2[, population := "HS F5 32"]
o2[, Mb := stop/1e6]
haplotypes_32 <- fread("dgrp_32.haps")
o2.long <- melt(o2, measure.vars=founders_32, variable.name="Founder", value.name="HARP freq output")

# bad
g.bad <- ggplot(o1.long[individual == 5], aes(x = Mb, y = Founder, fill=`HARP freq output`)) +
geom_tile() +
facet_grid(individual~CHROM) +
theme_few(8) +
theme( strip.background = element_blank(),
  strip.text.y = element_blank()
)
# good
g.good <- ggplot(o1.long[individual == 8], aes(x = Mb, y = Founder, fill=`HARP freq output`)) +
geom_tile() +
facet_grid(individual~CHROM) +
theme_few(8) +
theme( strip.background = element_blank(),
  strip.text.y = element_blank()
)

ggsave(g.bad, file="HARP_ambiguous.png", width=20, height=25, units="cm") 
ggsave(g.good, file="HARP_unambiguous.png", width=20, height=25, units="cm") 