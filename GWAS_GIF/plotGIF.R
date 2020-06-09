#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)

if(file.exists("GIF.tab")) {
 dat <- fread("GIF.tab")
} else {
  dat <- rbindlist(
    foreach(filename=list.files(pattern="*.out")) %do% {
      readFile <- fread(filename)
      setnames(readFile, c( "slurm_iteration",
                            "permutation_iteration",
                            "causal_iteration", 
                            "GIF_all", 
                            "GIF_unlinked", 
                            "GIF_linked", 
                            "sum_effect"
                          )
              )
      readFile[, fn := filename]  
    }
  )
  fwrite(dat, file="GIF.tab", quote=F, row.names=F, col.names=T, sep="\t")
}

dat[fn %like% "hybrid_swarm", population := "hybrid_swarm"]

dat[fn %like% "hybrid_swarm", c("n_founders","n_generations","n_loci","effect") := tstrsplit(fn, split="_")[-c(1,2)]]
dat[! fn %like% "hybrid_swarm", c("population","n_founders","n_generations","n_loci","effect") := tstrsplit(fn, split="_")]
dat[, effect := tstrsplit(effect, split=".out")]

dat[, label := paste(population, n_founders, n_generations, sep="_")]

dat[, fn := NULL]

dat.long <- melt(dat, measure.vars=c("GIF_all","GIF_unlinked","GIF_linked"), variable.name="GIF_type", value.name="GIF")



dat.long[, label := factor(label, levels=
                        c("inbredLines_32_F0",
                        "hybrid_swarm_32_F1",
                        "hybrid_swarm_32_F2",
                        "hybrid_swarm_32_F5",
                        "inbredLines_128_F0",
                        "hybrid_swarm_128_F1",
                        "hybrid_swarm_128_F2",
                        "hybrid_swarm_128_F5",
                        "RILs_8_F50",
                        "outbred_128_F50")
        )]
        
dat.long[, facetLabel := factor()]
dat.long[label=="inbredLines_32_F0",  facetLabel := "Inbred~Lines" ]
dat.long[label=="hybrid_swarm_32_F1",  facetLabel := "F[1]~Hybrid~Swarm" ]
dat.long[label=="hybrid_swarm_32_F2",  facetLabel := "F[2]~Hybrid~Swarm" ]
dat.long[label=="hybrid_swarm_32_F5",  facetLabel := "F[5]~Hybrid~Swarm " ]
dat.long[label=="inbredLines_128_F0",  facetLabel := "Inbred~Lines" ]
dat.long[label=="hybrid_swarm_128_F1",  facetLabel := "F[1]~Hybrid~Swarm" ]
dat.long[label=="hybrid_swarm_128_F2",  facetLabel := "F[2]~Hybrid~Swarm" ]
dat.long[label=="hybrid_swarm_128_F5",  facetLabel := "F[5]~Hybrid~Swarm" ]
dat.long[label=="RILs_8_F50",  facetLabel := "Rec.~Inbred~lines " ]
dat.long[label=="outbred_128_F50",  facetLabel := "Outbred" ]

dat.long[, n_loci_label := factor()]
dat.long[n_loci==1, n_loci_label := "Single~Locus~Trait" ]
dat.long[n_loci==10, n_loci_label := "Multi-Locus~Trait" ]

dat.long[, GIF_type_label := factor()]
dat.long[GIF_type=="GIF_all", GIF_type_label := "Genome-Wide"]
dat.long[GIF_type=="GIF_linked", GIF_type_label := "Linked Loci"]
dat.long[GIF_type=="GIF_unlinked", GIF_type_label := "Unlinked Loci"]


        
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

g1 <- ggplot(dat.long[n_founders != 32], aes(x=label, fill=GIF_type_label, y=GIF)) + 
theme_few(10) +
labs(x="Mapping Population", y=parse(text="Genomic~Inflation~Factor~(\u03bb[1000])")) +
facet_grid(n_loci_label~facetLabel, scales="free", switch="x", labeller = label_parsed) +
stat_summary(fun.data = quantiles_95, geom="boxplot", position=position_dodge(1)) +
scale_fill_manual(values=c("white","gray80","gray40")) +
geom_hline(yintercept=1.0, linetype="dashed", alpha=0.6) +
theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
guides(fill=guide_legend(""))

ggsave(g1, file="GIF_128.svg", width=24, height=10, units="cm")

g2 <- ggplot(dat.long[n_founders == 32], aes(x=label, fill=GIF_type_label, y=GIF)) + 
theme_few(10) +
labs(x="Mapping Population", y=parse(text="Genomic~Inflation~Factor~(\u03bb[1000])")) +
facet_grid(n_loci_label~facetLabel, scales="free", switch="x", labeller = label_parsed) +
stat_summary(fun.data = quantiles_95, geom="boxplot", position=position_dodge(1)) +
scale_fill_manual(values=c("white","gray80","gray40")) +
geom_hline(yintercept=1.0, linetype="dashed", alpha=0.6) +
theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
guides(fill=guide_legend(""))

ggsave(g2, file="GIF_32.svg", width=18, height=10, units="cm")

