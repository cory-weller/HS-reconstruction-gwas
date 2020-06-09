library(data.table)
library(epiR)
library(ggplot2)


dat <- fread("reconstruction_recombination_estimates.tab")

dat2 <- dat[! chromosome %in% c("X","chr5") ,list(totRecombTrue=sum(nRecombTrue), totRecombEst=sum(nRecombEst)), by=list(ind, chromosome,coverage,nFounders,Ne,mu,population_iteration)]

v1 <- dat2[is.na(Ne) & coverage=="0.05X" & nFounders==32]$totRecombTrue
v2 <- dat2[is.na(Ne) & coverage=="0.05X" & nFounders==32]$totRecombEst


dat3 <- dat2[, list("rho"=unlist(epi.ccc(totRecombEst, totRecombTrue)[[1]][1]), "meanDelta"=mean(totRecombEst-totRecombTrue), "sdDelta"=sd(totRecombEst-totRecombTrue)), by=list(coverage, nFounders, Ne, mu)]

fwrite(dat3, file="recombination_correlation_coefficients.tab", row.names=F, col.names=T, quote=F, sep="\t")
