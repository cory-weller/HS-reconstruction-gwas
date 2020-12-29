#!/usr/bin/env bash

# Download read counts from high coverage sequencing

mkdir -p comparison && cd comparison

wget -O highcov.readcounts.tar.gz "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118833&authkey=AKHmP_3zBkEvdtE"
tar -zxvf highcov.readcounts.tar.gz

tar -zxvf ../4CB-1-73_S73.tar.gz  && mv 4CB-1-73_S73.estimate.genotypes 73.estimate.genotypes
tar -zxvf ../4CB-1-77_S77.tar.gz  && mv 4CB-1-77_S77.estimate.genotypes 77.estimate.genotypes
tar -zxvf ../4CB-1-78_S78.tar.gz  && mv 4CB-1-78_S78.estimate.genotypes 78.estimate.genotypes
tar -zxvf ../4CB-1-79_S79.tar.gz  && mv 4CB-1-79_S79.estimate.genotypes 79.estimate.genotypes
tar -zxvf ../4CB-1-81_S81.tar.gz  && mv 4CB-1-81_S81.estimate.genotypes 81.estimate.genotypes
tar -zxvf ../4CB-1-82_S82.tar.gz  && mv 4CB-1-82_S82.estimate.genotypes 82.estimate.genotypes

#!/usr/bin/env R

library(data.table)
library(foreach)
library(ggplot2)
library(ggrepel)
library(ggthemes)

VCF_filename <- "../dgrp2.haplotypes.vcf.gz"
siteFreqs_filename <- "../siteFreqs.tab"

extractRegions <- function(sites, haplotypes) {
    foreach(i=1:nrow(haplotypes), .combine="rbind") %do% {
        return(sites[contig == haplotypes[i,chromosome] & position >= haplotypes[i,start] & position <= haplotypes[i,stop]])
        }
}

chromosomes <- c("2L","2R","3L","3R","X")

getHaplotypes <- function(ind, chromosome) {
    fread(paste(ind, ".", chromosome, ".estimate.haplotypes", sep=""))
}

importVCF <- function(VCF_filename) {
    VCF_filename_split <- strsplit(VCF_filename, "\\.")[[1]]
    VCF_extension <- VCF_filename_split[length(VCF_filename_split)]
    if(VCF_extension == "vcf") {
        VCF <- fread(VCF_filename)
    } else if(VCF_extension == "gz") {
        VCF <- fread(cmd = paste("zcat", VCF_filename))
    }
    return(VCF)
}

siteFreqs <- fread(siteFreqs_filename, col.names=c("contig", "position", "popRefCount", "popAltCount"))
setkey(siteFreqs, "contig", "position")


o <- foreach( ind = c(73, 77, 78, 79, 81, 82), .combine="rbind" ) %do% {
    haplotypes <- foreach(chromosome = chromosomes, .combine="rbind") %do% {
        dt <- getHaplotypes(ind, chromosome)
        data.table("contig" = chromosome, "n_recomb" = nrow(dt)-1)
    }
    setkey(haplotypes, contig)

    lowcov <- fread(paste(ind, ".estimate.genotypes", sep=""))
    setnames(lowcov, c("contig", "position", "allele1", "allele2"))
    setkey(lowcov, contig, position)
    highcov <- fread(paste(ind, ".highcov.readcounts", sep=""), 
                    select=c(1,2,6,7,8,12,13), 
                    col.names=c("contig","position","refCount","altCount","totalCount","otherBases","improperPairs")
    )
    setkey(highcov, contig,position)
    highcov <- highcov[!duplicated(highcov)]

    dat <- merge(lowcov,highcov)
    dat <- merge(dat, siteFreqs)
    dat[, "ind" := ind ]
    dat <- merge(dat, haplotypes)

    

    # filter out sites with unrecognized bases or improper pairs
    dat <- dat[otherBases==0 & improperPairs==0] # & contig != "X"]
    dat <- dat[which(refCount < quantile(refCount, 0.999) & altCount < quantile(altCount, 0.999))]

    dat[allele1==0 & allele2==0, estGenotype := "HomRef"]
    dat[allele1==1 & allele2==1, estGenotype := "HomAlt"]
    dat[allele1==0 & allele2==1, estGenotype := "Het"]
    dat[allele1==1 & allele2==0, estGenotype := "Het"]
    dat[allele1=="." | allele2==".", estGenotype := NA]

    dat[refCount >= 6 & altCount == 0, trueGenotype := "HomRef"]
    dat[refCount == 0 & altCount >= 6, trueGenotype := "HomAlt"]
    dat[refCount >= 3 & altCount >= 3 , trueGenotype := "Het"]
    dat <- dat[!is.na(trueGenotype) & ! is.na(estGenotype)]
    return(dat)
    #out <- dat[, list("N_Sites"=.N, "accurate_sites"=sum(estGenotype==trueGenotype)), by=list(contig, ind, n_recomb)]
    #out[, acc := accurate_sites / N_Sites ]

    #return(out)
}


o[, refFreq := popRefCount / (popRefCount + popAltCount)]
o[, altFreq := popAltCount / (popRefCount + popAltCount)]
o[, idx := 1:.N]
o[, MAF := min(refFreq, altFreq), by = idx]
o[, MAF_grp := cut(MAF, breaks=c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, Inf))]

dat <- o[, list("N_Sites"=.N, "accurate_sites"=sum(estGenotype==trueGenotype)), by=list(contig, n_recomb, ind)]
dat[, acc := accurate_sites / N_Sites]
dat[, ind := factor(ind)]

fwrite(dat, file="reconstruction_comparison.tab", quote=F, row.names=F, col.names=T, sep="\t")

g <- ggplot(dat, aes(x=n_recomb, y=acc, shape=ind)) + geom_point() + geom_text_repel(aes(label=contig)) + geom_vline(xintercept=10, linetype="dashed") +
theme_few(12) + 
labs(x="Estimated number of recombination events", y="Genotype estimate concordance", color="Individual") +
theme(legend.position = "none")

ggsave(g, file="concordance.svg", height=10, width=15, units="cm")

ggplot(dat, aes(x=n_recomb, y=acc, color=ind))) + facet_grid(trueGenotype~.) + geom_jitter()


## Chromosome Paintings

haplotype_files <- list.files(pattern="*.estimate.haplotypes")

dat <- foreach(fn = haplotype_files, .combine="rbind") %do% {
    ind = strsplit(fn, "\\.")[[1]][1]
    dt <- fread(fn)
    dt[, "ind" := ind]
}

dat <- melt(dat, measure.vars=c("par1","par2"), variable.name="haplotype")
dat[haplotype=="par1", haplotype2 := 1]
dat[haplotype=="par2", haplotype2 := 2]
dat[, haplotype := NULL]
setnames(dat, "haplotype2", "haplotype")
setnames(dat, "value", "lineID")

convert_to_wide <- function(DT) {
  # Collapse redundant breakpoints
  DT[, lineID_rleid := rleid(chromosome, haplotype, ind, lineID)]
  haps.collapsed <- DT[, list(chromosome,haplotype,start=min(start), stop=max(stop), ind, lineID), by=lineID_rleid]
  haps.collapsed[,lineID_rleid := NULL]
  DT <- haps.collapsed[!duplicated(haps.collapsed)]
}

dat <- convert_to_wide(dat)

dat[haplotype=="1", ystart := 0]
dat[haplotype=="1", ystop := 1]
dat[haplotype=="2", ystart := 1]
dat[haplotype=="2", ystop := 2]

mytheme <- theme_few() + theme(legend.position = "none") + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank(),  strip.background = element_blank(),
        strip.text.x = element_blank(), strip.text.y= element_blank())

g1 <- ggplot(dat[chromosome=="2L" & ind==73], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g2 <- ggplot(dat[chromosome=="2R" & ind==73], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g3 <- ggplot(dat[chromosome=="3L" & ind==73], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g4 <- ggplot(dat[chromosome=="3R" & ind==73], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g5 <- ggplot(dat[chromosome=="X"  & ind==73], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g6 <- ggplot(dat[chromosome=="2L" & ind==77], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g7 <- ggplot(dat[chromosome=="2R" & ind==77], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g8 <- ggplot(dat[chromosome=="3L" & ind==77], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g9 <- ggplot(dat[chromosome=="3R" & ind==77], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g10 <- ggplot(dat[chromosome=="X"  & ind==77], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g11 <- ggplot(dat[chromosome=="2L" & ind==78], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g12  <- ggplot(dat[chromosome=="2R" & ind==78], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g13  <- ggplot(dat[chromosome=="3L" & ind==78], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g14  <- ggplot(dat[chromosome=="3R" & ind==78], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g15  <- ggplot(dat[chromosome=="X"  & ind==78], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g16  <- ggplot(dat[chromosome=="2L" & ind==79], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g17 <- ggplot(dat[chromosome=="2R" & ind==79], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g18 <- ggplot(dat[chromosome=="3L" & ind==79], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g19 <- ggplot(dat[chromosome=="3R" & ind==79], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g20 <- ggplot(dat[chromosome=="X"  & ind==79], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g21 <- ggplot(dat[chromosome=="2L" & ind==81], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme

g22 <- ggplot(dat[chromosome=="2R" & ind==81], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g23 <- ggplot(dat[chromosome=="3L" & ind==81], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g24 <- ggplot(dat[chromosome=="3R" & ind==81], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g25 <- ggplot(dat[chromosome=="X"  & ind==81], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g26 <- ggplot(dat[chromosome=="2L" & ind==82], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme
g27 <- ggplot(dat[chromosome=="2R" & ind==82], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme

g28 <- ggplot(dat[chromosome=="3L" & ind==82], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme

g29 <- ggplot(dat[chromosome=="3R" & ind==82], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme

g30 <- ggplot(dat[chromosome=="X"  & ind==82], aes(xmin=start, xmax=stop, ymin=ystart, ymax=ystop, fill=lineID)) + geom_rect() + facet_grid(ind~chromosome) + mytheme

g.out <- plot_grid(plotlist=mget(paste0("g", 1:30)), ncol=5)

paintings <- tempfile("paintings", fileext = ".png")
save_plot(paintings, g.out)



# o <- foreach( ind = c(73, 77, 78, 79, 81, 82), .combine="rbind" ) %do% {
#     haplotypes <- foreach(chromosome = chromosomes, .combine="rbind") %do% {
#         dt <- getHaplotypes(ind, chromosome)
#         dt[, "n_recomb" := nrow(dt)-1]
#     }

#     n_recomb <- nrow(haplotypes)-1

#     lowcov <- fread(paste(ind, ".estimate.genotypes", sep=""))
#     setnames(lowcov, c("contig", "position", "allele1", "allele2"))
#     setkey(lowcov, contig, position)
#     highcov <- fread(paste(ind, ".highcov.readcounts", sep=""), 
#                     select=c(1,2,6,7,8,12,13), 
#                     col.names=c("contig","position","refCount","altCount","totalCount","otherBases","improperPairs")
#     )
#     setkey(highcov, contig,position)
#     highcov <- highcov[!duplicated(highcov)]

#     dat <- merge(lowcov,highcov)

#     goodSites <- extractRegions(dat, haplotypes[stop - start > 1e6])
#     goodSites[, mask := FALSE]

#     badSites <- extractRegions(dat, haplotypes[stop - start < 1e6])
#     badSites[, mask := TRUE]

#     dat <- rbindlist(list(goodSites, badSites))
#     dat[, "ind" := ind ]

#     # filter out sites with unrecognized bases or improper pairs
#     dat <- dat[otherBases==0 & improperPairs==0] # & contig != "X"]
#     dat <- dat[which(refCount < quantile(refCount, 0.999) & altCount < quantile(altCount, 0.999))]

#     dat[allele1==0 & allele2==0, estGenotype := "HomRef"]
#     dat[allele1==1 & allele2==1, estGenotype := "HomAlt"]
#     dat[allele1==0 & allele2==1, estGenotype := "Het"]
#     dat[allele1==1 & allele2==0, estGenotype := "Het"]
#     dat[allele1=="." | allele2==".", estGenotype := NA]

#     dat[refCount >= 6 & altCount == 0, trueGenotype := "HomRef"]
#     dat[refCount == 0 & altCount >= 6, trueGenotype := "HomAlt"]
#     dat[refCount >= 3 & altCount >= 3 , trueGenotype := "Het"]
#     dat <- dat[!is.na(trueGenotype) & ! is.na(estGenotype)]

#     return(dat[, list("n_recomb"=n_recomb, "N_Sites"=.N, "accurate_sites"=sum(estGenotype==trueGenotype)), by=list(contig, mask, ind)])
# }

# # make table after masking

# o.masked <- o[mask==TRUE]

# o.kept <- o[mask==FALSE]
# o.kept[, mask := NULL]
# col_order <- colnames(o.kept)


# o.kept.combined <- o.kept[, list("contig"="combined", N_Sites=sum(N_Sites), accurate_sites=sum(accurate_sites)), by=list(ind)]
# setcolorder(o.kept.combined, col_order)
# o.kept <- rbindlist(list(o.kept, o.kept.combined))

# o.nofilter.combined <- o[, list("contig"="combined", N_Sites=sum(N_Sites), accurate_sites=sum(accurate_sites)), by=list(ind)]
# o.nofilter <- o[, list( N_Sites=sum(N_Sites), accurate_sites=sum(accurate_sites)), by=list(ind, contig)]
# setcolorder(o.nofilter, col_order)
# setcolorder(o.nofilter.combined, col_order)
# o.nofilter <- rbindlist(list(o.nofilter.combined, o.nofilter))
# o.nofilter[, mask := "Unmasked"]
# o.kept[, mask := "With Masking"]

# dat <- rbindlist(list(o.nofilter, o.kept))
dat[, acc := accurate_sites / N_Sites]
dat[, ind := factor(ind)]
dat[, contig := factor(contig, levels=c("2L","2R","3L","3R","X", "combined"))]

ggplot(dat, aes(x=mask, y=acc, color=ind, group=ind)) + geom_point() + geom_line() + facet_grid(.~contig) +
labs(x="", y=

