liability <- function(x) {1/(1+exp(-1*x))}

# create dummy data set
# dd <- data.table(x=0.5)
# g1 <- ggplot(dd) + stat_function(fun=liability) + xlim(-5,5) + ylim(0,1) +
# labs(x="Liability / Risk score", y="Probability of 'case' assignment") +
# theme_few(16)

# single locus effect function, based on allele frequency
sle <- function(AF, effect_size) {
  midpoint <- effect_size - (effect_size * 2 *AF)
  rnorm(1, mean=midpoint, sd=0.005)
}

tabix_read <- function(DT) {
  # input of data.table with columns CHR and POS
  o <- foreach(CHR = DT[[1]], POS=DT[[2]], .combine="rbind") %do% {
      tabix_cmd <- paste("tabix -h ../../input_data/haplotypes.vcf.gz ", CHR, ":", POS, "-", POS, sep="")
      dat <- fread(cmd=tabix_cmd, sep="\t")
      return(dat)
  }
  return(o)
}


fillVectorNAs <- function(x) {
  na.omit(x)[cumsum(!is.na(x))]
}

fillDataTableNAs <- function(DT, cols) {
  DT[, (cols) := lapply(.SD, function(x) fillVectorNAs(x)), .SDcols=cols]
}

convert_to_wide <- function(DT) {
  # Collapse redundant breakpoints
  DT[, lineID_rleid := rleid(chromosome, haplotype, ind, lineID)]
  haps.collapsed <- DT[, list(chromosome,haplotype,start=min(start), stop=max(stop), ind, sex, lineID), by=lineID_rleid]
  haps.collapsed[,lineID_rleid := NULL]
  DT <- haps.collapsed[!duplicated(haps.collapsed)]
  

    unique_stops <- unique((sort(DT[,stop])))
    if(length(unique_stops) == 1) {
      unique_starts <- 1
    } else {
      unique_starts <- c(0,unique_stops[1:(length(unique_stops)-1)]) +1
    }
    
    DT.out <- data.table(start=unique_starts, stop=unique_stops)

    setkey(DT.out, start)
    setkey(DT, start)

    DT.out <- dcast(DT.out[DT], start+stop~haplotype, value.var="lineID")

    fillDataTableNAs(DT.out, c("1","2"))

    setcolorder(DT.out, c("start","stop","1","2"))
    setnames(DT.out, c("1","2"), c("par1","par2"))
    return(DT.out[])
}

load_haplotypes <- function(population, RIL) {
  haplotypes_file <- paste("../../01_forward_simulator/", population, ".haps", sep="")
  haplotypes_file_wide <- paste("../", population, ".haps.wide", sep="")
  if(! file.exists(haplotypes_file_wide)) {
    # generate haplotypes file
    cat(paste("building wide haplotypes file for population ", population, "...\n", sep=""))
    all_haplotypes <- fread(haplotypes_file)
    chr <- "2L"


    if(RIL == TRUE) {
      setnames(all_haplotypes, "RIL_ID", "ind")
      o <- foreach(n=1:max(all_haplotypes[,ind]), .combine="rbind") %dopar% {
        print(n)
        ind_wide.dt <- convert_to_wide(all_haplotypes[ind==n & chromosome==chr])
        ind_wide.dt[, "chromosome" := chr]
        ind_wide.dt[, "ind" := n]
        return(ind_wide.dt)
      }
    } else {
      o <- foreach(ind_i = 1:10000, .combine="rbind") %dopar% {
        print(ind_i)
        ind_wide.dt <- convert_to_wide(all_haplotypes[ind==ind_i & chromosome==chr])
        ind_wide.dt[, "chromosome" := chr]
        ind_wide.dt[, "ind" := ind_i]
        return(ind_wide.dt)
      }
    }
    fwrite(o, file=haplotypes_file_wide, quote=F, row.names=F, col.names=T, sep="\t")
    return(o)
  } else {
    cat("loading pre-calculated wide haplotypes file...\n")
    o <- fread(haplotypes_file_wide)
    return(o)
  }
}

load_VCF <- function(chromosomes) {
  vcf <- fread(cmd="zcat ../../input_data/haplotypes.polarized.vcf.gz", skip="CHROM")
  return(vcf[CHROM %in% chromosomes])
}

# Example of relationship between allele freq (AF) and effect size
# o <- foreach(AF=seq(0.01, 1, 0.01), .combine="rbind") %do% {
#   data.table("AF"=AF,
#             "effect"=sapply(1:100, function(x) sle(AF))
#           )
# }
# 
# ggplot(o, mapping=aes(x=factor(AF), y=effect)) + geom_boxplot()

#autosomes <- c("2L","2R","3L","3R")

make_sites_vcf_file <- function(chromosomes) {
  o <- foreach(chr = chromosomes, .combine="rbind") %do% {
    fread(paste("../../input_data/", chr, ".sites", sep=""), header=F, select=c(1,2), col.names=c("CHROM","POS"))
  }
  return(o[])
}

permuteHaplotypes <- function(haps.DT, origFounderIDs, allFounderIDs) {
    haplotypeTranslationDictionary <- data.table("originalID"=origFounderIDs, "newID"=sample(allFounderIDs, size=length(origFounderIDs), replace=FALSE))
    new.haps.DT <- merge(haps.DT, haplotypeTranslationDictionary, by.x="par1", by.y="originalID")
    new.haps.DT[,par1 :=NULL]
    setnames(new.haps.DT, "newID", "par1")
    new.haps.DT <- merge(new.haps.DT, haplotypeTranslationDictionary, by.x="par2", by.y="originalID")
    new.haps.DT[,par2 :=NULL]
    setnames(new.haps.DT, "newID", "par2")
    setkey(new.haps.DT, ind, chromosome, start, stop)
    setcolorder(new.haps.DT, c("start", "stop", "par1", "par2", "chromosome", "ind"))
    return(new.haps.DT[])
}

get_dosage_at_causative_loci <- function(haps.DT, causative.DT) {
  o2 <- foreach(ind_n=1:max(haps.DT[,ind]), .combine="rbind") %dopar% {
    print(ind_n)
    haplotypes <- haps.DT[ind==ind_n]
    foreach(i=1:nrow(haplotypes), .combine="rbind", .errorhandling="remove") %do% {
      o <- causative.DT[CHROM==haplotypes[i,chromosome] & POS %between% c(haplotypes[i,start], haplotypes[i,stop])][,c("CHROM", "POS", haplotypes[i,par1], haplotypes[i,par2]), with=F]
      setnames(o, c("CHROM","POS","hap1","hap2"))
      o[, "ind_n" := ind_n ]
      return(o)
    }
  }

  #o2[, altDosage := apply(.SD, 1, function(x) sum(x=="1/1", na.rm=T)), .SDcols=c("hap1","hap2")]
  #dosage <- merge(o2, causative.DT[,c("CHROM","POS","effect")], by=c("CHROM","POS"))
  #dosage <- dosage[, list("liability"=sum(effect*altDosage)), by=list(ind_n)]
  #dosage[, risk := liability(liability)]
  return(o2)
}


calcFreqs <- function(dat.haps, dat.vcf) {

    # Find every recombination breakpoint throughout all maps
        bps <- sort(unique(dat.haps[,stop]))
        highest.bp <- max(dat.haps[,stop])
        bps <- bps[1:length(bps)-1]

    # Convert all recombination breakpoints to a single table of base pair ranges
        if(length(bps) == 0) {
          segments.to.iter <- data.table(start=c(1), stop=c(highest.bp), containsSNP=c(TRUE))
        } else {
          segments.to.iter <- data.table(start=c(1, bps+1), stop=c(bps, highest.bp))
          # For faster computation, determine which ranges actually contain a SNP,
          # So that ranges without SNPs can be excluded downstream
          dat.sites <- copy(dat.vcf[,.(POS)])
          dat.sites[, nearest := POS]
          setkey(dat.sites, POS)
          setkey(segments.to.iter, start)
          segments.to.iter <- dat.sites[segments.to.iter, list(start, stop, nearest), roll=-Inf]
          segments.to.iter[, containsSNP := ifelse(nearest >= start & nearest <= stop, TRUE, FALSE)]
        }
        setkey(dat.haps, start, stop)

    # While counting the number of times each haplotype is present in the population
        range.haplotype.counts.long <- foreach(start.i = segments.to.iter[containsSNP==TRUE][,start], stop.i=segments.to.iter[containsSNP==TRUE][,stop], .combine="rbind") %do% {
          dat.tmp <- dat.haps[start <= start.i & stop >= stop.i]
          dat.tobind <- rbindlist(list(dat.tmp[, .N, by=list("par"=par1)], dat.tmp[, .N, by=list("par"=par2)]))[, sum(N), by=par]
          dat.tobind[, "start" := start.i]
          dat.tobind[, "stop" := stop.i]
          return(dat.tobind)
        }
        setnames(range.haplotype.counts.long, "V1", "N")  
        setnames(range.haplotype.counts.long, "par", "lineID")  
    # Convert to wide format to mimic organization of the founder VCF file
        range.haplotype.counts <- dcast(range.haplotype.counts.long, start+stop ~lineID, value.var="N", fun.aggregate=sum, fill=0)

    # Determine which columns (founder IDs) have 0 representation within the population
        missing.cols <- allFounderIDs[!(allFounderIDs %in% unique(range.haplotype.counts.long[,lineID]))]

    # Add columns, fill with frequency of 0
        if(length(missing.cols) != 0) {
            range.haplotype.counts[, (missing.cols) := 0]
        }

    # Rearrange columns in frequency same as the parental VCF
        setcolorder(range.haplotype.counts, c("start", "stop", allFounderIDs))

# Fast count snps in each window
    snps <- dat.vcf[,.(POS)]
    snps[, segment := cut(POS, breaks=unique(c(0,segments.to.iter[containsSNP==TRUE][,stop],+Inf)))]
    snps[, segment := as.numeric(as.factor(segment))]

    row.nums <- snps[, .N, by=segment][, list(N=rep(segment, N))][,N]

    snp.haplotype.counts <- range.haplotype.counts[row.nums]

    dat.finalCounts <- dat.vcf[,allFounderIDs, with=F] * snp.haplotype.counts[,allFounderIDs, with=F]

    dat.finalCounts[, sumRef := apply(.SD, 1, function(x) sum(x[(x>0)])), .SDcols=allFounderIDs]

    dat.finalCounts[, sumAlt := apply(.SD, 1, function(x) -1*sum(x[(x<0)])), .SDcols=allFounderIDs]

    dat.finalCounts[, POS := dat.vcf[,POS]]
    dat.finalCounts[, "CHROM" := dat.vcf[,"CHROM"]]

    dat.finalCounts[, c("CHROM","POS","sumRef","sumAlt")][]
}


get.lambdaGC1000 <- function(chisq.array, N_case_individuals, N_control_individuals) {
        expected.median.chisq <- qchisq(0.5, df=1)
        unadjusted.lambda <- median(chisq.array) / expected.median.chisq
    lambda.1000 <- 1 + (unadjusted.lambda-1) * (1/N_case_individuals + 1/N_control_individuals)/(1/1000 + 1/1000)
        return(lambda.1000)
}
