#!/bin/bash
#
#
#
#
# BEGIN parameter definitions

# The parameter file is a tab-delimited file containing a list of individuals to reconstruct within a SLURM job array
parameterFile=${1}
manual_line_entry=${2}

if [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
    if [[ ! -z "$manual_line_entry" ]]; then
        # If both slurm array ID and manual entry exist 
        echo "Error: Both job array ID and manual line number specified"
        echo "Either run a single job with manual line number, or submit a job array spanning desired line numbers to impute."
        exit 1
    elif [[ -z "$manual_line_entry" ]]; then
        # If slurm array ID exists, no manual entry
        lineNumber=${SLURM_ARRAY_TASK_ID}
    fi    
elif [[ -z "$SLURM_ARRAY_TASK_ID" ]]; then
    # if slurm array ID does not exist
    if [[ ! -z "$manual_line_entry" ]]; then
        # if only manual entry exists
        lineNumber=${manual_line_entry}
    elif [[ -z "$manual_line_entry" ]]; then
        # If neither job array nor manual line number exist
        echo "Error: no manual line number specified, and job array not used"
        echo "Either run a single job with manual line number, or submit a job array spanning desired line numbers to impute."
        exit 1
    fi
fi

reconstructionGroup=$(sed "${lineNumber}q;d" $parameterFile | cut -f 1 )
ind_id=$(sed "${lineNumber}q;d" $parameterFile | cut -f 2 )
nGenerations=$(sed "${lineNumber}q;d" $parameterFile | cut -f 3 )

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"
topDirectory=${SCRIPTPATH%scripts}
workDir="${topDirectory}/${reconstructionGroup}/"
gatk=${topDirectory}/etc/GenomeAnalysisTK.jar

cd $workDir
# Check for presence of "chromosomes.txt" for imputing
if [[ -f "chromosomes.txt" ]]; then
    readarray -t chromosomes < chromosomes.txt
    echo "imputing over all chromosomes in file chromosome.txt:"
    echo "${chromosomes[@]}"
else
    # use all chromosomes contained within the VCF file
    chromosomes=($(awk '{ if (!a[$1]++ ) print $1;}' haplotypes.vcf | grep -v "^#"))
    echo "WARNING: imputing over all chromosomes in VCF file"
    echo "This may have unintended effects if your VCF file contains extra chromosomes not being imputed, e.g. mitochrondrial genomes, Wolbachia, Y chromosomes."
    echo "If this is unintended behavior, add a text file "chromosome.txt" containing a list of chromosomes to impute, one per line."
    echo "Chromosomes being imputed, as listed in VCF file:"
    echo "${chromosomes[@]}"
fi
# Get array of chromosomes to reconstruct, based on unique chromosome values within haplotypes.vcf

#RABBIT Parameters
eps="0.005"
epsF="0.005"
RABBITmodel="jointModel"
RABBITestfun="origViterbiDecoding"
RABBITpackageLocation="${topDirectory}/etc/RABBIT/RABBIT_Packages/"
topNfounders=10     #topNfounders indicates that X most likely founders are chosen for performing higher-resolution imputation
maxNSites=5000      # Maximum number of SNPs per chromosome to use for imputation

echo "exit line for dev"
exit 0

# END parameter definitions
#
#
#
#
#
# BEGIN calculate .bam file stats

if [[ ! -f "$ind_id".bam.stats ]]; then
    echo "calculating .bam file stats"
    ${topDirectory}/etc/samtools idxstats ${ind_id}.bam > ${ind_id}.bam.stats
fi

# END calculate .bam file stats
#
#
#
#
#
# BEGIN calculating readcounts

if [[ ! -f $ind_id.readcounts.gz ]]; then
    echo "calculating readcounts for $ind_id"

    java \
    -Xmx6G \
    -jar $gatk \
    -R reference.fasta \
    -T ASEReadCounter \
    -I $ind_id.bam \
    -o $ind_id.readcounts \
    -sites haplotypes.vcf
    gzip $ind_id.readcounts
fi

# END calculating readcounts
#
#
#
#
#
# BEGIN calculate haplotype freqs with HARP

function getFreqs {
    topDirectory=$1
    reconstructionGroup=$2
    ind_id=$3
    chromosome=$4


    workDir=$topDirectory/$reconstructionGroup/
    
    referenceGenome=$workDir/reference.fasta
    
    window_step=100000
    window_width=100000
    
    picard=$topDirectory/etc/picard.jar
    pear=$topDirectory/etc/pear
    bwa=$topDirectory/etc/bwa/bwa
    harp=$topDirectory/etc/harp


# iterate HARP through chromosomes

if [[ ! -f "$workDir/$ind_id.freqs.gz" ]]; then
# if final freqs output of all chromosomes does not exist

    # continue if THIS chromosome is done
    if [[ -f "$workDir/$ind_id.$chromosome.freqs" ]]; then
        echo "chromosome $chromosome already complete"
    else
    echo working on chr $chromosome
        
        
        echo starting harp on chromosome $chromosome
        # change to ramdisk directory for fast read/write of harp operations
        ramdiskDir=/dev/shm/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID/$reconstructionGroup/$ind_id/$chromosome
        mkdir -p $ramdiskDir && cd $ramdiskDir
         
        length=$(awk -v chr=$chromosome '$1==chr {print $2}' $workDir/$ind_id.bam.stats)
        
        # Run harp like
        echo running harp like
        $harp like \
        --bam $workDir/$ind_id.bam \
        --region $chromosome:1-$length \
        --refseq $referenceGenome \
        --snps $workDir/$chromosome.priors.csv \
        --stem $ind_id.$chromosome
        
        # Run harp freq
        echo running harp freq
        $harp freq \
        --bam $workDir/$ind_id.bam \
        --region $chromosome:1-$length \
        --refseq $referenceGenome \
        --snps $workDir/$chromosome.priors.csv \
        --stem $ind_id.$chromosome \
        --window_step $window_step \
        --window_width $window_width \
        --em_min_freq_cutoff 0.0001
        
        mv $ind_id.$chromosome.freqs $workDir && cd $workDir && rm -rf $ramdiskDir
    fi
else
    echo "harp freq file already exists for $chromosome"
fi
}

export -f getFreqs

parallel -j 1 getFreqs ::: $topDirectory ::: $reconstructionGroup ::: $ind_id ::: ${chromosomes[@]}

echo "all chromosomes done for ind $ind_id, cleaning up"

cd $workDir && rm -rf /dev/shm/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
cat $workDir/$ind_id.*.freqs | gzip -c > $workDir/$ind_id.freqs.gz && rm $workDir/$ind_id.*.freqs

# END calculate haplotype freqs with HARP
#
#
#
#
#
# BEGIN calculate most likely parents

module load R/3.3.0
cd ${workDir}

Rscript - <<EOF > $ind_id.mlp
#!/usr/bin/env Rscript

library(data.table)

harp.freqs <- fread("zcat ${ind_id}.freqs.gz", header=FALSE, showProgress=FALSE, na.strings="-nan")
founders.list <- fread('founders.txt', header=FALSE)
setnames(harp.freqs, c("chromosome","start","stop", founders.list[,V1]))

harp.freqs.long <- melt(harp.freqs, measure.vars = colnames(harp.freqs)[4:length(colnames(harp.freqs))], variable="lineID", value="freq")
harp.freqs.long[, q99 := quantile(freq, 0.99, na.rm=TRUE), by=chromosome]

mlp <- harp.freqs.long[, .N, by=list(chromosome, lineID, freq>=q99)][freq==TRUE][order(chromosome,-N)][,c("chromosome","lineID","N")]
write.table(mlp, file="", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# output: $ind_id.mlp
EOF


# END calculate haplotype freqs with HARP
#
#
#
#
#
# BEGIN preparation of harp input file

Rscript - <<EOF

#!/usr/bin/env Rscript

library(data.table)
library(foreach)


# in R, combine information from three sources:
# 1. Readcounts, to find out which SNPs is the whole set we're considering
# 2. population VCF file, for all variants
# 3. alleleFreqs, for information content
# 4. recombination .bed file

# 3. Of remaining SNPs, merge in read depth and allele frequency
# 4. Subset by top N, according to criteria deemed important
# 5. Re-organize by chromosome + POS
# 6. Iterate over columns and print to file (effectively transposing without holding in memory)


# 1. read in population VCF, readcounts, and freqs for this individual

ind_id <- "${ind_id}"
topNfounders <- ${topNfounders}
maxNSites <- ${maxNSites}


founders.list <- fread('founders.txt', header=FALSE, showProgress=FALSE)

vcf <- fread("zcat haplotypes.vcf.gz", skip="#CHROM", header=TRUE, na.strings="./.", showProgress=FALSE)
vcf <- vcf[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", founders.list[,V1]), with=FALSE]
setnames(vcf, "#CHROM", "CHROM")
setkey(vcf, "CHROM", "POS", "ID")

readcounts <- fread(paste("zcat ", ind_id, ".readcounts.gz", sep=""), header=TRUE, showProgress=FALSE)
setkey(readcounts, "contig", "position", "variantID")

mlp <- fread(paste(ind_id, ".mlp", sep=""), header=TRUE, showProgress=FALSE)
mlp[,rank := frank(-N, ties.method="random"), by=chromosome]
mlp.chosen <- mlp[rank <= topNfounders]


# Read in .bed file containing recombination rates
# Rquired column headers are chr (chromosome); start; stop; c (recombination rate, in units of cM/Mb)
bed <- fread("recombination.bed", header=TRUE, showProgress=FALSE)
# Key .bed file by chromosome and window start
setkey(bed, chr, start)

# Generate column for cM (centiMorgans)
bed[, cM := c*((stop-start)/1e6)]

# Calculate cumulative cM at the END of each window
bed[, cumulative_cM := cumsum(cM), by=chr]

# Generate functions (cM ~ BP) to translate base pair position (in VCF) to units of cM
recombination_function <- new.env()

# Create recombination function based on .bed file
for(chromosome in unique(bed[,chr])) {
    recombination_function[[as.character(chromosome)]] <- approxfun(c(0, bed[chr==chromosome][,stop]), c(0,bed[chr==chromosome][,cumulative_cM]))
}

#options(scipen=999)

for(chr.i in unique(bed[,chr])) {
    cat(chr.i)
    cat("\n")
    chosen.founders <- sort(mlp.chosen[chromosome==chr.i][,lineID])
    
    # If no most likely founders for this chromosome, then skip. This should not happen.
    if(length(chosen.founders)==0) {
        cat("No most likely founders for chromosome ")
        cat(chr.i)
        cat("\n")
        next
    }

    # if only one founder was chosen for the entire chromosome (e.g., homozygous for entire chromosome for a single founder)
    if(length(chosen.founders)==1) {
        writeLines(paste("Nonrecombinant Homozygous ", chosen.founders, sep=""), con=paste(ind_id, ".", chr.i, ".RABBIT.out.csv", sep=""))
        next
    }

    vcf.sub <- copy(vcf[CHROM==chr.i][,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", chosen.founders), with=FALSE])
    # if in 0/0 & 1/1 format instead of 0 or 1:
    if(any(vcf.sub=="0/0") | any(vcf.sub=="1/1")) {
        # Convert to factor with level 1 = "0/0", level 2 = "1/1"
        vcf.sub[, (chosen.founders) := lapply(.SD, factor, levels=c("0/0","1/1")), .SDcols=chosen.founders]
        # Convert to numeric
        vcf.sub[, (chosen.founders) := lapply(.SD, as.numeric), .SDcols=chosen.founders]
        # Subtract 1, such that "0/0" is now 0, "1/1" is now 1
        vcf.sub[, (chosen.founders) := lapply(.SD, "-", 1), .SDcols=chosen.founders]
    }
    
    vcf.sub[, nRef := apply(.SD, 1, function(x) sum(x == 0, na.rm=TRUE)), .SDcols=chosen.founders]
    vcf.sub[, nAlt := apply(.SD, 1, function(x) sum(x == 1, na.rm=TRUE)), .SDcols=chosen.founders]
    vcf.sub[, refFreq := nRef/(nRef+nAlt)]
    vcf.sub[, altFreq := nAlt/(nRef+nAlt)]
    setkey(vcf.sub, ID)
    
    # Merge in sequenced allele freq
    vcf.sub.merge <- merge(vcf.sub, readcounts, by.x="ID", by.y="variantID")[refFreq != 1 & altFreq != 1 & otherBases == 0 & improperPairs == 0]
    vcf.sub.merge[, c("totalCount","lowMAPQDepth","lowBaseQDepth","rawDepth","otherBases","improperPairs") := NULL]
    vcf.sub.merge[, "imputed_ind" := ifelse(refCount > 0 & altCount > 0, "12",
                               ifelse(altCount==0, "1N",
                               ifelse(refCount==0, "2N", "NN")))]
                               
    # Select up to maximum number of SNPs
    vcf.sub.merge[,indx := 1:.N]
    
    # Retain ALL heterozygous sites
    vcf.sub.merge[refCount>0 & altCount>0, marker := TRUE]
    vcf.sub.merge[, freq := min(refFreq, altFreq), by=indx]
    nSites <- dim(vcf.sub.merge[is.na(marker)])[1]
    samplePercent <- 0.1
    sampleNSites <- trunc(nSites * samplePercent)
    retainmentPercent <- 0.002
    retainNSites <- trunc(nSites * retainmentPercent)

    

    if(dim(vcf.sub.merge)[1] <= maxNSites) {
        vcf.sub.merge[, marker := TRUE]
    } else {
        while(dim(vcf.sub.merge[marker == TRUE])[1] < maxNSites ) {
            indicesToMark <- vcf.sub.merge[is.na(marker)][sample(.N, size=sampleNSites, replace=TRUE)][order(freq)][1:retainNSites][,indx]
            vcf.sub.merge[indicesToMark, marker := TRUE]
        }
    }

    markers <- vcf.sub.merge[marker==TRUE][sample(.N, size=min(maxNSites, .N), replace=FALSE)]
    setkey(markers, "CHROM", "POS")
    markers[, c("contig", "position", "refAllele", "altAllele", "refCount", "altCount", "freq", "indx", "marker", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "refFreq", "altFreq", "nRef", "nAlt") := NULL] 
    setnames(markers, "CHROM", "Chromosome")
    setnames(markers, "ID", "SNP")
    markers[, cM := recombination_function[[chr.i]](POS)][]
    markers[, POS := NULL]

    markers[, (chosen.founders) := lapply(.SD, "+", 1), .SDcols=chosen.founders]
    line.names=colnames(markers)[3:(length(colnames(markers))-1)]
    nFoundersUsed <- length(colnames(markers))-4

    setcolorder(markers, c("SNP","Chromosome","cM", line.names))
    writeLines(paste("#founders,",nFoundersUsed, sep=""), con=paste(ind_id, ".", chr.i, ".RABBIT.in", sep=""))
    write.table(t(as.matrix(markers)), file=paste(ind_id, ".", chr.i, ".RABBIT.in", sep=""), quote=FALSE, row.names=TRUE, col.names=FALSE, sep=",", na="N", append=TRUE)
}
EOF

# END preparation of harp input files
#
#
#
#
#
# START RABBIT 

module load mathematica
# Generate mathematica script
for filename in $ind_id.*.RABBIT.in; do
chromosome=$(echo $filename | cut -d "." -f 2)
python - <<EOF > $ind_id.$chromosome.RABBIT.m

print """SetDirectory["%s"]""" % "${RABBITpackageLocation}"
print """Needs["MagicReconstruct\`"]"""
print """SetDirectory["%s"]""" % "$(pwd)"
print """popScheme = Table["RM1-E", {%s}]""" % "${nGenerations}"
print 'epsF = %s' % "${epsF}"
print 'eps = %s' % "${eps}"
print 'model = "%s"' % "${RABBITmodel}"
print 'estfun = "%s"' % "${RABBITestfun}"
print 'inputfile = "%s"' % "$(pwd)/${ind_id}.${chromosome}.RABBIT.in"
print 'resultFile = "%s.txt"' % "$(pwd)/${ind_id}.${chromosome}.RABBIT.out"
print """magicReconstruct[inputfile, model, epsF, eps, popScheme, resultFile, HMMMethod -> estfun, PrintTimeElapsed -> True]"""
print 'summaryFile = StringDrop[resultFile, -4] <> ".csv"'
print 'saveAsSummaryMR[resultFile, summaryFile]'
print 'Exit'
EOF

math -noprompt -script $ind_id.$chromosome.RABBIT.m

done

module remove mathematica

# Remove extra files
ls ${ind_id}.*RABBIT.m | xargs rm
ls ${ind_id}.*RABBIT.in | xargs rm
ls ${ind_id}.*RABBIT.out.txt | xargs rm

# END RABBIT
#
#
#
#
#
# BEGIN conversion of RABBIT output to haplotype map

echo -e "chromosome\tstart\tstop\tpar1\tpar2" > "${ind_id}.estimate.haps"

for filename in ${ind_id}.*.out.csv; do

nLines=$(wc -l $filename | awk '{print $1}')
if [[ ${nLines} -eq 1 ]]; then
    par=$(awk '{print $3}' $filename)
    chromosome=$(echo $filename | cut -d "." -f 2)
    echo -e "${chromosome}\t1\t25000000\t${par}\t${par}" >> "${ind_id}.estimate.haps"
else
python - ${filename} <<EOF >> "${ind_id}.estimate.haps"
#!/usr/bin/env python

import re
import sys

filename = sys.argv[1]
ind_id, chromosome = filename.split(".")[0:2]

def getLines(filename):
    mylines = []
    with open(filename, 'r') as infile:
        for line in infile:
            if not line.startswith(('Chromosome', 'cM')):
                mylines.append(line.rstrip())
                if 'ViterbiPath' in line:
                    viterbiLine = len(mylines)
        return [mylines, viterbiLine]

def getDiplotypes(mylines):
    diplotypes = {}
    pattern = re.compile('^diplotype[0-9]+')
    for line in mylines:
        if re.match(pattern, line):
            diplotype, parentCodes, founders = line.split(', ')
            diplotype = diplotype.split('diplotype')[-1]
            founders = founders.split('---')
            diplotypes[int(diplotype)] = founders
    return diplotypes

def getNucleotidePositions(mylines):
    positions = []
    for line in mylines:
        if line.startswith('SNP'):
            #return [x.split('_')[0:2] for x in line.rstrip().split(', ')] if returning Chromosome and position
            return [int(x.split('_')[1]) for x in line.rstrip().split(', ')[1:]] #exclude 1st item in line because it's simply 'SNP'

def getSampleDiplotypes(mylines, viterbiLine):
    paths = {}
    for i in range(viterbiLine, len(mylines)):
        sample, path = mylines[i].rstrip().split(', ')
        path = path.split('-')
        paths[sample] = path
    return paths

def phasePaths(paths):
    phasedPaths = {}
    for path in paths:
        entry = paths[path]
        pairs = []
        for i in range(0,len(entry),2):
            pairs.append(entry[i:i+2])
        for pair in pairs[:-1]:
            pair[1] = [x for x in diplotypes[int(pair[1])]]
        phasedPaths[path] = pairs
    return phasedPaths

def convertToPositions(phasedPaths):
    convertedPaths = {}
    for path in phasedPaths:
        segments = []
        for segment in phasedPaths[path][:-1]:
            segments.append([positions[int(segment[0]) - 1], segment[1]])
        for i in range(0,len(segments)-1):
            segments[i][0] = [int(segments[i][0]), int(segments[i+1][0])-1]
        segments[-1][0] = [int(segments[-1][0]), chromosomeLength]
        convertedPaths[path] = segments
    return convertedPaths

mylines, viterbiLine = getLines(filename)

#viterbiPath = [int(x) for x in mylines[viterbiLine].split()[1].split("-")]

# Get dictionary of all possible diplotype pairs and their numeric code
diplotypes = getDiplotypes(mylines)

# Get nucleotide positions from SNP ID
positions = getNucleotidePositions(mylines)

# Get the viterbi path from RABBIT output
paths = getSampleDiplotypes(mylines, viterbiLine)

phasedPaths = phasePaths(paths)

chromosomeLength = int(positions[-1])

convertedPaths = convertToPositions(phasedPaths)

for i in convertedPaths:
    for j in convertedPaths[i]:
        print '\t'.join([chromosome] + [str(x) for x in j[0]]) + '\t' + '\t'.join([str(x) for x in j[1]])

EOF
fi
done

# END conversion of RABBIT output to haplotype map
#
#
#
#
#
# BEGIN conversion of haplotype map to VCF file

# Set up associative array of lineID : column, to be used by TABIX
vals=($(zgrep -m 1 "^#CHROM" haplotypes.vcf.gz | xargs | cut -f 10-))
declare -A founderIndices

for i in $(seq 0 "${#vals[@]}"); do
    let j=$i+1
    founderIndices[[${vals[$i]}]]=$j
done

# Read through diploid paths, extracting genotypes with TABIX, appending to estimate.vcf file
# Here is where you would change the bgzipped vcf filename to whatever one has all the sites you're trying to pull out based on the path
while read chromosome start stop par1 par2; do
    col1=${founderIndices[[$par1]]}
    col2=${founderIndices[[$par2]]}
    /scratch/$USER/genome-reconstruction/etc/htslib/bin/tabix haplotypes.vcf.gz ${chromosome}:${start}-${stop} | awk -v col1=$col1 -v col2=$col2 '{print $1,$2,$col1,$col2}' >> ${ind_id}.estimate.vcf
done < <(awk 'NR > 1 {print}' ${ind_id}.estimate.haps)
gzip ${ind_id}.estimate.vcf 

ls ${ind_id}.*RABBIT.out.csv | xargs rm

# END conversino of paths to VCF file
#
#
#
#
#
# Finished with individual
