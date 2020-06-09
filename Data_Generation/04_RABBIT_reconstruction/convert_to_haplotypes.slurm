#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# This script converts RABBIT output into haplotype map format (which can then be converted
# to site-by-site genotypes via tabix, in convert_to_genotypes.slurm

population_stem="${1}"
chromosome="${2}"
first_ind=${3}
last_ind=${4}
generation="${SLURM_ARRAY_TASK_ID}"

topDirectory="${PWD}"

population="${population_stem}_F${generation}"
#haplotypesFile="${population}/${chromosome}.RABBIT.vcf"

mkdir -p "${population}/imputed_haplotypes/"

for ind_id in $(seq ${first_ind} ${last_ind}); do
  
rabbit_csv="${population}/RABBIT_output/${ind_id}.${chromosome}.RABBIT.out.csv"
haplotypes_out_file="${population}/imputed_haplotypes/${ind_id}.${chromosome}.estimate.haps"
echo -e "chromosome\tstart\tstop\tpar1\tpar2" \
  > "${haplotypes_out_file}"

# if nonrecombinant (i.e. no rabbit-generated output, only manuala output),
# parse manual output to make nonrecombinant haplotype map file
nLines=$(wc -l $rabbit_csv | awk '{print $1}')
if [[ ${nLines} -eq 1 ]]; then
    par=$(awk '{print $3}' $rabbit_csv)
    chromosome=$(echo $rabbit_csv | cut -d "." -f 2)
    echo -e "${chromosome}\t1\t25000000\t${par}\t${par}" >> "${haplotypes_out_file}"
else
python - ${rabbit_csv} <<EOF >> "${haplotypes_out_file}"
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
