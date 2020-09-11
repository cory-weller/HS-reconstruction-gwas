#!/usr/bin/env py

import sys

sites_vcf = sys.argv[1]

alleles = {
'0/0' : '0',
'1/1' : '1',
'./.' : 'NA',
}

previous_chromosome = 'BEGIN_FILE'
output = ''
with open(sites_vcf, 'r') as infile:
    for line in infile:
        if not line.startswith("##"):
            splitline = line.rstrip().split()
            if line.startswith("#"):
                line_ids = splitline[9:]
                n_lines = len(line_ids)
                continue
            else:
                chromosome = splitline[0]
                position = splitline[1]
                genotypes = [alleles[x] for x in splitline[9:]]
            if previous_chromosome == 'BEGIN_FILE' :
                previous_chromosome = chromosome
                outfile = open('''%s.RABBIT.csv''' % (chromosome), 'w')
                print("Working on %s.RABBIT.csv" % (chromosome))
                header = ','.join([str(x) for x in (["#CHROM", 'POS'] + line_ids)]) + '\n'
                outfile.write(header)
            if previous_chromosome == chromosome:
                vals = [chromosome, position] + genotypes
                row = ','.join([str(x) for x in vals]) + '\n'
                outfile.write(row)
            else:
                outfile.close()
                previous_chromosome = chromosome
                outfile = open('''%s.RABBIT.csv''' % (chromosome), 'w')
                print("Working on %s.RABBIT.csv" % (chromosome))
                header = ','.join([str(x) for x in (["#CHROM", 'POS'] + line_ids)]) + '\n'
                outfile.write(header)
    outfile.close()

