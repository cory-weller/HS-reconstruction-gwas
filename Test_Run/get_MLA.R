#!/usr/bin/env Rscript
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

sample_id <- args[1]
chromosome <- args[2]

freqs_filename <- paste(sample_id, ".", chromosome, ".freqs", sep="")
harp_csv <- paste(chromosome, ".harp.csv", sep="")

# Read in header of harp csv file (for column header values)
con <- file(harp_csv,"r")
first_line <- readLines(con,n=1)
close(con)
csv_col_names <- unlist(strsplit(first_line, ","))
line_IDs <- csv_col_names[3:(length(csv_col_names)-1)]


harp_freqs <- fread(input = freqs_filename,
                    header = FALSE,
                    showProgress = FALSE,
                    na.strings = "-nan")



setnames(harp_freqs, c("chromosome","start","stop", line_IDs))

# Subset to only include known founders

harp_freqs.long <- melt(harp_freqs, measure.vars = colnames(harp_freqs)[4:length(colnames(harp_freqs))], variable="lineID", value="freq")
harp_freqs.long[, q99 := quantile(freq, 0.99, na.rm=TRUE), by=chromosome]

mla <- harp_freqs.long[, .N, by=list(chromosome, lineID, freq >= q99)][freq == TRUE][order(chromosome,-N)][,c("chromosome","lineID","N")]

fwrite(mla, file=paste(sample_id, ".", chromosome, ".mla", sep=""),
            sep="\t",
            quote=F,
            row.names=F,
            col.names=T
      )
