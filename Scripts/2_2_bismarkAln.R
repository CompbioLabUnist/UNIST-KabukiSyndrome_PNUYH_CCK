#!/usr/bin/Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(methylKit))

option_list <- list(
	make_option(c("--input"), action="store", type="character", help="1"),
	make_option(c("--sample_id"), action="store", type="character", help="2"),
	make_option(c("--output"), action="store", type="character", help="3")
)

opt <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=T)

input <- opt$input
sample_id <- opt$sample_id
output <- opt$output

my.methRaw = processBismarkAln( location = input, sample.id=sample_id, assembly="hg38", read.context="CpG", save.folder=output)
