###### MethylKit ######
library(data.table)
library(dplyr)
library(methylKit)
library(gplots) 
library(genomation)
library(GenomicRanges)
library(pheatmap)
library(GenomicRanges)
library(ggplot2)
library(readxl)

#### 1. calculate all site beta value ####
#### import processBismarkAln output files for all samples
file.list = Sys.glob("/path/to/your/bismark/output/sample_name/*_CpG.txt")

x.outdir = "/path/to/your/output/directory"
dir.create(x.outdir)

sample.list = c(basename(dirname(file.list)))

####
type = rep(1, each=length(file.list))

myobj <- methRead(as.list(file.list), sample.id = as.list(basename(dirname(file.list))), treatment= as.numeric(type), assembly="hg38", context = "CpG", resolution = "base")
meth=methylKit::unite(myobj)
########### exclude XY chromosomes 
patterns = paste0("chr", seq(1:22))
meth = meth[meth$chr %in% patterns,]
meth.df <- getData(meth)

## 1-1. read count ##
meth.colname = c("chr", "start", "end", "strand")
for (sample.id in sample.list) {
  sample.colnames <- c(paste0("coverage_", sample.id), paste0("numCs_", sample.id), paste0("numTs_", sample.id))
  meth.colname <- c(meth.colname, sample.colnames)
}
colnames(meth.df) <- meth.colname


out.file = paste0(x.outdir, "/MethylKit.AllCpGSite.readcount.txt")
write.table(meth.df, out.file, sep = "\t" , quote = F, row.names = F)

## 1-2. betavalue ##
percM=percMethylation(meth)

percM.df <- data.frame(percM)
beta.df <- percM.df/100
colnames(beta.df) <- c(basename(dirname(file.list)))

out.file = paste0(x.outdir, "/MethylKit.AllCpGSite.betavalue.txt")
write.table(beta.df, out.file, sep = "\t" , quote = F, row.names = F)

## 1-3. read count + betavalue ##
colnames(beta.df) <- paste0("beta_", colnames(beta.df))
meth.beta.df <- cbind(meth.df, beta.df)

out.file = paste0(x.outdir, "/MethylKit.AllCpGSite.readcount.betavalue.txt")
write.table(meth.beta.df, out.file, sep = "\t" , quote = F, row.names = F)

#### 2. gene annotation ####
meth.file = paste0(x.outdir, "/MethylKit.AllCpGSite.readcount.betavalue.txt")
meth.df <- fread(meth.file, header = T, data.table = F)

##### gene symbols using pre-made promoter +/-1000bp bed file: GENCODE/v41/gencode.v41.promoter_region.genelevel.1000.bed
bed.file = "/path/to/your/BED/file"
bed.df <- fread(bed.file, header = T, data.table = F)
colnames(bed.df) = c("chr", "prom.start", "prom.end", "gene_name", "idk", "strand")

meth.obj <- makeGRangesFromDataFrame(meth.df, start.field="start", end.field="end", strand.field="strand")
bed.obj <- makeGRangesFromDataFrame(bed.df, start.field="prom.start", end.field="prom.end", strand.field="strand")

overlap.idx.df <- data.frame(findOverlaps(meth.obj, bed.obj, ignore.strand = TRUE))
overlap.df <- cbind(meth.df[overlap.idx.df$queryHits,], bed.df[overlap.idx.df$subjectHits,c("gene_name", "prom.start", "prom.end", "strand")])
rownames(overlap.df) <- NULL

out.file = paste0(x.outdir, "/MethylKit.CpGSite.prom_1000.readcount.betavalue.txt")
write.table(overlap.df, out.file, sep = "\t" , quote = F, row.names = F)

#### 3. calculate mean beta value ####
overlap.file = paste0(x.outdir, "/MethylKit.CpGSite.prom_1000.readcount.betavalue.txt")
overlap.df <- fread(overlap.file, header = T, data.table = F)
colnames(overlap.df)[which(colnames(overlap.df) == "gene_name")] <- "SYMBOL"

merge_by = "SYMBOL" ## c("SYMBOL", "TranscripId") select one 
ids <- unique(overlap.df[,merge_by])

sum.df = data.frame()
for (i in 1:length(ids)) {
  id = ids[i]
  if(i %% 1000==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", i, "/", length(ids),"\n"))
  }
  sub.df <- subset(overlap.df, overlap.df[,merge_by] == id)
  
  sub.coverage.df <- unique(sub.df[,c(grep("coverage", colnames(sub.df)),grep("num", colnames(sub.df)))])
  sum.coverage.row <- data.frame(matrix(colSums(sub.coverage.df), nrow = 1))
  colnames(sum.coverage.row) <- paste0("sum_", colnames(sub.coverage.df))
  
  pos.list <- c(sub.df$start, sub.df$end)
  
  chr.row <- data.frame(chr = unique(sub.df$chr), start = min(pos.list), end = max(pos.list), strand = "+" , by = id)
  colnames(chr.row)[colnames(chr.row) == "by"] <- merge_by
  sub.row <- cbind(chr.row, sum.coverage.row)
  
  sum.df <- rbind(sum.df, sub.row)
}

meth.gene.mean.df = sum.df[,c("chr", "start", "end", "strand", merge_by)]
for (sample.id in sample.list) {
  sub.df <- sum.df[,grep(sample.id, colnames(sum.df))]
  sub.df$meanbeta <- sub.df[,grep("numCs", colnames(sub.df))]/sub.df[,grep("coverage", colnames(sub.df))]
  colnames(sub.df)[colnames(sub.df) == "meanbeta"] <- paste0("meanbeta_", sample.id)
  
  meth.gene.mean.df <- cbind(meth.gene.mean.df, sub.df)
}

out.file = paste0(x.outdir, "/MethylKit.CpGSite.prom_1000.genelevel.sum_readcount.mean_betavalue.txt")
write.table(meth.gene.mean.df, out.file, sep = "\t" , quote = F, row.names = F)
















