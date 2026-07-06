rm(list = ls())

library(data.table)
library(tximport)
library(readxl)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(xlsx)
library(biomaRt)
library(EnhancedVolcano)
library(Cairo)
library(ggrepel)
library(openxlsx)
library(pheatmap)

###### 1. input_RSEM data ######
x.dir = "/path/to/your/rsem/data"
x.files = Sys.glob(paste0(x.dir, "/*.genes.results"))


x.outdir = "/path/to/your/output/directory"
dir.create(x.outdir)

##### patient information
patient.df = read_excel()
patient.df = as.data.frame(patient.df)


##### gene symbols using "GENCODE/v41/gencode.v41.annotation.gtf"
gtf= rtracklayer::import(paste0("/path/to/your/GENCODE/GTF/file")
gtf.df = as.data.frame(gtf)
'%ni%' = Negate('%in%')
gtf.df = subset(gtf.df, seqnames %ni% c("chrX", "chrY", "chrM"))

symbol.df = gtf.df[,c("gene_name","gene_id")]
symbol.df = symbol.df[!duplicated(symbol.df$gene_id),]
rownames(symbol.df) = symbol.df$gene_id
symbol.df$gene_id = NULL
colnames(symbol.df) = "SYMBOL"

###### import RSEM data
rsem.df = tximport(x.files, type = "rsem", txIn = FALSE, txOut = FALSE)
rsem.df = as.data.frame(rsem.df$counts)
rsem.df <- na.omit(rsem.df)
##### filter genes
rsem.df <- rsem.df[!rowSums(rsem.df) == 0,]

patient_ids = gsub("RNA_", "", basename(x.files))
patient_ids = gsub(".genes.results", "", patient_ids)
patient_ids = sapply(strsplit(patient_ids ,"-"), `[`, 1)
colnames(rsem.df) = patient_ids

rsem.df = merge(rsem.df, symbol.df, by = "row.names")
rownames(rsem.df) = rsem.df$Row.names
rsem.df$Row.names = NULL


########## mean value of duplicated genes
unique.ids = unique(rsem.df$SYMBOL)
f.x = function(x){
	unique.df = subset(rsem.df, SYMBOL == x)
	mean.row = as.data.frame(t(colMeans(unique.df[, -ncol(rsem.df)])))
	rownames(mean.row) = x
	return(mean.row)
}
rsem.df = do.call(rbind, lapply(unique.ids, f.x))


###### make coldata and run DESEQ2 for each phenotypic condition
conditions_list = conditions_list <- c("Phenotype_A", "Phenotype_B", "Phenotype_C")


for (i in 1:length(conditions_list)) {
  condition_name = conditions_list[i]
  print(condition_name)
  condition.outdir = paste0(x.outdir, "/", condition_name)
  dir.create(condition.outdir)
  
  ###### 2. input_coldata ######
  coldata = patient.df[,c("Patient_ID", condition_name)]
  coldata = coldata[!grepl("NA", coldata[,c(condition_name)]),]
  
  
  coldata$condition[coldata[,c(condition_name)] == "NO"] <- "0_control"
  coldata$condition[coldata[,c(condition_name)] == "YES"] <- "1_case"
  coldata <- na.omit(coldata)
  
  
	### run deseq2
	cts <- rsem.df[,coldata$Patient_ID]
	cts <- cts[rowSums(is.na(cts)) != ncol(cts), ]
	cts = round(cts)
	
	
	rownames(coldata) = coldata$Patient_ID
	coldata$Patient_ID = NULL
	cts <- cts[, rownames(coldata)]
  all(rownames(coldata) == colnames(cts))
	

	dds <- DESeqDataSetFromMatrix(countData = cts,
	                              colData = coldata,
	                              design= ~ condition)
               
	dds <- DESeq(dds)
	res <- as.data.frame(results(dds))
	res = res[!is.na(res$padj),]

	openxlsx::write.xlsx(res, paste0(condition.outdir,"/DEG_", condition_name,"_unfiltered_genes.xlsx"), row.names = TRUE)
}




