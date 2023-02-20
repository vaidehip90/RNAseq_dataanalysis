# Installing packages and loading the packages needed for RNAseq data anlysis
install.packages("R.utils")
library(R.utils)
BiocManager::install("Rsubread")
library(Rsubread)

install.packages("data.table")
library(data.table)

BiocManager::install("RUVSeq")
library(RUVSeq)

BiocManager::install("DESeq2")
library(DESeq2)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

install.packages("pheatmap")
library(pheatmap)

install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("ggplot2")
library(ggplot2)

BiocManager::install("Rqc")
library(Rqc)


#downloading the file 
dataf <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143630/suppl/GSE143630_RCC_htseq_counts.txt.gz"
file<-"GSE143630_RCC_htseq_counts.txt.gz"
download.file(dataf,file)

#Read the file

candat <- read.csv("GSE143630_RCC_htseq_counts.txt.gz",sep = "",row.names = 1)

#count how many patients have stage 1 and stage 2 cancer
count_of_T1=length(grep(x=colnames(candat),pattern = "^T1."))
count_of_T2=length(grep(x=colnames(candat),pattern = "^T2."))

#filtering the zeros 

filter <- apply(candat,1,function(x) length(x[x>0])>=2)
filtered <- candat[filter,]
genes <-rownames(filtered)

phase1 <- rep(c("T1"),each=count_of_T1)
phase2 <- rep(c("T2"),each=count_of_T2)
x<-c(phase1,phase2)
x<-as.factor(x)

set <- newSeqExpressionSet(as.matrix(filtered),phenoData = data.frame(x,row.names = colnames(filtered) ))
set

#setting the color
colors <- brewer.pal(3,"Set2")

pdf(file="allplot.pdf")
# Number of reads across the samples
plotRLE(set,outline=FALSE, ylim=c(-4,4),col=colors[x]) + title("Number of reads across the samples") 


#Normalization
set <- betweenLaneNormalization(set, which="upper")
plotRLE(set,outline=FALSE, ylim=c(-4,4),col=colors[x]) + title("Plot after Normalization")

differentgroups <- makeGroups(x)
set3 <- RUVs(set, genes,k=1,differentgroups)

#differential Expression

Diff <- DESeqDataSetFromMatrix(countData = counts(set3), colData = pData(set3),design = ~ W_1 + x)
Diff <- DESeq(Diff)
diffres <- results(Diff)
head(diffres)

write.table(diffres, file = "RUVseq.csv", sep = ",",col.names = NA, qmethod = "double")
dfdata<- read.csv("RUVseq.csv")

# statiscal significant
dfdata <- subset(dfdata,dfdata$padj <= 0.05)

#creating a volcano plot

EnhancedVolcano(diffres, lab = rownames(diffres),x = 'log2FoldChange',y= 'pvalue')

#creating a heatmap
#Top 10 genes
candat1 <- candat[1:10,]
heatmap(candat1, cluster_cols = FALSE)

#box plot for phase1 and phase2

#Transpose first
candat <- t(candat)
candat <- as.data.frame(candat)

Tumorstages <- c(rep(c("T1"),each=count_of_T1), rep(c("T2"),each=count_of_T2))
candat <- cbind(Tumorstages,candat)                 

gg <- ggplot(candat, aes(x=Tumorstages, y=ABAT, fill=Tumorstages)) +geom_boxplot() + geom_jitter(shape=16,position=position_jitter(0.01)) +ggtitle("gene-ABAT for phase 1 and phase 2")                 
gg

dev.off()
