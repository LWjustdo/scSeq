#!/bin/Rscript 
#################################################################################
#          Author: liang.wan, 2019-11-22                                        #
#          Contact: lwan@dragondropinn.com                                      #
#################################################################################

##################### Arguments Input #################
rm(list=ls())
#library(optparse)
args=commandArgs(T)
argNum <- length(args)

if(argNum <= 6){
  cat("\nArguments is not right!\n\n")
  cat("Usage: Rscript seurat_multiSample.R <arg_workDir> <prefix> <arg_minGene> <arg_nGene> <arg_mito> <arg_check> <infile_1> <name_1> <infile_2> <name_2> ...\n\n")
  cat("=============================================\n")
  cat("Note: Maximum of 5 groups were analyzed together!\n")
  cat("=============================================\n")
  cat("Contact:lwan@dragondropinn.com, 2019-03-13\n")
  quit()
}
workDir = args[1]
prefix = args[2]
arg_minGene = as.numeric(args[3])
arg_nGene = as.numeric(args[4])
arg_mito = as.numeric(args[5])
arg_check = args[6]
################## Load packages ###################
setwd(workDir)
library(tidyverse)
library(cowplot)
library(Matrix)
library(Seurat)
#setup ppi and log file
ppi <- 300
sink(paste(prefix, ".seuratAnalysis4Muti.txt", sep=""), append = T)
################## Handle different group input ####
if(argNum == 8){ #1 group to comparison 
  sample_1 = args[7]; name_1 = args[8]
  data_1 <- read.table(file = sample_1, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_1) <- paste(name_1, colnames(data_1), sep = "_")
  mergeData <- data_1
}else if(argNum == 10){ #2 groups to comparison
  sample_1 = args[7]; name_1 = args[8]; sample_2 = args[9]; name_2 = args[10]
  data_1 <- read.table(file = sample_1, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_1) <- paste(name_1, colnames(data_1), sep = "_")
  data_2 <- read.table(file = sample_2, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_2) <- paste(name_2, colnames(data_2), sep = "_")
  mergeData <- cbind(data_1, data_2)
}else if(argNum == 12){ #3 groups to comparison
  sample_1 = args[7]; name_1 = args[8]; sample_2 = args[9]; name_2 = args[10]; sample_3 = args[11]; name_3 = args[12]
  data_1 <- read.table(file = sample_1, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_1) <- paste(name_1, colnames(data_1), sep = "_")
  data_2 <- read.table(file = sample_2, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_2) <- paste(name_2, colnames(data_2), sep = "_")
  data_3 <- read.table(file = sample_3, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_3) <- paste(name_3, colnames(data_3), sep = "_")
  mergeData <- cbind(data_1, data_2, data_3)
}else if(argNum == 14){ #4 groups to comparison
  sample_1 = args[7]; name_1 = args[8]; sample_2 = args[9]; name_2 = args[10]; sample_3 = args[11]; name_3 = args[12]; sample_4 = args[13]; name_4 = args[14]
  data_1 <- read.table(file = sample_1, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_1) <- paste(name_1, colnames(data_1), sep = "_")
  data_2 <- read.table(file = sample_2, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_2) <- paste(name_2, colnames(data_2), sep = "_")
  data_3 <- read.table(file = sample_3, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_3) <- paste(name_3, colnames(data_3), sep = "_")
  data_4 <- read.table(file = sample_4, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_4) <- paste(name_4, colnames(data_4), sep = "_")
  mergeData <- cbind(data_1, data_2, data_3, data_4)
}else if(argNum == 16){ #5 groups to comparison
  sample_1 = args[7]; name_1 = args[8]; sample_2 = args[9]; name_2 = args[10]; sample_3 = args[11]; name_3 = args[12]; sample_4 = args[13]; name_4 = args[14]; sample_5 = args[15]; name_5 = args[16]
  data_1 <- read.table(file = sample_1, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_1) <- paste(name_1, colnames(data_1), sep = "_")
  data_2 <- read.table(file = sample_2, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_2) <- paste(name_2, colnames(data_2), sep = "_")
  data_3 <- read.table(file = sample_3, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_3) <- paste(name_3, colnames(data_3), sep = "_")
  data_4 <- read.table(file = sample_4, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_4) <- paste(name_4, colnames(data_4), sep = "_")
  data_5 <- read.table(file = sample_5, row.names = 1, header = TRUE,sep = "\t")
  colnames(data_5) <- paste(name_5, colnames(data_5), sep = "_")
  mergeData <- cbind(data_1, data_2, data_3, data_4, data_5)
}else{
  cat("\nArgumentis is not right!\n\n")
}

#��������
mergeObj <- CreateSeuratObject(mergeData, project = 'SYDZ', min.cells = 10, min.features = 500,  names.field = 1, names.delim = "_", meta.data= NULL)
#���QC percent
mergeObj[["percent.mt"]] <- PercentageFeatureSet(mergeObj, pattern = "^MT-")
mergeObj[["percent.ribo"]] <- PercentageFeatureSet(mergeObj, pattern = "^RP[SL]")
png(paste(prefix, "_QCPlot.png", sep = ""),width = 15*ppi, height = 10*ppi, res = ppi)
VlnPlot(mergeObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0.1)
dev.off()
#��ͬά�������ȽϿ���������״̬
plot1 <- FeatureScatter(mergeObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mergeObj, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot3 <- FeatureScatter(mergeObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(paste(prefix, "_ScatterPlot.png", sep = ""), width = 15*ppi, height = 3*ppi, res = ppi)
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
dev.off()
#�鿴����֮ǰcell/featurel
mergeObj
#��һ�ι���ϸ��,���ݲ⵽�Ļ�����Ŀ�������庬���������庬����ά�Ƚ��й���
mergeObj <- subset(mergeObj, subset = nFeature_RNA > 2000 & percent.ribo < 40 & percent.mt < 15)
#�鿴����֮��cell/featurel
mergeObj
#���ݱ�׼��
mergeObj <- NormalizeData(mergeObj, normalization.method = "LogNormalize", scale.factor = 10000)
#Ѱ�Ҳ�����������Щ���������η�����fearurel
mergeObj <- FindVariableFeatures(mergeObj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mergeObj), 10)
plot1 <- VariableFeaturePlot(mergeObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(paste(prefix, "_VariableFeaturePlot.png", sep = ""), width = 15*ppi, height = 3*ppi, res = ppi)
CombinePlots(plots = list(plot1, plot2))
dev.off()
#PCA֮ǰ����������
all.genes <- rownames(mergeObj)
mergeObj <- ScaleData(mergeObj, features = all.genes)
#����PCA��ά
mergeObj <- RunPCA(mergeObj, features = VariableFeatures(object = mergeObj))
print(mergeObj[["pca"]], dims = 1:9, nfeatures = 10)
png(paste(prefix, "_PCA.png", sep = ""), width = 12*ppi, height = 15*ppi, res = ppi)
VizDimLoadings(mergeObj, dims = 1:9, reduction = "pca")
dev.off()
#PCA��ά������ӻ�
png(paste(prefix, "_PCACluster.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
DimPlot(mergeObj, reduction = "pca")
dev.off()
#PCA��ͼ
png(paste(prefix, "_PCHeatmap.png", sep = ""), width = 8*ppi, height = 9*ppi, res = ppi)
DimHeatmap(mergeObj, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
#����JackStraw�ز����㷨ȷ�����ݵ�ά��
mergeObj <- JackStraw(mergeObj, num.replicate = 100)
mergeObj <- ScoreJackStraw(mergeObj, dims = 1:20)
png(paste(prefix, "_JackStraw.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
JackStrawPlot(mergeObj, dims = 1:15)
dev.off()
#��������ͼѡ��PCA�������η���
png(paste(prefix, "_PCAElbow.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
ElbowPlot(mergeObj)
dev.off()
#KNN����ϸ����ע�����÷ֱ��ʣ�ֵ�����ӵ��¼�Ⱥ���������ӡ�
#���Ƿ��֣����˲���������0.4-1.2֮��ͨ����Ϊ��Լ3K����Ԫ�ĵ�Ԫ���ݼ��������õĽ����
#���ڽϴ�����ݼ�����ѷֱ���ͨ�������ӡ�
mergeObj <- FindNeighbors(mergeObj, dims = 1:15)
mergeObj <- FindClusters(mergeObj, resolution = 0.5)
#UMAP�����Խ�ά������ӻ�
mergeObj <- RunUMAP(mergeObj, dims = 1:15)
png(paste(prefix, "_UMAPCluster.png", sep = ""), width = 12*ppi, height = 4*ppi, res = ppi)
p1 <- DimPlot(mergeObj, pt.size = 0.1, group.by = "orig.ident")
p2 <- DimPlot(mergeObj, pt.size = 0.1, reduction = "umap")
CombinePlots(plots = list(p1, p2))
dev.off()
#��������ܼ�����ļ�
saveRDS(mergeObj, file = paste0(prefix, ".RDS"))
#Ѱ��marker Gene��ע��min.pct����������Ԫ���κ�һ��������С�ٷֱȼ������,logfc.thresholdΪ������ﱶ����ֵ
#���һ��cluster
#cluster1.markers <- FindMarkers(mergObj, ident.1 = 1, min.pct = 0.25)
#�������cluste֮��Ƚ�
#cluster5.markers <- FindMarkers(mergObj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
mergeObj.markers <- FindAllMarkers(mergeObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mergeObj.markers, file = paste(prefix, "_markerGene.csv", sep = ""),  row.names = F)
#��Top10��marker gene������ͼ��FeaturePlot��VlnPlot����չʾ
top10 <- mergeObj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#��ͼ
png(paste(prefix, "_ClusterHeatmap.png", sep = ""), width = 12*ppi, height = 12*ppi, res = ppi)
DoHeatmap(mergeObj, features = top10$gene, size = 5, lines.width = 50)
dev.off()
#Top����VlnPlot
dir.create("./Top10MarkerGene_VlnPlot")
setwd("./Top10MarkerGene_VlnPlot")
for(g in top10$gene){
    png(paste0(g, "_VlnPlot.png"), width = 6*ppi, height = 4*ppi, res = ppi)
    print(VlnPlot(mergeObj, features = c(g), pt.size = 0.1))
    dev.off()
}
#Top����featurePlot
dir.create("../Top10MarkerGene_FeaturePlot")
setwd("../Top10MarkerGene_FeaturePlot")
for(g in top10$gene){
    png(paste0(g, "_FeaturePlot.png"), width = 6*ppi, height = 4*ppi, res = ppi)
    print(FeaturePlot(mergeObj, features = c(g), pt.size = 0.1))
    dev.off()
}
setwd("../")
#write meta Data
write.csv(mergeObj@meta.data, file = paste(prefix, "_metaData.csv", sep = ""),  row.names = T)
#����clusterϸ����Ŀ
write.csv(table(mergeObj@active.ident), file = paste(prefix, "_clusterCellNum.csv", sep = ""), row.names = F)
########################################
#����singleR����ϸ������ע�ͷ���




















