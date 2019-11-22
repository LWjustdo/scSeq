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

#创建对象
mergeObj <- CreateSeuratObject(mergeData, project = 'SYDZ', min.cells = 10, min.features = 500,  names.field = 1, names.delim = "_", meta.data= NULL)
#获得QC percent
mergeObj[["percent.mt"]] <- PercentageFeatureSet(mergeObj, pattern = "^MT-")
mergeObj[["percent.ribo"]] <- PercentageFeatureSet(mergeObj, pattern = "^RP[SL]")
png(paste(prefix, "_QCPlot.png", sep = ""),width = 15*ppi, height = 10*ppi, res = ppi)
VlnPlot(mergeObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0.1)
dev.off()
#不同维度两两比较看样本质量状态
plot1 <- FeatureScatter(mergeObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mergeObj, feature1 = "nCount_RNA", feature2 = "percent.ribo")
plot3 <- FeatureScatter(mergeObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(paste(prefix, "_ScatterPlot.png", sep = ""), width = 15*ppi, height = 3*ppi, res = ppi)
CombinePlots(plots = list(plot1, plot2, plot3), ncol = 3)
dev.off()
#查看过滤之前cell/featurel
mergeObj
#第一次过滤细胞,根据测到的基因数目，线粒体含量，核糖体含量多维度进行过滤
mergeObj <- subset(mergeObj, subset = nFeature_RNA > 2000 & percent.ribo < 40 & percent.mt < 15)
#查看过滤之后cell/featurel
mergeObj
#数据标准化
mergeObj <- NormalizeData(mergeObj, normalization.method = "LogNormalize", scale.factor = 10000)
#寻找差异表达基因，这些基因是下游分析的fearurel
mergeObj <- FindVariableFeatures(mergeObj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mergeObj), 10)
plot1 <- VariableFeaturePlot(mergeObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(paste(prefix, "_VariableFeaturePlot.png", sep = ""), width = 15*ppi, height = 3*ppi, res = ppi)
CombinePlots(plots = list(plot1, plot2))
dev.off()
#PCA之前的缩放数据
all.genes <- rownames(mergeObj)
mergeObj <- ScaleData(mergeObj, features = all.genes)
#数据PCA降维
mergeObj <- RunPCA(mergeObj, features = VariableFeatures(object = mergeObj))
print(mergeObj[["pca"]], dims = 1:9, nfeatures = 10)
png(paste(prefix, "_PCA.png", sep = ""), width = 12*ppi, height = 15*ppi, res = ppi)
VizDimLoadings(mergeObj, dims = 1:9, reduction = "pca")
dev.off()
#PCA降维聚类可视化
png(paste(prefix, "_PCACluster.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
DimPlot(mergeObj, reduction = "pca")
dev.off()
#PCA热图
png(paste(prefix, "_PCHeatmap.png", sep = ""), width = 8*ppi, height = 9*ppi, res = ppi)
DimHeatmap(mergeObj, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()
#运用JackStraw重采样算法确定数据的维度
mergeObj <- JackStraw(mergeObj, num.replicate = 100)
mergeObj <- ScoreJackStraw(mergeObj, dims = 1:20)
png(paste(prefix, "_JackStraw.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
JackStrawPlot(mergeObj, dims = 1:15)
dev.off()
#运用手肘图选择PCA进入下游分析
png(paste(prefix, "_PCAElbow.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
ElbowPlot(mergeObj)
dev.off()
#KNN聚类细胞，注意设置分辨率，值的增加导致集群数量的增加。
#我们发现，将此参数设置在0.4-1.2之间通常可为大约3K个单元的单元数据集返回良好的结果。
#对于较大的数据集，最佳分辨率通常会增加。
mergeObj <- FindNeighbors(mergeObj, dims = 1:15)
mergeObj <- FindClusters(mergeObj, resolution = 0.5)
#UMAP非线性降维聚类可视化
mergeObj <- RunUMAP(mergeObj, dims = 1:15)
png(paste(prefix, "_UMAPCluster.png", sep = ""), width = 12*ppi, height = 4*ppi, res = ppi)
p1 <- DimPlot(mergeObj, pt.size = 0.1, group.by = "orig.ident")
p2 <- DimPlot(mergeObj, pt.size = 0.1, reduction = "umap")
CombinePlots(plots = list(p1, p2))
dev.off()
#保存计算密集结果文件
saveRDS(mergeObj, file = paste0(prefix, ".RDS"))
#寻找marker Gene，注意min.pct是在两组像元的任何一组中以最小百分比检测特征,logfc.threshold为差异表达倍数阈值
#针对一个cluster
#cluster1.markers <- FindMarkers(mergObj, ident.1 = 1, min.pct = 0.25)
#针对两两cluste之间比较
#cluster5.markers <- FindMarkers(mergObj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
mergeObj.markers <- FindAllMarkers(mergeObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(mergeObj.markers, file = paste(prefix, "_markerGene.csv", sep = ""),  row.names = F)
#对Top10的marker gene进行热图和FeaturePlot，VlnPlot进行展示
top10 <- mergeObj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#热图
png(paste(prefix, "_ClusterHeatmap.png", sep = ""), width = 12*ppi, height = 12*ppi, res = ppi)
DoHeatmap(mergeObj, features = top10$gene, size = 5, lines.width = 50)
dev.off()
#Top基因VlnPlot
dir.create("./Top10MarkerGene_VlnPlot")
setwd("./Top10MarkerGene_VlnPlot")
for(g in top10$gene){
    png(paste0(g, "_VlnPlot.png"), width = 6*ppi, height = 4*ppi, res = ppi)
    print(VlnPlot(mergeObj, features = c(g), pt.size = 0.1))
    dev.off()
}
#Top基因featurePlot
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
#各个cluster细胞数目
write.csv(table(mergeObj@active.ident), file = paste(prefix, "_clusterCellNum.csv", sep = ""), row.names = F)
########################################
#运用singleR进行细胞类型注释分析





















