#!/Rscript
#################################################################################
#          Author: liang.wan, 2019-03-11                                        #
#          Contact: lwan@dragondropinn.com                                      #
#          Reference: https://satijalab.org/seurat/pbmc3k_tutorial.html         #
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

###################################################################################
#          setup the seurat object                                                #
###################################################################################
mergeObj <- CreateSeuratObject(raw.data = mergeData, min.cells = 5, min.genes = arg_minGene, project = prefix, names.field = 1, names.delim = "_")
mergeObj
###################################################################################
#          QC and selecting cells for further analysis                            #
###################################################################################
#add mito info and QC plot
mito.genes <- grep(pattern = "^MT-", x = rownames(x = mergeObj@data), value = TRUE)
percent.mito <- Matrix::colSums(mergeObj@raw.data[mito.genes, ])/Matrix::colSums(mergeObj@raw.data)
mergeObj <- AddMetaData(object = mergeObj, metadata = percent.mito, col.name = "percent.mito")

png(paste(prefix, "_QCPlot.png", sep = ""),width = 12*ppi, height = 10*ppi, res = ppi)
VlnPlot(object = mergeObj, features.plot = c("nGene", "nUMI"), nCol = 2, point.size.use = 1, x.lab.rot = TRUE)
dev.off()
png(paste(prefix, "_QCmitoPlot.png", sep = ""),width = 6*ppi, height = 10*ppi, res = ppi)
VlnPlot(object = mergeObj, features.plot = "percent.mito", point.size.use = 1, x.lab.rot = TRUE, y.max = 1)
dev.off()

#Gene Linear plot
png(paste(prefix, "_MitoPlot.png", sep = ""), width = 4*ppi, height = 4*ppi, res = ppi)
GenePlot(object = mergeObj, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.3)
dev.off()
png(paste(prefix, "_GenePlot.png", sep = ""), width = 4*ppi, height = 4*ppi, res = ppi)
GenePlot(object = mergeObj, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.3)
dev.off()

#filter cells
mergeObj <- FilterCells(object = mergeObj, subset.names = c("percent.mito", "nGene"), low.thresholds = c(-Inf, arg_nGene), high.thresholds = c(arg_mito, 2500))
mergeObj

####################################################################################
#           normalization data and find variable genes                             #
####################################################################################
mergeObj <- NormalizeData(object = mergeObj, normalization.method = "LogNormalize", scale.factor = 10000)
png(paste(prefix, "_VariableGenes.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
mergeObj <- FindVariableGenes(object = mergeObj, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 0.5)
dev.off()
cat(paste("The total variable gene number is:", length(x = mergeObj@var.genes), sep = " "))

####################################################################################
#          Scaling the data and removing unwanted sources of variation             #
####################################################################################
mergeObj <- ScaleData(object = mergeObj, vars.to.regress = c("nUMI", "percent.mito"))

####################################################################################
#           Perform linear dimensional reduction                                   #
####################################################################################
mergeObj <- RunPCA(object = mergeObj, pc.genes = mergeObj@var.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 10)
png(paste(prefix, "_PCA.png", sep = ""), width = 6*ppi, height = 9*ppi, res = ppi)
VizPCA(object = mergeObj, pcs.use = 1:9)
dev.off()

mergeObj <- ProjectPCA(object = mergeObj, do.print = FALSE)
#PC Heatmap plot
png(paste(prefix, "_PCHeatmap.png", sep = ""), width = 6*ppi, height = 9*ppi, res = ppi)
PCHeatmap(object = mergeObj, pc.use = 1:9, cells.use = 500, do.balanced = TRUE, label.columns = FALSE,  use.full = FALSE)
dev.off()

####################################################################################
#           Determine statistically significant principal components               #
####################################################################################
#Jack Straw plot
png(paste(prefix, "_JackStraw.png", sep = ""), width = 6*ppi, height = 9*ppi, res = ppi)
mergeObj <- JackStraw(object = mergeObj, num.replicate = 100, display.progress = FALSE)
JackStrawPlot(object = mergeObj, PCs = 1:12)
dev.off()
#PCElbowPlot
png(paste(prefix, "_PCElbow.png", sep = ""), width = 6*ppi, height = 4*ppi, res = ppi)
PCElbowPlot(object = mergeObj)
dev.off()

####################################################################################
#            Cluster the cells                                                     #
####################################################################################
#Find cluster
mergeObj <- FindClusters(object = mergeObj, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = mergeObj)

####################################################################################
#            Run Non-linear dimensional reduction (tSNE)                           #
####################################################################################
#run TSNE
mergeObj <- RunTSNE(object = mergeObj, dims.use = 1:20, do.fast = TRUE, check_duplicates = arg_check)
png(paste(prefix, "_TSNEPlot.png", sep = ""), width = 12*ppi, height = 4*ppi, res = ppi)
p1 <- TSNEPlot(mergeObj, do.return = T, pt.size = 1, group.by = "orig.ident")
p2 <- TSNEPlot(mergeObj, do.label = T, do.return = T, pt.size = 1, label.size = 4)
plot_grid(p1, p2)
dev.off()
#save RDS
saveRDS(mergeObj, file = paste(workDir, "/", prefix, ".rds", sep = ""))

####################################################################################
#            Finding differentially expressed genes (cluster biomarkers)           #
####################################################################################

#find all marker Gene
mergeObj.markers <- FindAllMarkers(object = mergeObj, min.pct = 0.15, only.pos = TRUE)
write.csv(mergeObj.markers, file = paste(prefix, "_markerGene.csv", sep = ""),  row.names = F)
#find Top10 marker Gene
top10 <- mergeObj.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
write.csv(top10, file = paste(prefix, "_top10_markerGene.csv", sep = ""),  row.names = F)
#heatmapPlot of top10 marker Gene
png(paste(prefix, "_ClusterHeatmap.png", sep = ""), width = 8*ppi, height = 4*ppi, res = ppi)
DoHeatmap(object = mergeObj, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE, cex.row = 3, cex.col = 3)
dev.off()

####################################################################################
#           Write CSV record                                                       #
####################################################################################
#write cluster cell number
write.csv(table(mergeObj@ident), file = paste(prefix, "_clusterCellNum.csv", sep = ""), row.names = F)
#write meta Data
write.csv(mergeObj@meta.data[,1:3], file = paste(prefix, "_metaData.csv", sep = ""),  row.names = T)

cat("Job is finished!")
sink()


















































































