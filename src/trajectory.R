library(Seurat)
library(monocle3)
library(SingleR)

ppi<-300

prefix='test'
sc.obj <- readRDS("test.rds")


gene_annotation <- as.data.frame(rownames(sc.obj@reductions[["pca"]]@feature.loadings), row.names = rownames(sc.obj@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
cell_metadata <- as.data.frame(sc.obj@assays[["RNA"]]@counts@Dimnames[[2]], row.names = sc.obj@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- sc.obj@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(sc.obj@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds <- new_cell_data_set(expression_matrix,
                        cell_metadata = cell_metadata,
                        gene_metadata = gene_annotation)

recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
cds@clusters@listData[["TSNE"]][["partitions"]] <- recreate.partition

list_cluster <- sc.obj@meta.data$seurat_clusters
names(list_cluster) <- sc.obj@assays[["RNA"]]@data@Dimnames[[2]]

cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
# cds@clusters@listData[["TSNE"]][["clusters"]] <- list_cluster

cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
# cds@clusters@listData[["TSNE"]][["louvain_res"]] <- "NA"

cds@reducedDims@listData[["UMAP"]] <-sc.obj@reductions[["umap"]]@cell.embeddings
# cds@reducedDims@listData[["TSNE"]] <-sc.obj@reductions[["tsne"]]@cell.embeddings

cds@preprocess_aux$gene_loadings <- sc.obj@reductions[["pca"]]@feature.loadings

# singler.obj = CreateSinglerObject(sc.obj@assays$RNA@counts, annot = NULL, prefix, min.genes = 0,
#                             technology = "DropSeq", species = "Human", citation = "",
#                             ref.list = list(), normalize.gene.length = F, variable.genes = "de",
#                             fine.tune = F, do.signatures = T, clusters = sc.obj@active.ident, do.main.types = T,
#                             reduce.file.size = T, numCores = 64)
# singler.obj$meta.data$orig.ident = sc.obj@meta.data$orig.ident

# singler.obj$meta.data$xy = sc.obj@reductions$umap@cell.embeddings
# singler.obj$meta.data$clusters = sc.obj@active.ident

# colData(cds)$cell.type <- as.vector(singler.obj$other)

cds <- learn_graph(cds, use_partition = T)

png(paste0(prefix,'_trajectoryCluster.png'),width = 4*ppi, height = 4*ppi, res = ppi)
clus <- plot_cells(cds,
                    color_cells_by = 'cluster',
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)
print(clus)
dev.off()

png(paste0(prefix,'_trajectoryCellType.png'),width = 4*ppi, height = 4*ppi, res = ppi)
ctype <- plot_cells(cds,
                    color_cells_by = 'cell.type',
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)
print(ctype)
dev.off()

colData(cds)$clusters <- as.vector(cds@clusters$UMAP$clusters)
cell_ids <- which(colData(cds)[, "clusters"] == "0")

closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

root_pr_nodes
cds = order_cells(cds, root_pr_nodes=root_pr_nodes)

png(paste0(prefix,'_pseudotime.png'),width = 5*ppi, height = 4*ppi, res = ppi)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()

# for(cls in levels(cds@clusters$UMAP$clusters)){
#     clsCells <- which(colData(cds)[, "clusters"] == cls)
#     closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#     closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#     root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[clsCells,]))))]
#     cds = order_cells(cds, root_pr_nodes=root_pr_nodes)
#     print(sprintf('%s_pseudotime_cluster%s.png',prefix,cls))
#     png(sprintf('%s_pseudotime_cluster%s.png',prefix,cls),width = 5*ppi, height = 4*ppi, res = ppi)
#     plot_cells(cds,
#             color_cells_by = "pseudotime",
#             label_cell_groups=FALSE,
#             label_leaves=FALSE,
#             label_branch_points=FALSE,
#             graph_label_size=1.5)
#     dev.off()
# }




