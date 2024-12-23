library(monocle3)

DefaultAssay(combined2) <- "SCT"
gene_metadata <- data.frame(gene_short_name = rownames(combined2@assays$SCT@data))
rownames(gene_metadata) <- rownames(combined2@assays$SCT@data)

cell_metadata <- combined2@meta.data
expression_matrix <- as(combined2@assays$SCT@counts, "sparseMatrix")

cds <- new_cell_data_set(expression_data = expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)

cds <- preprocess_cds(cds, num_dim = 5)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds, reduction_method = "UMAP")

colData(cds)$seurat_clusters <- as.character(combined2$seurat_clusters)

cds <- cluster_cells(cds)

plot_cells(cds)
plot_cells(cds, color_cells_by = "seurat_clusters", label_branch_points = FALSE, label_leaves = FALSE, label_cell_groups = FALSE,label_roots = F)
plot_cells(cds, color_cells_by = "partition", label_branch_points = FALSE, label_leaves = FALSE, label_cell_groups = FALSE)
plot_cells(cds, color_cells_by = "orig.ident", label_branch_points = FALSE, label_leaves = FALSE, label_cell_groups = FALSE,label_roots = F)
plot_cells(cds, color_cells_by = "seurat_clusters") + facet_wrap(~orig.ident)

cds <- learn_graph(cds, learn_graph_control = list(geodesic_distance_ratio = 0.1, minimal_branch_len = 10))

root_cluster <- "Primary FAPs"
root_cell <- colnames(cds)[which(colData(cds)$seurat_clusters == root_cluster)[1]]

cds <- order_cells(cds, root_cells = root_cell)

pt <- pseudotime(cds)

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.3,label_roots = F)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = TRUE) + facet_wrap(~orig.ident)