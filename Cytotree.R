####################################################
# Showing an example of analysing a single FCS file, a similar approach was used for all sample
# CytoTree analysis for sample 23 #
# No. Of cell counts: 50699 #

####################################################

## Loading packages ## 
library(ggplot2)
library(CytoTree)
library(flowCore)
library(stringr)
library(colorRamps)

# Set a working directory
setwd("~/Suzanne/Cytotree/NSCLC_samples/")

## 1. Import of FCS files ## 
# get the path of the file
fcs.path = "~/Suzanne/Cytotree/NSCLC_samples/data/early_cohort/Preprocessed_FCSfiles/"
# get the FCS file
fcs_file = paste0(fcs.path, "CLEAN_Mix 4_Tumor_NSCLC_with_comp_median_023.fcs")

# Extract the expression matrix from the FCS file
fcs = runExprsExtract(fcs_file, comp = FALSE, transformMethod = "none")
#refine the rownames for easaier use
rownames(fcs) = gsub("CLEAN_Mix 4_Tumor_NSCLC_with_comp_median_", "", rownames(fcs), fixed=TRUE)

# Refine colnames of fcs data
# Note: delete or omit the columns that will not be used 
recol = c(`FSC-A<NA>` = "FSC-A",
          `SSC-A<NA>` = "SSC-A",
          `FJComp-APC-A<NA>` = "CD68",
          `FJComp-Alexa Fluor 700-A<NA>` = "CD33",
          `FJComp-Brilliant Violet 421-A<NA>` = "CD11b", 
          `FJComp-Brilliant Violet 510-A<NA>` = "CD3", 
          `FJComp-Brilliant Violet 605-A<NA>` = "CD15",
          `FJComp-FITC-A<NA>` = "HLA-DR",
          `FJComp-PE-A<NA>` = "PDL-1",
          `FJComp-PE-Cy7-A<NA>` = "CD14", 
          `FJComp-PerCP-Cy5-5-A<NA>` = "CD11c", 
          `FJComp-Qdot 585-A<NA>` = "CD1a",
          `FJComp-Qdot 625-A<NA>` = "CD16", 
          `FJComp-Qdot 705-A<NA>` = "CD19/CD20")
colnames(fcs)[match(names(recol), colnames(fcs))] = recol
fcs = fcs[, recol]

## 2. Create metadata ## 
# Patient_id 
donor_list = "023"
# creating a metadata table with cells, markers, and patient ID
meta.data = data.frame(cell = rownames(fcs), #cells
                       stage = substr(rownames(fcs),1,3)) #patient ID
meta.data$stage = factor(as.character(meta.data$stage), levels = donor_list) 
markers = colnames(fcs) #markers used 

## 3. Clustering and Dimensionlaity reduction ##
# Build CYT object 
cyt = createCYT(raw.data = fcs,
                   markers = markers,
                   meta.data = meta.data,
                   normalization.method = "none",
                   verbose = TRUE)


# Cluster cells by SOM algorithm 
set.seed(12345)
cyt = runCluster(cyt, cluster.method = "som", xdim = 5, ydim = 5) #25 clusters

# Processing Clusters
cyt = processingCluster(cyt)

# Dimensionality reduction
# Run Principal Component Analysis (PCA)
cyt = runFastPCA(cyt)

# Run Uniform Manifold Approximation and Projection (UMAP)
set.seed(12345)
cyt = runUMAP(cyt)


## 4. Visualization ## 
# Visulaize the distribution of cells based on their marker expression in a 2D UMAP plot

# Plot 2D UMAP, and cells are colored by cluster ID
plot2D(cyt_23, item.use = c("UMAP_1", "UMAP_2"), color.by = "cluster.id",show.cluser.id = TRUE,
       alpha = 0.5, main = "NSCLC:Clusters", category = "categorical", plot.theme = theme_classic()) + theme(aspect.ratio = 1)

# Plot 2D UMAP, and cells are colored by patient donor
plot2D(cyt_23, item.use = c("UMAP_1", "UMAP_2"), color.by = "stage",
       alpha = 1, main = "NSCLC:Stage/samples", category = "categorical",plot.theme = theme_classic()) + theme(aspect.ratio = 1)

# defining colours for 2D UMAP plot
col_fun = colorRamp2(c(-3, -1.5, 0, 1.5, 3), c("#053061","#579EC9","#FFFFFF","#D96651","#67001F"))

#extracting the expression matrix from CYT object for plotting
cyt_exp = cyt@raw.data 

# Plot 2D UMAP, and cells are colored according to CD68 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD68",
       alpha = 1, main = "NSCLC:CD68", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD68"]),max(cyt_exp[,"CD68"]),by=0.01)))
                         
# Plot 2D UMAP, and cells are colored according to CD33 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD33",
       alpha = 1, main = "NSCLC:CD33", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD33"]),max(cyt_exp[,"CD33"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD11b expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD11b",
       alpha = 1, main = "NSCLC:CD11b", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD11b"]),max(cyt_exp[,"CD11b"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD3 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD3",
       alpha = 1, main = "NSCLC:CD3", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD3"]),max(cyt_exp[,"CD3"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD15 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD15",
       alpha = 1, main = "NSCLC:CD15", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD15"]),max(cyt_exp[,"CD15"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to HLA-DR expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "HLA-DR",
       alpha = 1, main = "NSCLC:HLA-DR", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"HLA-DR"]),max(cyt_exp[,"HLA-DR"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to PD-L1 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "PDL-1",
       alpha = 1, main = "NSCLC:PDL-1", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"PDL-1"]),max(cyt_exp[,"PDL-1"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD14 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD14",
       alpha = 1, main = "NSCLC:CD14", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD14"]),max(cyt_exp[,"CD14"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD11c expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD11c",
       alpha = 1, main = "NSCLC:CD11c", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD11c"]),max(cyt_exp[,"CD11c"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD1a expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD1a",
       alpha = 1, main = "NSCLC:CD1a", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD1a"]),max(cyt_exp[,"CD1a"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD16 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD16",
       alpha = 1, main = "NSCLC:CD16", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD16"]),max(cyt_exp[,"CD16"]),by=0.01)))

# Plot 2D UMAP, and cells are colored according to CD19/CD20 expression 
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD19/CD20",
       alpha = 1, main = "NSCLC:CD19/CD20", category = "numeric") +
  scale_colour_gradientn(colors = col_fun(seq(min(cyt_exp[,"CD19/CD20"]),max(cyt_exp[,"CD19/CD20"]),by=0.01)))

# Plot heatmap of clusters
plotClusterHeatmap(cyt)

## 5.Building trajectory ## 
# Build tree using expression matrix and UMAP
cyt = buildTree(cyt, dim.type = "umap", dim.use = 1:2)
cyt@meta.data$branch.id = paste0("B", cyt_23@meta.data$branch.id)

# Fetch plot meta information for each cluster
cluster_cyt = fetchClustMeta(cyt)
knitr::kable(head(cluster_cyt))

## 6. Visualization ## 
# Tree plot using branch.id
plotTree(cyt, color.by = "branch.id", show.node.name = T, cex.size = 1)

# Plot heatmap using marker expression on branches
plotBranchHeatmap(cyt, colorRampPalette(c("#00599F", "#FFFFFF", "#FF3222"))(100),
                  clustering_method = "complete")

# Expression of CD68 overlaid on tree plot
plotTree(cyt, color.by = "CD68", show.node.name = T, cex.size = 1) +
  ggtitle("CD68")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD68"]),max(cluster_cyt[,"CD68"]),by=0.01)))

# Expression of CD33 overlaid on tree plot
plotTree(cyt, color.by = "CD33", show.node.name = T, cex.size = 1) +
  ggtitle("CD33")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD33"]),max(cluster_cyt[,"CD33"]),by=0.01)))

# Expression of CD11b overlaid on tree plot
plotTree(cyt, color.by = "CD11b", show.node.name = T, cex.size = 1) +
  ggtitle("CD11b")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD11b"]),max(cluster_cyt[,"CD11b"]),by=0.01)))

# Expression of CD3 overlaid on tree plot
plotTree(cyt, color.by = "CD3", show.node.name = T, cex.size = 1) +
  ggtitle("CD3")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD3"]),max(cluster_cyt[,"CD3"]),by=0.01)))

# Expression of CD15 overlaid on tree plot
plotTree(cyt, color.by = "CD15", show.node.name = T, cex.size = 1) +
  ggtitle("CD15")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD15"]),max(cluster_cyt[,"CD15"]),by=0.01)))

# Expression of HLA-DR overlaid on tree plot
plotTree(cyt, color.by = "HLA-DR", show.node.name = T, cex.size = 1) +
  ggtitle("HLA-DR")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"HLA-DR"]),max(cluster_cyt[,"HLA-DR"]),by=0.01)))

# Expression of PD-L1 overlaid on tree plot
plotTree(cyt, color.by = "PDL-1", show.node.name = T, cex.size = 1) +
  ggtitle("PDL-1")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"PDL-1"]),max(cluster_cyt[,"PDL-1"]),by=0.01)))

# Expression of CD14 overlaid on tree plot
plotTree(cyt, color.by = "CD14", show.node.name = T, cex.size = 1) +
  ggtitle("CD14")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD14"]),max(cluster_cyt[,"CD14"]),by=0.01)))

# Expression of CD11c overlaid on tree plot
plotTree(cyt, color.by = "CD11c", show.node.name = T, cex.size = 1) +
  ggtitle("CD11c")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD11c"]),max(cluster_cyt[,"CD11c"]),by=0.01)))

# Expression of CD1a overlaid on tree plot
plotTree(cyt, color.by = "CD1a", show.node.name = T, cex.size = 1) +
  ggtitle("CD1a")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD1a"]),max(cluster_cyt[,"CD1a"]),by=0.01)))

# Expression of CD16 overlaid on tree plot
plotTree(cyt, color.by = "CD16", show.node.name = T, cex.size = 1) +
  ggtitle("CD16")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD16"]),max(cluster_cyt[,"CD16"]),by=0.01)))

# Expression of CD19/CD20 overlaid on tree plot
plotTree(cyt, color.by = "CD19/CD20", show.node.name = T, cex.size = 1) +
  ggtitle("CD19/CD20")+
  scale_colour_gradientn(colors = col_fun(seq(min(cluster_cyt[,"CD19/CD20"]),max(cluster_cyt[,"CD19/CD20"]),by=0.01)))

## 7. Optimization ##
# Run differential expression analysis for markers in clusters (source code modified to find the differentially expressed markers in clusters.,rather than using branch.id)
# For the modified code refer "NSCLC_Project/runDiff_source_code.R"
diff.cyt = runDiff(cyt)

# Branches are assigned to the immune cells population 
# Re-naming the branch.id based based on the markers expression
branch.id[branch.id %in% c("B2","B5")] = "T cells"
branch.id[branch.id %in% c("B4")] = "Myeloid cells"
branch.id[branch.id %in% c("B1")] = "B cells"
branch.id[branch.id %in% c("B3")] = "Undefined"

# From the biological analysis, found some clusters are assigned to the wrong branch/immune population; so after checking clusters
# were re-assigned to their respective immune cell population 
branch.id[cyt@meta.data$cluster.id %in% c(17,4,5,2,13,8,10,19)] = "T cells"
branch.id[cyt@meta.data$cluster.id %in% c(9,16,23,6,20,21,25,1,14)] = "B cells"
branch.id[cyt@meta.data$cluster.id %in% c(18)] = "Neutrophils"
branch.id[cyt@meta.data$cluster.id %in% c(12,24)] = "Monocytes"
branch.id[cyt@meta.data$cluster.id %in% c(3,15)] = "DCs"
branch.id[cyt@meta.data$cluster.id %in% c(11)] = "NKT cells"
branch.id[cyt@meta.data$cluster.id %in% c(7)] = "NK cells"
branch.id[cyt@meta.data$cluster.id %in% c(22)] = "Undefined"

# Refine branch.id
cyt@meta.data$branch.id = branch.id

# Plot tree plot after re-assigning the clusters
plotTree(cyt, color.by = "branch.id", show.node.name = TRUE, cex.size = 1)

# Plot 2D UMAP plot for branches
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "branch.id",
       alpha = 1,main = "NSCLC sample 23", category = "categorical", show.cluser.id = F,
       plot.theme = theme_classic(),show.cluser.id.size = 3)



