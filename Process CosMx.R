suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(SeuratDisk)
  library(SingleR)
  library(celldex)
  library(ggrastr)
  library(ggridges)
  library(RColorBrewer)
  library(limma)
  library(edgeR)
  library(Glimma)
  library(sccomp)
  library(ComplexHeatmap)
  library(cluster)
  library(factoextra)
  library(gplots)
  library(corrplot)
  library(Polychrome)
  library(parallel)
  library(GSA)
  library(scDblFinder)
  library(BiocParallel)
  library(harmony)
  library(cowplot)
  library(pdftools)
  library(ggrepel)
  library(progressr)
  library(gridExtra)
})

# Ideas
# Try mapping spatial onto MMR UMAP with sctransform counts

# Set a homedir for the analysis so as to no use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R") 

# Set output folders
plots <- "./Results_paper/Plots/"

md_old <- read_csv("./raw/md_liam.csv")%>%
  mutate(fov = as.numeric(fov))%>%
  group_by(Slide)%>%
  mutate(fov = (fov+1)- min(fov))%>%
  ungroup()%>%
  mutate(Slide = gsub(" ", "_", Slide))%>%
  mutate(Barcode = paste0(Slide, "_", cell_ID, "_", fov))

tb <- md_old%>%
  filter(Slide == "Tumour_B")

ta <- md_old%>%
  filter(Slide == "Tumour_A")

# Read in the metadata
fov_line_up <- read_csv("./raw/Run5850_2186_2/Sample_to_FOV.csv")

unique(fov_line_up$Slide)

# Read in the new metadata
md_new_1 <- read_csv("./raw/Run5850_2186_1/Run5850_2186_1_metadata_file.csv")%>%
  mutate(Slide = "Run5850.2186.1")

md_new_2 <- read_csv("./raw/Run5850_2186_2/Run5850_2186_2_metadata_file.csv")%>%
  mutate(Slide = "Run5850.2186.2")

md_new_all <- bind_rows(md_new_1, md_new_2)%>%
  mutate(fov = as.numeric(fov))%>%
  left_join(fov_line_up)%>%
  group_by(Slide)%>%
  mutate(fov = (fov+1)- min(fov))%>%
  ungroup()%>%
  mutate(Slide = gsub(" ", "_", Slide))%>%
  mutate(Barcode = paste0(Slide, "_", cell_ID, "_", fov))

# Combine all the metadata
md_combined <- bind_rows(md_old, md_new_all)

# Read the counts only out of the nanostring objects
Tumor_A <- get_nano_counts(data.dir = "./raw/Run5654_Tumor_A/",sample_name = "Tumour_A")
Tumor_B <- get_nano_counts(data.dir = "./raw/Run5629_Tumor_B/",sample_name = "Tumour_B")
Normal_A <- get_nano_counts(data.dir = "./raw/Run5654_Normal_A/",sample_name = "Normal_A")
new_1 <- get_nano_counts(data.dir = "./raw/Run5850_2186_1/",sample_name = "Run5850.2186.1")
new_2 <- get_nano_counts(data.dir = "./raw/Run5850_2186_2/",sample_name = "Run5850.2186.2")

# Merge in tumour counts
merged <-  merge(x = Tumor_A, y = c(Tumor_B, Normal_A, new_1, new_2))%>%
  JoinLayers()

qsave(merged, "./Intermediate/Merged_raw_unfiltered.qs")
#merged <- qread("./Intermediate/Merged_raw_unfiltered.qs")

# Make a plot of the top expressed genes
gene_exprs <- rowMeans(merged@assays$RNA$counts)

# Make a plot that compares single probes to single cell gene expression
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

# Get the genes from the MMR atlas
gene_exprs_mmr <- rowMeans(mmr_atlas@assays$RNA@counts)

gene_exprs_df_mmr <- data.frame(gene_exprs_mmr)%>%
  rownames_to_column("Gene")

gene_exprs_df <- data.frame(gene_exprs)%>%
  rownames_to_column("Gene")%>%
  arrange(gene_exprs)%>%
  mutate(nums = 1:n())%>%
  left_join(gene_exprs_df_mmr)%>%
  mutate(labs = replace(Gene, nums < 960, NA))

# Plot just the most expressed genes
gene_exprs_df_top <- gene_exprs_df%>%
  filter(nums > 950)

# Plot the most highly expressed probes
figure_s_0_a <- ggplot(data = gene_exprs_df_top, aes(x = nums, y= gene_exprs, label = labs))+
  geom_point()+
  geom_text_repel(max.overlaps = 20,force = 10, colour = "blue", fontface = "italic", size =2)+
  labs(x = "Probe index", y = "Mean probes per cell")+
  blank_theme

figure_s_0_a

ggsave(plot = figure_s_0_a,"./Results_paper/Plots/CosMx panel v1 top expressed genes per cell.pdf", 
       width = 80, height = 50, units = "mm")

cor <- cor.test(gene_exprs_df$gene_exprs, gene_exprs_df$gene_exprs_mmr)

figure_s_0_b <- ggplot(data = gene_exprs_df, aes(x = log2(gene_exprs), y= log2(gene_exprs_mmr), label = labs))+
  geom_point(size =0.3)+
  annotate("text", x = 3, y =-2, label = paste0("Pearson's r = ", round(cor$estimate,2)), size = 2.5)+
  labs(x =  expression('Log'[2]*' mean CosMx probe count'), y = expression('Log'[2]*' mean 10x read count'))+
  blank_theme+
  theme(aspect.ratio = 1)

figure_s_0_b

# Plot the expression of the negative probes
negprb_plt <- merged@assays$RNA$counts[grepl("NegPrb", rownames(merged)),]
negprb_plt[1:5,1:5]

negprb_plt <- CreateSeuratObject(negprb_plt)

negprb_plt$Barcode <- colnames(negprb_plt)
# Join the datasets on 'Barcode'
negprb_plt@meta.data <- left_join(negprb_plt@meta.data, md_combined)
# Reset the rownames
rownames(negprb_plt@meta.data) <- negprb_plt@meta.data$Barcode

slide_short <- data.frame(Slide= negprb_plt$Slide)%>%
  left_join(slide_alias)

negprb_plt$Slide_short <- slide_short$slide_names
negprb_plt$Sample_Slide <- paste0(negprb_plt$Sample, " ", negprb_plt$Slide_short)

probes <- rownames(negprb_plt)

# Make a dot plot of the negative probes
qc1_a <- DotPlot(object = negprb_plt, features = probes, group.by = "Sample_Slide")+
  blank_theme+
  labs(y = "Sample + slide", x = "Probes")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(4, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))

qc1_a

# Drop the negative probes and lncRNAs that don't help cell type ID
negprb <- rownames(merged)[grepl("NegPrb", rownames(merged))]
merged <- merged@assays$RNA$counts
# https://kb.10xgenomics.com/hc/en-us/articles/360004729092-Why-do-I-see-high-levels-of-Malat1-in-my-gene-expression-data-
# 10x says dead and dying cells might have high MALAT1
# These genes do not help with celltype ID
merged <- merged[!rownames(merged) %in% c(negprb, "NEAT1", "MALAT1"),]

# Make a seurat object, only keeping genes detected in > 50 cells
merged <- CreateSeuratObject(merged, min.cells = 50)
colnames(merged)[1:5]
merged$Barcode <- colnames(merged)
unique(merged$orig.ident)
# Join the datasets on 'Barcode'
merged@meta.data <- left_join(merged@meta.data, md_combined)
# Reset the rownames
rownames(merged@meta.data) <- merged@meta.data$Barcode

# Save the object containing all cells
qsave(merged, "./Intermediate/All_cells.qs")

# merged <- qread("./Intermediate/All_cells.qs")

# Make a higher quality object to 
# Get the total number of cells
total_start <- length(merged$Barcode)

# Read in the results of per-sample clustering/annotation
# Keep only HQ samples/slides from silhouette analysis 

# Loop over each sample and calculate a silhouette score
files <- list.files("./Results_paper/analyses/single_sample_analyses/", pattern = "*.qs", recursive = T, full.names = T)

# A list for sil QC stats
sil_list <- list()
# Make a list of the cells that passed QC
good_cells_list <- list()
plotlist <- list()
for (i in 1:length(files)){
  
  print(i)
  
  file <- files[i]
  
  sample_name <- gsub("_seurat_processed.qs", "", basename(file))
  
  sample <- qread(file)
  
  dim(sample)
  
  sil_df <- data.frame(Sample = sample_name, median_sil = median(sample$sil), 
                       medcount = median(sample$nCount_RNA),
                       clusters = max(as.numeric(sample$seurat_clusters)))
  
  sil_score <- round(sil_df$median_sil,2)
  
  sample_short <- gsub(".*_", "",sample_name)
  slide <- gsub(paste0("_", sample_short), "",sample_name)
  slide <- slide_alias$slide_names[slide_alias$Slide == slide]
  
  plt <- DimPlot(sample, reduction = "umap", label = T, group.by = "Manual_toplevel_pred", 
                 cols = cell_type_colors, raster = T, label.size = 1)+
    ggtitle(paste0(sample_short, " slide ", slide, ", med sil = ", sil_score))+
    blank_theme+
    NoLegend()+
    labs(x = NULL, y = NULL)+
    theme(axis.text.x = element_text(size = 4),
                  axis.title.x = element_text(size = 5),
                  axis.title.y = element_text(size = 5),
                  axis.text.y = element_text(size = 4),
          plot.title = element_text(size = 5.5))
  
  plotlist[[i]] <- plt
  
  sil_list[[i]] <- sil_df
  
  good_cells_list[[i]] <- sample@meta.data
  
}

# Read in differentiation scores
differentiation <- read_csv("./raw/Differentiation.csv")%>%
  mutate(Sample = gsub("ORG", "", Identifier))%>%
  dplyr::select(-Identifier)

# Read in all the clinical data
clinical <- read_csv("./raw/Patients for spatial clinical info.csv")%>%
  filter(!is.na(`ORG #`))%>%
  select(Sample = `ORG #`, `Age at Dx`, `Primary tumour location`,
         Sex, Sidedness, `MMR status`, `Surgery date`)%>%
  mutate(Sample = gsub("ORG", "", Sample))%>%
  left_join(differentiation)%>%
  mutate(Donor = gsub("T|N","", Sample))%>%
  write_csv("./Results_paper/Tables/Table S1. Sample information.csv")

length(unique(clinical$Donor))

all_sil <- bind_rows(sil_list)%>%
  mutate(Slide_sample = Sample)%>%
  mutate(Sample = gsub(".*_", "", Slide_sample))%>%
  left_join(clinical)%>%
  mutate(Slide = gsub("_[^_]+$", "", Slide_sample))

sil_plot <- all_sil%>%
  filter(!is.na(`Surgery date`))%>%
  mutate(`Surgery date` = as.Date(`Surgery date`,format = "%d/%m/%Y"))%>%
  left_join(slide_alias)%>%
  mutate(Differentiation = replace(Differentiation, is.na(Differentiation), "NS"))

date_num <- as.numeric(sil_plot$`Surgery date`)
date_cor <- cor.test(date_num, sil_plot$median_sil)

r2_date <- date_cor$estimate^2

# Read in the results of the per-FOV silhouette analysis
per_fov_sil <- read_csv("./Results_paper/Tables/Per FOV silhouette scores.csv")

good_samples <- per_fov_sil%>%
  filter(median_sil > -0.025)%>%
  mutate(Slide_sample = paste0(Slide, "_", Sample))

count_cor <- cor.test(sil_plot$medcount, sil_plot$median_sil)
count_cor

r2_count <- count_cor$estimate^2
r2_count

single_sample_plots <- plot_grid(plotlist = plotlist,ncol = 4)

ggsave(plot = single_sample_plots,"./Results_paper/Plots/Single sample celltypes.png", 
       width = 170, height = 340, units = "mm", dpi = 600)
ggsave(plot = single_sample_plots,"./Results_paper/Plots/Single sample celltypes.pdf", 
       width = 170, height = 340, units = "mm")

sil_plot_comp <- ggplot(data = sil_plot, aes(x = medcount, y = median_sil, label = Sample, colour = slide_names))+
  geom_point()+
  geom_text(nudge_x = 20,size = 2)+
  geom_hline(yintercept = -0.025, linetype = 2,linewidth = 0.5)+
  blank_theme+
  annotate("text",size=3, x = 150, y = 0.1, label = paste0("Pearson's r = ", round(count_cor$estimate,2)))+
  labs(x = "Median probe count", y = "Median silhouette score", colour = "Slide", title = "")+
  guides(colour = "none")

sil_plot_comp

ggsave(plot = sil_plot_comp, "./Results_paper/Extra_plots/Silhouette plot.pdf", width = 6, height = 5)

sil_plot$Differentiation <- factor(sil_plot$Differentiation, levels = c("NS", "poor", "mod", "well"))

# Plot sil score vs differentiation looking at tumours
diff_sil <- ggplot(data = sil_plot, aes(x = Differentiation, y = median_sil))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(aes(colour = Slide_sample), width = 0.3, height = 0)+
  guides(colour = "none")+
  blank_theme+
  labs(y = "Median silhouette score (sample/slide)")

sil_summaries <- plot_grid(diff_sil,sil_plot_comp ,ncol = 2, labels = c("a", "b"), label_size = 8, align = "h")

sil_summaries

ggsave(plot = sil_summaries,"./Results_paper/Plots/Silhouette summaries.png", 
       width = 170, height = 70, units = "mm", dpi = 600)
ggsave(plot = sil_summaries,"./Results_paper/Plots/Silhouette summaries.pdf", 
       width = 170, height = 70, units = "mm")

# Possibly drop FOV 24 from 280T as it may be a normal
# fov <- qread("/oldvol/Data/Spatial/Results_paper/analyses/single_sample_analyses/120131_280T/120131_280T_seurat_processed.qs")
# t4 <- fov[,fov$fov == 24]
# DimPlot(t4, group.by = "fov", label = T)
# DimPlot(fov, group.by = "fov", label = T)

# Combine the metadata from the good cells
good_cell_md <- bind_rows(good_cells_list)%>%
  mutate(Slide_sample = paste0(Slide, "_", Sample))%>%
  filter(Slide_sample %in% good_samples$Slide_sample)

# Drop any cells not in 'good' samples
merged <- merged[,merged$Barcode %in% good_cell_md$Barcode]

# Record after this filter step
total_good_samples <- length(merged$Barcode)

# Drop any cells with low counts 
merged <- merged[, merged$nCount_RNA >=50]

# Record after this filter step
total_nCount_RNA <- length(merged$Barcode)

# Change slide names
slide_short <- data.frame(Slide= merged$Slide)%>%
  left_join(slide_alias)

merged$Slide_short <- slide_short$slide_names
merged$Sample_Slide <- paste0(merged$Sample, " ", merged$Slide_short)

VlnPlot(merged, features = c("nCount_RNA"),
        group.by = "Sample_Slide", pt.size = -1)

# Calculate the MAD values for counts features
# Do this for each individual sample
md <- merged@meta.data%>%
  group_by(Sample)%>%
  mutate(m = median(nFeature_RNA))%>%
  mutate(s = mad(nFeature_RNA))%>%
  mutate(robzscore_nFeature_RNA = abs((nFeature_RNA - m) / (s)))%>%
  mutate(m = median(nCount_RNA))%>%
  mutate(s = mad(nCount_RNA))%>%
  mutate(robzscore_nCount_RNA = abs((nCount_RNA - m) / (s)))%>%
  ungroup()%>%
  data.frame()

# Reset the rownames
rownames(md) <- md$Barcode
merged@meta.data <- md

min_QC_robz <- 3

# Subset down based on robust Z score cutoffs for each sample
merged <- subset(merged, subset = robzscore_nFeature_RNA < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

total_robz_cutoff <- length(merged$Barcode)

# Doublet detection/filtering done at the single sample level

# Use scDblFinder on the data
# The tool by defualt expects a 10x like doublet rate.
# I'm not sure this will be the case for nanostring but there might
# be cells with overlapping cytoplasm, so worth trying to clean up

# Run without cluster info
dbs <- scDblFinder(merged@assays$RNA$counts,
                   samples = merged$Sample_Slide,
                   BPPARAM=MulticoreParam(20))


# Look at the doublet rate
table(dbs$scDblFinder.class)

# Drop doublets
merged$scDblFinder.class <- dbs$scDblFinder.class
merged <- merged[,merged$scDblFinder.class == "singlet"]

# ~10% of cells removed
# 21168/(21168 +189768)

# Track the total after cell count filtering
total_doublet_cutoff <- length(merged$Barcode)

qc_df <- data.frame(`Starting cells` = total_start,
                    `Cells in good samples` = total_good_samples,
                    `Cells passing RobZ` = total_robz_cutoff,
                    check.names = F)%>%
  write_csv("./Results_paper/Tables/QC_filtering.csv")

# Make a plot of QC filtering steps
qc_filt <- qc_df%>%
  gather(Step, Cells)%>%
  mutate(Step = factor(Step, levels = Step))

qc1_b <- ggplot(data  = qc_filt, aes(x = Step, y = Cells))+
  geom_bar(stat = "identity")+
  blank_theme+
  labs(x = "Filtering step", y = "Total cells")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Make a QC figure for the supp 
qc1_c <- VlnPlot(merged, features = c("nCount_RNA"),
                 group.by = "Sample_Slide", pt.size = -1)+
  blank_theme+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample + slide", y = "Total probes", title = NULL)

# Set narrower line widths
qc1_c$layers[[1]]$aes_params$size = 0.25

qc1_d <- VlnPlot(merged, features = c("nFeature_RNA"),
                 group.by = "Sample_Slide", pt.size = -1)+
  blank_theme+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Sample + slide", y = "Unique probes (genes detected)",title = NULL)

qc1_d$layers[[1]]$aes_params$size = 0.25

# Log normalise
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst")

# Identify the 20 most highly variable genes
t20 <- head(VariableFeatures(merged), 20) 

# Show the most varible features
plot1 <- VariableFeaturePlot(merged)
qc1_e <- LabelPoints(plot = plot1, points = t20, repel = TRUE, size = 2, fontface = "italic")+
  blank_theme+
  theme(legend.position = "top")+
  NoLegend()

top_row <- plot_grid(qc1_a, qc1_b, labels = c("a", "b"), label_size = 8, rel_widths = c(2.3,1), ncol = 2)
bottom_row <- plot_grid(qc1_c,qc1_d,qc1_e, labels = c("c", "d", "e"), label_size = 8, align = "h", axis = "bt", nrow = 1)

Figure_s1 <- plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1.2,1))
ggsave(plot = Figure_s1,"./Results_paper/Plots/Figure S1.pdf", width = 170, height = 170, units = "mm")

# Scale data and run PCA
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))

# A lot of B cell genes in the variable features
print(merged[["pca"]], dims = 1:5, nfeatures = 5)

# 30 dims looks like plenty
qc2_a <- ElbowPlot(merged, ndims = 30)+
  blank_theme

# Find neighbour cells (in PCA space, not real space)
merged <- FindNeighbors(merged, dims = 1:30)

# Run with default cluster params
merged <- FindClusters(merged)

# Run UMAP
merged <- RunUMAP(merged, dims = 1:30)

DimPlot(merged, reduction = "umap", label = T, group.by = "seurat_clusters")

qc2_b <- DimPlot(merged, reduction = "umap", group.by = "Slide_short")+
  blank_theme+
  labs(colour = "Slide", title = NULL, x = "UMAP\n1", y = "UMAP\n2")

FeaturePlot(merged, features = "nCount_RNA")

# Run harmony to remove treatment effect and donor effect
merged <- RunHarmony(merged, c("Sample", "Slide"), reduction="pca", reduction.save="harmony")
merged <- FindNeighbors(merged, reduction="harmony", dims=1:30)
merged <- FindClusters(merged, cluster.name = "harmony_clusters")
merged <- RunUMAP(merged, reduction="harmony", dims=1:30, reduction.name="umap_harmony")

FeaturePlot(merged, features = "nCount_RNA", max.cutoff = "q95", reduction = "umap_harmony")

# Add on T/N and Donor annotation 
merged$Sample_type <- ifelse(grepl("T", merged$Sample), "T", "N")
merged$Donor <- gsub("T|N","", merged$Sample)

qc2_c <- DimPlot(merged, reduction="umap_harmony", group.by="harmony_clusters", label = T, label.size = 2)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc2_c

qc2_d <- DimPlot(merged, reduction="umap_harmony", group.by="Slide_short")+
  blank_theme+
  theme()+
  labs(colour = "Slide", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc2_d

qc2_e <- DimPlot(merged, reduction="umap_harmony", group.by="Sample", label = F, label.size = 2)+
  blank_theme+
  theme(
    legend.key.size = unit(4, 'mm'))+
  labs(colour = "Donor", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc2_e

qc2_f <- DimPlot(merged, reduction="umap_harmony", group.by="Sample_type", label = F, label.size = 2)+
  blank_theme+
  labs(colour = "Sample type", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc2_f

top_row <- plot_grid(qc2_a, qc2_b, labels = c("a", "b"), label_size = 8, rel_widths = c(1,1.2), ncol = 2,align = "h", axis = "bt")
middle_row <- plot_grid( qc2_c, qc2_d, labels = c("c", "d"), label_size = 8, rel_widths = c(1,1.2), align = "h", axis = "bt", nrow = 1)
bottom_row <- plot_grid( qc2_e,qc2_f, labels = c("e", "f"), label_size = 8, rel_widths = c(1,1), align = "h", axis = "bt", nrow = 1)

# Save the figure
Figure_s2 <- plot_grid(top_row, middle_row, bottom_row, nrow = 3, rel_heights = c(1,1,1),align = "hv", axis = "bt")
ggsave(plot = Figure_s2,"./Results_paper/Plots/Figure S2.pdf", width = 170, height = 180, units = "mm")

# Start doing some cell type annotations

# Start with the human primary cell atlas
ref.data <- HumanPrimaryCellAtlasData()

# Run singleR and compare against the human primary cell atlas
predictions <- SingleR(test=merged@assays$RNA$counts,
                       ref=ref.data, labels=ref.data$label.main, num.threads = 30)
merged$PrimaryCellAtlas_pred_main <- predictions$labels

# Get the fine labels from the primary cell atlas
predictions <- SingleR(test=merged@assays$RNA$counts,
                       ref=ref.data, labels=ref.data$label.fine, num.threads = 30)
merged$PrimaryCellAtlas_pred_fine <- predictions$labels

# Add on the annotations for the MMR atlas
# Run singleR and compare against the MMR reference
# Aggreate the reference
predictions <- SingleR(test=merged@assays$RNA$data, 
                       ref=mmr_atlas@assays$RNA@data, labels=mmr_atlas$clMidwayPr,
                       num.threads = 30, aggr.ref = T)

merged$MMR_pred_midlevel <- predictions$labels

# Run the MMR atlas at the toplevel
predictions <- SingleR(test=merged@assays$RNA$data, 
                       ref=mmr_atlas@assays$RNA@data, labels=mmr_atlas$clTopLevel,
                       num.threads = 30, aggr.ref = T)

merged$MMR_pred_toplevel <- predictions$labels

# Run the MMR atlas at the toplevel
predictions <- SingleR(test=merged@assays$RNA$data, 
                       ref=mmr_atlas@assays$RNA@data, labels=mmr_atlas$cl295v11SubFull,
                       num.threads = 30, aggr.ref = T)

merged$MMR_pred_lowlevel <- predictions$labels

qc_3_a <- DimPlot(merged, reduction="umap_harmony", group.by="PrimaryCellAtlas_pred_main", label = T, label.size = 2, repel = T)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "Main primary cell atllas pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc_3_a

qc_3_b <- DimPlot(merged, reduction="umap_harmony", group.by="MMR_pred_toplevel", label = T, label.size = 2)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "MMR atlas toplevel pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc_3_b

qc_3_c <- DimPlot(merged, reduction="umap_harmony", group.by="MMR_pred_midlevel", label = T, label.size = 2, cols = cell_type_colors)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "MMR atlas midlevel pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc_3_c

# I don't think resoltion is there for the lowlevel preds
qc_3_d <- DimPlot(merged, reduction="umap_harmony", group.by="MMR_pred_lowlevel", label = T, 
                  label.size = 2, repel = T)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "MMR atlas lowlevel pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc_3_d

# Save the figure
Figure_s3 <- plot_grid(qc_3_a, qc_3_b, qc_3_c, qc_3_d, labels = c("a", "b", "c", "d"), label_size = 8, align = "h", axis = "bt", nrow = 2)
ggsave(plot = Figure_s3,"./Results_paper/Plots/Figure S3.pdf", width = 170, height = 170, units = "mm")

# Do celltype annotation
# Start with the MMR atlas midlevel pred
merged$Manual_toplevel_pred <- merged$MMR_pred_midlevel

# Get the celltype counts
ct_counts <- table(merged$MMR_pred_midlevel)%>%
  data.frame()%>%
  arrange(-Freq)

# Take a look at some clusters in fine deatil
# Try with the fibroblasts
fibro <- merged[, merged$MMR_pred_midlevel %in% c("Endo", "Fibro")]

# Take a look at the object
fibro

# Run UMAP clustering on the dataset
fibro <- RunPCA(fibro, verbose = FALSE)
ElbowPlot(fibro, ndims = 30)+
  blank_theme
fibro <- RunHarmony(fibro, c("Sample", "Slide"), reduction="pca", reduction.save="harmony")
fibro <- FindNeighbors(fibro, reduction="harmony", dims=1:15)
fibro <- FindClusters(fibro, cluster.name = "harmony_clusters")
fibro <- RunUMAP(fibro, reduction="harmony", dims=1:15, reduction.name="umap_harmony")
# Plot the results
DimPlot(fibro,label = T,reduction = "umap_harmony")
c5 <- FindMarkers(fibro, ident.1 = 5)
c1 <- FindMarkers(fibro, ident.1 = 1)
DimPlot(fibro,label = T, group.by = "MMR_pred_midlevel",reduction = "umap_harmony")
DimPlot(fibro,label = T, group.by = "PrimaryCellAtlas_pred_main",reduction = "umap_harmony")+
  NoLegend()
DimPlot(fibro,label = T, group.by = "Sample",reduction = "umap_harmony")+
  NoLegend()
FeaturePlot(fibro, features = c("ACTA2") ,reduction = "umap_harmony")
FeaturePlot(fibro, features = c("FN1") ,reduction = "umap_harmony")
FeaturePlot(fibro, features = c("EPCAM") ,reduction = "umap_harmony")
FeaturePlot(fibro, features = c("ANGPT2") ,reduction = "umap_harmony")

# Try with mac/mono
mac_neu <- merged[, merged$MMR_pred_midlevel %in% c("Mono", "Macro", "Granulo", "DC")]
unique(merged$PrimaryCellAtlas_pred_main)

# Take a look at the object
mac_neu
# Run UMAP clustering on the dataset
mac_neu <- RunPCA(mac_neu, verbose = FALSE)
ElbowPlot(mac_neu, ndims = 30)+
  blank_theme
mac_neu <- RunHarmony(mac_neu, c("Sample", "Slide"), reduction="pca", reduction.save="harmony")
mac_neu <- FindNeighbors(mac_neu, reduction="harmony", dims=1:15)
mac_neu <- FindClusters(mac_neu, cluster.name = "harmony_clusters")
mac_neu <- RunUMAP(mac_neu, reduction="harmony", dims=1:15, reduction.name="umap_harmony")
# Plot the results
DimPlot(mac_neu,label = T,reduction = "umap_harmony")
#c5 <- FindMarkers(mac_neu, ident.1 = 5, ident.2 = 1)
DimPlot(mac_neu,label = T, group.by = "MMR_pred_midlevel",reduction = "umap_harmony")
DimPlot(mac_neu,label = T, group.by = "Sample",reduction = "umap_harmony")
DimPlot(mac_neu,label = T, group.by = "Sample_type",reduction = "umap_harmony")
DimPlot(mac_neu,label = T, group.by = "Donor",reduction = "umap_harmony")
DimPlot(mac_neu,label = T, group.by = "PrimaryCellAtlas_pred_main",reduction = "umap_harmony")+
  NoLegend()
FeaturePlot(mac_neu, features = c("S100A8") ,reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("S100A9") ,reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("LYZ") ,reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("HCAR2"),reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("CD3E"),reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("MS4A4A"),reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("IGHG1"),reduction = "umap_harmony")
FeaturePlot(mac_neu, features = c("CXCL8"),reduction = "umap_harmony")
FeaturePlot(merged, features = c("MKI67"),reduction = "umap_harmony")
FeaturePlot(merged, features = c("TOP2A"),reduction = "umap_harmony")

qc4_a <- FeaturePlot(merged, features = "EPCAM", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("EPCAM")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_a

qc4_b <- FeaturePlot(merged, features = "FN1", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("FN1")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_b

qc4_c <- FeaturePlot(merged, features = "IGHM", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("IGHM")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_c

qc4_d <- FeaturePlot(merged, features = "SPP1", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("SPP1")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_e <- FeaturePlot(merged, features = "CD3E", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("CD3E")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_e

qc4_f <- FeaturePlot(merged, features = "CD4", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("CD4")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_g <- FeaturePlot(merged, features = "CD8A", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("CD8A")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_h <- FeaturePlot(merged, features = "ACTA2", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("ACTA2")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_i <- FeaturePlot(merged, features = "IGHG1", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("IGHG1")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")


Figure_s4 <- plot_grid( qc4_a, qc4_b, qc4_c, qc4_d, qc4_e, qc4_f, qc4_g, qc4_h,qc4_i, labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"), label_size = 8, 
                        align = "hv", axis = "bt", nrow = 3)
ggsave(plot = Figure_s4,"./Results_paper/Plots/Figure S4.pdf", width = 170, height = 130, units = "mm")

# Make a plot of the immune, epi and stromal compartments
merged <- qread("./Intermediate/lognorm_merged_integrated_annotated.qs")
# Griffith data generated in 'Process CosMx V2 only.R'
merged_g <- qread("./Intermediate/lognorm_merged_integrated_annotated_g.qs")

dim(merged)
dim(merged_g)

# Total cells
#219836+171149

# Get some stats on the cells
mean(merged$nCount_RNA)
mean(merged$nFeature_RNA)
mean(merged_g$nCount_RNA)
mean(merged_g$nFeature_RNA)

# Make some QC plots of total UMIs vs unique to give an idea of saturation
set.seed(42)
sampled <- sample(nrow(merged@meta.data),10000)
merged_sampled_sat <- merged@meta.data[sampled,] %>%
  select(nCount_RNA, nFeature_RNA, Slide_short, Manual_toplevel_pred)

set.seed(42)
sampled <- sample(nrow(merged_g@meta.data),10000)
merged_g_sampled_sat <- merged_g@meta.data[sampled,] %>%
  select(nCount_RNA, nFeature_RNA, Slide_short, Manual_toplevel_pred)%>%
  bind_rows(merged_sampled_sat)%>%
  mutate(Slide_short = paste0("Slide ", Slide_short))

# This plot shows we have not hit saturation
# Plamsa cells also tend to be outliers in terms of total probes detected
plt <- ggplot(data = merged_g_sampled_sat, aes(x = nCount_RNA, y = nFeature_RNA, colour = Manual_toplevel_pred))+
  geom_point(size =0.5)+
  facet_wrap(~Slide_short)+
  blank_theme+
  labs(x = "Total probes detected", y = "Unique probes detected", colour  = "Cell type")+
  scale_colour_manual(values =cell_type_colors)

plt

ggsave(plot = plt,"./Results_paper/Plots/Sampled per-slide saturation.png", width = 170, 
       height = 150, units = "mm", dpi = 600)

ggsave(plot = plt,"./Results_paper/Plots/Sampled per-slide saturation.pdf", width = 170, 
       height = 150, units = "mm")
 

# Min dist
# This controls how tightly the embedding is allowed compress points together. 
# Larger values ensure embedded points are more evenly distributed, while smaller 
# values allow the algorithm to optimise more accurately with regard to local structure.
# Spread- The effective scale of embedded points. In combination with min.dist this determines how clustered/clumped the embedded points are.
merged_umap_params <- RunUMAP(merged, reduction="harmony", dims=1:30, reduction.name="umap_harmony", 
                              min.dist = 0.01, spread = 5)

merged_g_umap_params <- RunUMAP(merged_g, reduction="harmony", dims=1:30, reduction.name="umap_harmony", 
                              min.dist = 0.01, spread = 5)

# These UMAPs do look slighlty better but they also highlight segmentation issues a bit more
umap_better <- DimPlot(merged_umap_params, reduction="umap_harmony", group.by="Manual_toplevel_pred", label = T, repel = T, cols = cell_type_colors)
#DimPlot(merged_umap_params, reduction="umap_harmony", group.by="MMR_pred_lowlevel", label = T, label.size = 1, repel = T)+NoLegend()
umap_better_g <-DimPlot(merged_g_umap_params, reduction="umap_harmony", group.by="Manual_toplevel_pred", label = T, repel = T, cols = cell_type_colors)
ggsave(plot = umap_better,"./Results_paper/Extra_plots/TAP_better_UMAP.pdf", width = 170, 
       height = 150, units = "mm")

ggsave(plot = umap_better_g,"./Results_paper/Extra_plots/Griffith_better_UMAP.pdf", width = 170, 
       height = 150, units = "mm")

Figure_1_a <- DimPlot(merged, reduction="umap_harmony", group.by="Manual_toplevel_pred", label = T, label.size = 1, repel = T, cols = cell_type_colors)+
  blank_theme+
  labs(title = "Cell type prediction", x = NULL, y = NULL,
       colour = NULL)+
  theme(legend.key.size = unit(2, 'mm'))+
  guides(color = guide_legend(override.aes = list(size = 0.5))) 

Figure_1_a
#qsave(merged, "./Intermediate/lognorm_merged_integrated_annotated.qs")
#merged <- qread("./Intermediate/lognorm_merged_integrated_annotated.qs")

# Change the NAs to 0 for missing stains
# merged$Mean.CD68 <- replace(merged$Mean.CD68, is.na(merged$Mean.CD68 ), 0)
# merged$Mean.MembraneStain <- replace(merged$Mean.MembraneStain, is.na(merged$Mean.MembraneStain ), merged$Mean.MembraneStain_B2M [is.na(merged$Mean.MembraneStain )])
# merged$Mean.CD3 <- replace(merged$Mean.CD3, is.na(merged$Mean.CD3 ), 0)
Figure_1_c_1 <- FeaturePlot(merged, features = c("Mean.PanCK"), max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+ 
  labs(colour = "PanCK", title = NULL, x = NULL, y = NULL)

Figure_1_c_2 <- FeaturePlot(merged, features = c("Mean.DAPI"), max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "DAPI", title = NULL, x = NULL, y = NULL)

Figure_1_c_3 <- FeaturePlot(merged, features = c("Mean.CD45"), max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "CD45", title = NULL, x = NULL, y = NULL)

cd68 <- merged[,!is.na(merged$Mean.CD68)]

Figure_1_c_4 <- FeaturePlot(cd68, features = c("Mean.CD68"), reduction = "umap_harmony",max.cutoff = "q95", raster = T)+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "CD68", title = NULL, x = NULL, y = NULL)

rm(cd68)

mem <- merged[,!is.na(merged$Mean.MembraneStain)]

Figure_1_c_5 <- FeaturePlot(mem, features = c("Mean.MembraneStain"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "Membrane", title = NULL, x = NULL, y = NULL)

mem_b2m <- merged[,!is.na(merged$Mean.MembraneStain_B2M)]

cd3 <- merged[,!is.na(merged$Mean.CD3)]

Figure_1_c_6 <- FeaturePlot(cd3, features = c("Mean.CD3"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "CD3", title = NULL, x = NULL, y = NULL)

Figure_1_c_6

Figure_1_c_7 <- FeaturePlot(mem_b2m, features = c("Mean.MembraneStain_B2M"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "Membrane (B2M)", title = NULL, x = NULL, y = NULL)

rm(mem_b2m)

Figure_1_c_7

# Assuming we just go with MMR midlevel
cell_type_counts <- merged@meta.data %>%
  group_by(Sample, Sample_type, Manual_toplevel_pred, Slide)%>%
  summarise(count = n())%>%
  group_by(Sample)%>%
  mutate(total = sum(count))%>%
  ungroup()%>%
  mutate(Pct = count/total*100)%>%
  arrange(-Pct)%>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))%>%
  dplyr::rename(Count = count, Total_cells_per_sample = total, Percent_total_cells = Pct)

# Cell type composition plot
Figure_1_b <- ggplot(data = cell_type_counts, aes(x = Sample, y = Percent_total_cells, fill =Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  facet_grid(. ~Sample_type, scales = "free_x", space='free') +
  labs(x = "Sample", y = "% of total cells", fill = "Cell type")+
  blank_theme+
  scale_fill_manual(values =cell_type_colors)+
  # Single column legend
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(2, 'mm'),
        axis.title.x = element_text(vjust = -1))

# Test out propellor or the dataset to get some DA results 
# https://bioconductor.org/packages/release/bioc/vignettes/speckle/inst/doc/speckle.html#test-for-differences-in-cell-line-proportions-in-the-three-technologies

# Set up the data from propellor
merged$cluster <- merged$Manual_toplevel_pred
merged$sample <- merged$Sample 
merged$group <- merged$Sample_type

merged$group <-  factor(merged$group, levels = c("N", "T")) 

Idents(merged) <- merged$cluster
# Save the object
# prop_results <- propeller(merged)
# 
# props <- getTransformedProps(merged$cluster, merged$sample, transform="logit")
# barplot(props$Proportions, col = cell_type_colors,legend=TRUE, 
#         ylab="Proportions")

# Run sccomp as well
sc_result <- merged |>
  sccomp_glm( 
    formula_composition = ~ group, 
    .sample = Sample,
    .cell_group = Manual_toplevel_pred, 
    bimodal_mean_variability_association = T,
    cores = 10
  )

plots <- sccomp_test(sc_result) 

# Seems to be big differences in B cells and neutrophils
Figure_1_e <- plots |> 
  sccomp_boxplot(factor = "group")

Figure_1_e <- Figure_1_e+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4))

# Read in some image PDFs
Tumour_composite <- "/oldvol/apattison/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Tumour.jpg"
normal_composite <- "/oldvol/apattison/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Normal.jpg"
tumour_cells_f <- "/oldvol/apattison/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Figure_1_cells_tumour.png"
normal_cells_f <- "/oldvol/apattison/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Figure_1_cells_normal.png"

tumour_comp <- ggdraw() + 
  draw_image(
    Tumour_composite, scale = 1, width = 1.5,halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

# Run GC to clear the rsession image cache
gc()
normal_comp <- ggdraw() + 
  draw_image(
    normal_composite,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

normal_cells <- ggdraw() + 
  draw_image(
    normal_cells_f,scale = 1,width = 1.5,halign = 0.2, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

tumour_cells <- ggdraw() + 
  draw_image(
    tumour_cells_f,scale = 1,width = 1.5,halign = 0.2, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

cluster_markers <- read_csv("./Results_paper/Tables/V1 cell type markers.csv")%>%
  arrange(-avg_log2FC)%>%
  group_by(cluster)%>%
  mutate(Order = 1:n())%>%
  # Take the top 3 marker genes of each cluster
  filter(Order < 4)%>%
  ungroup()%>%
  arrange(cluster)

DimPlot(merged, cols = cell_type_colors, reduction = "umap_harmony")

Idents(merged) <- merged$seurat_clusters

seu_clust_marks <-  FindAllMarkers(merged)

write_csv(seu_clust_marks, "./Results_paper/Tables/V1 cell type markers seurat clusters.csv")

DimPlot(merged, reduction = "umap_harmony", label = T)
DimPlot(merged, reduction = "umap_harmony", group.by = "MMR_pred_lowlevel", label = T)+NoLegend()
DimPlot(merged, reduction = "umap_harmony", group.by = "seurat_clusters", label = T)

Figure_1_dot <- DotPlot(object = merged, features = unique(cluster_markers$gene), group.by = "cluster",dot.scale = 3)+
  blank_theme+
  labs(y = "Cell type", x = "Probe")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressing"))+
  guides(color = guide_colorbar(title = "Expression"))

Figure_1_dot

Epi <- merged[, merged$Manual_toplevel_pred %in% c("Epi")]

# Take a look at the object
Epi
# Run UMAP clustering on the dataset
Epi <- FindNeighbors(Epi, reduction="harmony", dims=1:30)
Epi <- FindClusters(Epi, cluster.name = "harmony_clusters")
Epi <- RunUMAP(Epi, reduction="harmony", dims=1:30, reduction.name="umap_harmony")

Figure_1_g_i <- DimPlot(Epi, label = T,reduction = "umap_harmony", group.by = "MMR_pred_midlevel", 
                        cols = cell_type_colors, label.size = 2, raster = T)+
  labs(title = "Epi", x = NULL, y = NULL)+
  blank_theme+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  NoLegend()

DimPlot(Epi, label = T,reduction = "umap_harmony")
FeaturePlot(Epi, features = "EPCAM", reduction = "umap_harmony")

# Run with a stromal population
stromal <- merged[, merged$MMR_pred_midlevel %in% c("Fibro", "Endo", "Peri", "SmoothMuscle", "Schwann")]
# Take a look at the object
stromal
# Run UMAP clustering on the dataset
stromal <- FindNeighbors(stromal, reduction="harmony", dims=1:30)
stromal <- FindClusters(stromal, cluster.name = "harmony_clusters")
stromal <- RunUMAP(stromal, reduction="harmony", dims=1:30, reduction.name="umap_harmony")

DimPlot(stromal, label = T,reduction = "umap_harmony")
strom_marks <- FindAllMarkers(stromal)
DimPlot(stromal, label = T,reduction = "umap_harmony", group.by = "MMR_pred_midlevel")

Figure_1_g_ii <- DimPlot(stromal, label = T,reduction = "umap_harmony", 
                         group.by = "MMR_pred_midlevel", cols = cell_type_colors,
                         label.size = 2, raster = T)+
  labs(title = "Stromal", x = NULL, y = NULL)+
  blank_theme+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  NoLegend()
  
FeaturePlot(stromal, features = c("EPCAM", "CCL19", "PECAM1"), reduction = "umap_harmony")

immune <- merged[,! merged$MMR_pred_midlevel %in% c("Fibro", "Endo", "Peri", "SmoothMuscle", "Schwann", "Epi")]

# Take a look at the object
immune
# Run UMAP clustering on the dataset
immune <- FindNeighbors(immune, reduction="harmony", dims=1:30)
immune <- FindClusters(immune, cluster.name = "harmony_clusters")
immune <- RunUMAP(immune, reduction="harmony", dims=1:30, reduction.name="umap_harmony")

Figure_1_g_iii <- DimPlot(immune, label = T,reduction = "umap_harmony", 
                          group.by = "MMR_pred_midlevel", cols = cell_type_colors, label.size = 2,
                          raster = T)+
  labs(title = "Immune", x = NULL, y = NULL)+
  blank_theme+
  guides(color = guide_legend(override.aes = list(size = 0.5)))+
  NoLegend()
  
Figure_1_g_iii

Figure_1_g <- plot_grid(Figure_1_g_i,Figure_1_g_ii,Figure_1_g_iii, nrow = 1)

# Make the last 2 parts of figure 1 a supp figure
Figure_celltypes <- plot_grid(Figure_1_dot, Figure_1_g, labels = c("a", "b"), label_size = 8,nrow = 2)
Figure_celltypes

ggsave(plot = Figure_celltypes,"./Results_paper/Plots/Figure celltypes.png", width = 170, 
       height = 120, units = "mm", dpi = 600)
ggsave(plot = Figure_celltypes,"./Results_paper/Plots/Figure celltypes.pdf", width = 170, 
       height = 120, units = "mm")

left <- plot_grid(Figure_1_a,Figure_1_b, labels = c("a", "c"), label_size = 8, 
                  rel_heights = c(1.6,1), nrow = 2)
right <- plot_grid(Figure_1_c_1, Figure_1_c_2, Figure_1_c_3, Figure_1_c_4, Figure_1_c_6,Figure_1_c_5, labels = c("b"), label_size = 8, 
                   ncol = 2, scale = 1, align = "hv",axis = "bt")
top <- plot_grid(left, right, rel_widths = c(1,1.1))
middle_left <- plot_grid(normal_cells, normal_comp, tumour_cells, tumour_comp, nrow = 2, labels = "d", label_size = 8)
middle_right <- plot_grid(Figure_1_e, nrow = 1, labels = c("e"),label_size = 8)
middle <- plot_grid(middle_left,middle_right, rel_widths = c(1,1.2))
Figure_1 <- plot_grid(top, middle, nrow = 2, rel_heights = c(1.34,1))
ggsave(plot = Figure_1,"./Results_paper/Plots/Figure 1.png", width = 170, 
       height = 170, units = "mm", dpi = 600)
ggsave(plot = Figure_1,"./Results_paper/Plots/Figure 1.pdf", width = 170, 
       height = 170, units = "mm")


# Repeat for Griffith data
Figure_1_a <- DimPlot(merged_g, reduction="umap_harmony", group.by="Manual_toplevel_pred", label = T, label.size = 1, repel = T, cols = cell_type_colors)+
  blank_theme+
  labs(title = "Cell type prediction", x = NULL, y = NULL,
       colour = NULL)+
  theme(legend.key.size = unit(2, 'mm'))+
  guides(color = guide_legend(override.aes = list(size = 0.5))) 

Figure_1_c_1 <- FeaturePlot(merged_g, features = c("Mean.PanCK"), max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+ 
  labs(colour = "PanCK", title = NULL, x = NULL, y = NULL)

Figure_1_c_2 <- FeaturePlot(merged_g, features = c("Mean.DAPI"), max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "DAPI", title = NULL, x = NULL, y = NULL)

Figure_1_c_3 <- FeaturePlot(merged_g, features = c("Mean.CD45"), max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "CD45", title = NULL, x = NULL, y = NULL)

Figure_1_c_4 <- FeaturePlot(merged_g, features = c("Mean.CD68"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "CD68", title = NULL, x = NULL, y = NULL)

Figure_1_c_5 <- FeaturePlot(merged_g, features = c("Mean.Membrane"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "Membrane", title = NULL, x = NULL, y = NULL)

Figure_1_c_6 <- FeaturePlot(merged_g, features = c("Area.um2"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "Area.um2", title = NULL, x = NULL, y = NULL)

# Assuming we just go with MMR midlevel
cell_type_counts <- merged_g@meta.data %>%
  group_by(Sample, Sample_type, Manual_toplevel_pred, Slide)%>%
  summarise(count = n())%>%
  group_by(Sample)%>%
  mutate(total = sum(count))%>%
  ungroup()%>%
  mutate(Pct = count/total*100)%>%
  arrange(-Pct)%>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))%>%
  dplyr::rename(Count = count, Total_cells_per_sample = total, Percent_total_cells = Pct)

# Cell type composition plot
Figure_1_b <- ggplot(data = cell_type_counts, aes(x = Sample, y = Percent_total_cells, fill =Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  facet_grid(. ~Sample_type, scales = "free_x", space='free') +
  labs(x = "Sample", y = "% of total cells", fill = "Cell type")+
  blank_theme+
  scale_fill_manual(values =cell_type_colors)+
  # Single column legend
  #guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(2, 'mm'),
        axis.title.x = element_text(vjust = -1))

# Run sccomp
merged_g$cluster <- merged_g$Manual_toplevel_pred
merged_g$sample <- merged_g$Sample 
merged_g$group <- merged_g$Sample_type

merged_g$group <-  factor(merged_g$group, levels = c("N", "T")) 

Idents(merged_g) <- merged_g$cluster

sc_result <- merged_g |>
  sccomp_glm( 
    formula_composition = ~ group, 
    .sample = Sample,
    .cell_group = Manual_toplevel_pred, 
    bimodal_mean_variability_association = T,
    cores = 10
  )

# Plot the results
plots <- sccomp_test(sc_result) 

# Seems to be big differences in B cells and neutrophils
Figure_1_e <- plots |> 
  sccomp_boxplot(factor = "group")

Figure_1_e <- Figure_1_e+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4))

left <- plot_grid(Figure_1_a,Figure_1_b, labels = c("a", "b"), label_size = 8,
                  rel_heights = c(1.6,1), nrow = 2)
right <- plot_grid(Figure_1_c_1, Figure_1_c_2, Figure_1_c_3, Figure_1_c_4, Figure_1_c_5, Figure_1_c_6, labels = c("c"), label_size = 8,
                   ncol = 2, scale = 1, align = "hv",axis = "bt")
top <- plot_grid(left, right, rel_widths = c(1,1.1))
#bottom_left <- plot_grid(normal_cells, normal_comp, tumour_cells, tumour_comp, nrow = 2, labels = "c", label_size = 8)
Figure_1 <- plot_grid(top, Figure_1_e, nrow = 2, rel_heights = c(1.34,1))
ggsave(plot = Figure_1,"./Results_paper/Plots_G/Figure 1.png", width = 170,
       height = 170, units = "mm", dpi = 600)
ggsave(plot = Figure_1,"./Results_paper/Plots_G/Figure 1.pdf", width = 170,
       height = 170, units = "mm")


# Run niche analysis

# Load in the nanostring data again
g1 <- LoadNanostring.FIX (data.dir = "./raw/10501A/",fov = "fov_all")
g2 <- LoadNanostring.FIX (data.dir = "./raw/120131/",fov = "fov_all")
g3 <- LoadNanostring.FIX (data.dir = "./raw/148191A/",fov = "fov_all")
g4 <- LoadNanostring.FIX (data.dir = "./raw/3461G/",fov = "fov_all")
Tumor_B <- LoadNanostring(data.dir = "./raw/Run5629_Tumor_B/", fov = "fov_all")
Tumor_A <- LoadNanostring(data.dir = "./raw/Run5654_Tumor_A/", fov = "fov_all")
Normal_A <- LoadNanostring(data.dir = "./raw/Run5654_Normal_A/", fov = "fov_all")
new_1 <- LoadNanostring(data.dir = "./raw/Run5850_2186_1/",fov = "fov_all")
new_2 <- LoadNanostring(data.dir = "./raw/Run5850_2186_2/",fov = "fov_all")

merged_g <- qread("./Intermediate/lognorm_merged_integrated_annotated_g.qs")

# Grab the metadata from the integrated object
md_all <- bind_rows(merged@meta.data, merged_g@meta.data)

# Check the slide names for the next steps
unique(md_all$Slide)

slidelist <- list(g1, g2, g3, g4, Tumor_B, Tumor_A, Normal_A, new_1, new_2)
slidenames <- c("10501A_", "120131_", "148191A_", "3461G_", "Tumour_B_","Tumour_A_", "Normal_A_", "Run5850.2186.1_", "Run5850.2186.2_")

slidelist_knn <- list()
for(i in 1:length(slidelist)){
  
  slide <- slidelist[[i]]
  
  slide_name <- slidenames[i]
  
  # Grab the knn matrix from each slide
  slidecounts <- grab_nichecounts(slide, slide_name, md_all, k_param = 50)
  
  # Get the niche object
  slidelist_knn[[i]] <- slidecounts[[2]]
  
}
slidelist_knn[[1]][1:5,1:5]

rownames(slidelist_knn[[1]]) == rownames(slidelist_knn[[2]])

# Loop over and match up rownames to the first object
for(i in 1:length(slidelist_knn)){
  
  slidelist_knn[[i]] <- slidelist_knn[[i]][rownames(slidelist_knn[[1]]),]
  
}

rownames(slidelist_knn[[1]]) == rownames(slidelist_knn[[2]])

# Combine the results into a KNN matrix
combined_niche <- do.call(cbind, slidelist_knn)

dim(combined_niche)

# Save the object
qsave(combined_niche, "./Intermediate/Combined_niche_knn_all.qs")
#combined_niche <- qread("./Intermediate/Combined_niche_knn_all.qs")
hist(combined_niche)

unique(md_all$Slide)

sampled <- sample(ncol(combined_niche), 10000)

combined_niche_samp <- combined_niche[,sampled]

set.seed(42)

# Try flowSOM to see if I can see neighbourhoods
fSOM <- FlowSOM(t(combined_niche),
                # Input options:
                compensate = FALSE,
                transform = FALSE,
                scale = F,
                # Metaclustering options:
                nClus = 9,
                seed = 42)

# Plot a summary of the data
FlowSOMmary(fSOM, plotFile = "./Results_paper/Extra_plots/FlowSOMmary.pdf")

# Each cell commes out mapped to its 50 nearest neighbors (in terms of actual location)
combined_niche[1:5,1:5]
colSums(combined_niche[,1:5])

combined_niche_seu <- CreateSeuratObject(combined_niche)
combined_niche_seu <- NormalizeData(combined_niche_seu)

combined_niche_seu <- FindVariableFeatures(combined_niche_seu)

top10 <- head(VariableFeatures(combined_niche_seu), 19)

plot1 <- VariableFeaturePlot(combined_niche_seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

combined_niche_seu <- ScaleData(combined_niche_seu)
combined_niche_seu <- RunPCA(combined_niche_seu, features = VariableFeatures(object = combined_niche_seu))
DimPlot(combined_niche_seu, reduction = "pca") + NoLegend()

# Scale the data before KNN
scaled <- ScaleData(combined_niche)
# Sample the object for the knn estimation
set.seed(42)

# Sample about 10% of cells
sampled <- sample(ncol(scaled), 1000)
# Estimate the optimal number of niches
# wss" (for total within sum of square) 

Figure_2_c <- fviz_nbclust(t(scaled[,sampled]), kmeans, method = "wss", k.max = 20)

Figure_2_c <- print(Figure_2_c)+
  blank_theme

Figure_2_c <- Figure_2_c+
  labs(y = "Total WSS")

#fviz_nbclust(t(scaled[,sampled]), kmeans, method = "silhouette", k.max = 20)

# It is hard to see a clear elbow
# https://www.statology.org/elbow-method-in-r/
niches.k <-9

set.seed(42)

# Kmeans is fast
# Set 20 iterations to hopfeully converge
results <- kmeans(x = t(scaled), centers = niches.k, iter.max = 20,
                  nstart = 30)

#results$centers
niche_df <- data.frame(Barcode = colnames(scaled), niches = results$cluster)%>%
  mutate(niches = as.character(niches))

# Map the niches back onto the merged integrated object
niche_match <- md_all%>%
  #select(-niches)%>%
  left_join(niche_df)

nichetab <- niche_match%>%
  group_by(niches)%>%
  mutate(Total = n())%>%
  group_by(Manual_toplevel_pred, niches, Total, Sample_type)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Niche_pct = Count/Total *100)%>%
  arrange(niches)%>%
  mutate(niches = factor(niches, levels = unique(niches)))%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels= (cell_types)))

ggplot(data = nichetab, aes(x = niches, y = Niche_pct, fill = Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = cell_type_colors)+
  labs(x = "Niche", y = "Celltype % of niche", fill = "Cell type")+
  blank_theme+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Reset the niches have more descriptive names
niche_match <- md_all%>%
  #select(-niches)%>%
  left_join(niche_df)%>%
  mutate(niches = replace(niches, niches == 1, "Normal like"))%>%
  mutate(niches = replace(niches, niches == 2, "Granulocyte rich"))%>%
  mutate(niches = replace(niches, niches == 3, "Macrophage rich"))%>%
  mutate(niches = replace(niches, niches == 4, "Normal like 2"))%>%
  mutate(niches = replace(niches, niches == 5, "SMC epithelial rich"))%>%
  mutate(niches = replace(niches, niches == 6, "Vasculature"))%>%
  mutate(niches = replace(niches, niches == 7, "TLS"))%>%
  mutate(niches = replace(niches, niches == 8, "Stroma"))%>%
  mutate(niches = replace(niches, niches == 9, "Normal immune infil"))

# Sae the metadata
niche_match%>%
  write_csv("./Results_paper/Tables/High quailty cells metadata.csv")

unique(niche_match$niches)

nichetab <- niche_match%>%
  group_by(niches)%>%
  mutate(Total = n())%>%
  group_by(Manual_toplevel_pred, niches, Total)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Niche_pct = Count/Total *100)%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels= (cell_types)))

# Cellular niche composition
Figure_2_a <- ggplot(data = nichetab, aes(x = niches, y = Niche_pct, fill = Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = cell_type_colors)+
  labs(x = "Niche", y = "Celltype % of niche", fill = "Cell type")+
  blank_theme+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(2, 'mm'))

Figure_2_a

# Map the lower quality samples back to a cell/niche

# Read in the probe for both datasets
all_names <- read_csv("./Intermediate/Probenames.csv")

merged_orig_all <- qread("./Intermediate/Merged_raw_unfiltered.qs")
merged_g_all <- qread("./Intermediate/Merged_raw_unfiltered_g.qs")

md_combined <- read_csv("./Results/Tables/all_metadata.csv")

sample_slide_summary <- md_combined%>%
  group_by(Sample, Slide)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  arrange(Slide, Sample)%>%
  mutate(Slide_sample = paste0(Slide, "-", Sample))%>%
  write_csv("./Results/Tables/Sample_slide.csv")

merged_orig_all@meta.data <- merged_orig_all@meta.data%>%
  rownames_to_column("Barcode")%>%
  left_join(md_combined)

merged_g_all@meta.data <- merged_g_all@meta.data%>%
  rownames_to_column("Barcode")%>%
  left_join(md_combined)

rownames(merged_orig_all@meta.data ) <- merged_orig_all$Barcode
rownames(merged_g_all@meta.data ) <- merged_g_all$Barcode

# Get the number of cells per FOV
cells_per_fov <- md_combined%>%
  group_by(Slide, fov)%>%
  summarise(Unfiltered_total_cells = n())%>%
  write_csv("./Results_paper/Tables/Cells_per_FOV.csv")
  
# Annotate the HQ objects with niches
niche_only <- niche_match%>%
  select(Barcode, niches)

merged_g@meta.data <- left_join(merged_g@meta.data, niche_only)
merged@meta.data <- left_join(merged@meta.data, niche_only)

# Reset rownames
rownames(merged_g@meta.data) <-merged_g@meta.data$Barcode
rownames(merged@meta.data) <-merged@meta.data$Barcode

# Create a cell/niche label 
merged$cell_niche <- paste0(merged$Manual_toplevel_pred, "_", merged$niches)
merged_g$cell_niche <- paste0(merged_g$Manual_toplevel_pred, "_", merged_g$niches)

# Get only the low quality cells from the TAP run
low_qual <- merged_orig_all[,!merged_orig_all$Barcode %in% merged$Barcode]

# Run singleR for the low quailty cells to give them an annotaion to visualise
predictions <- SingleR(test=low_qual@assays$RNA$counts,
                       ref=merged@assays$RNA$counts, 
                       labels=merged$cell_niche, 
                       num.threads = 30,
                       aggr.ref = T)

# Add a prediction to the low quailty cells
low_qual$cell_niche <- predictions$labels

# Read the annotation of the high quailty cells
hq_labs <- merged@meta.data%>%
  select(Barcode, cell_niche)

lq_labs <- low_qual@meta.data%>%
  select(Barcode,cell_niche_match= cell_niche)

# Fill in the low quality cell annotations
anno_matched <- merged_orig_all@meta.data%>%
  left_join(hq_labs)%>%
  left_join(lq_labs)%>%
  mutate(cell_niche = replace(cell_niche, is.na(cell_niche), cell_niche_match[is.na(cell_niche)]))%>%
  mutate(Manual_toplevel_pred = gsub("_.*", "", cell_niche))%>%
  # Non greedy match
  mutate(niches = gsub("^.*?_", "", cell_niche))

# Add the annos on
merged_orig_all$Manual_toplevel_pred <- anno_matched$Manual_toplevel_pred
merged_orig_all$niches <- anno_matched$niches

# All matched successfully 
sum(is.na(merged_orig_all$niches))

slide_short <- data.frame(Slide= merged_orig_all$Slide)%>%
  left_join(slide_alias)

merged_orig_all$Slide_short <- slide_short$slide_names
merged_orig_all$Sample_Slide <- paste0(merged_orig_all$Sample, " ", merged_orig_all$Slide_short)

# Add on T/N and Donor annotation 
merged_orig_all$Sample_type <- ifelse(grepl("T", merged_orig_all$Sample), "T", "N")
merged_orig_all$Donor <- gsub("T|N","", merged_orig_all$Sample)
merged_orig_all$Donor <- paste0("Donor", merged_orig_all$Donor)

# I think a min cutoff of 20 mols is at least required to have some
# idea of cell type
merged_orig_all$Manual_toplevel_pred <- replace(merged_orig_all$Manual_toplevel_pred, merged_orig_all$nCount_RNA <20, "QC_fail")
merged_orig_all$niches <- replace(merged_orig_all$niches, merged_orig_all$nCount_RNA <20, "QC_fail")

# Save the annotated objects
qsave(merged_g, "./Intermediate/lognorm_merged_integrated_annotated_g.qs")
qsave(merged, "./Intermediate/lognorm_merged_integrated_annotated.qs")

# Save the annotated all cells object for the TAP
qsave(merged_orig_all, "./Intermediate/All_cells_orig_annotated.qs")

# Get only the low quality cells from the Griffith run
low_qual <- merged_g_all[,!merged_g_all$Barcode %in% merged_g$Barcode]

# Run singleR for the low quailty cells to give them an annotaion to visualise
predictions <- SingleR(test=low_qual@assays$RNA$counts,
                       ref=merged_g@assays$RNA$counts, 
                       labels=merged_g$cell_niche, num.threads = 30,
                       aggr.ref = T) 

# Add a prediction to the low quailty cells
low_qual$cell_niche <- predictions$labels

# Read the annotation of the high quailty cells
hq_labs <- merged_g@meta.data%>%
  select(Barcode, cell_niche)

lq_labs <- low_qual@meta.data%>%
  select(Barcode,cell_niche_match= cell_niche)

anno_matched <- merged_g_all@meta.data%>%
  left_join(hq_labs)%>%
  left_join(lq_labs)%>%
  mutate(cell_niche = replace(cell_niche, is.na(cell_niche), cell_niche_match[is.na(cell_niche)]))%>%
  mutate(Manual_toplevel_pred = gsub("_.*", "", cell_niche))%>%
  # Non greedy match
  mutate(niches = gsub("^.*?_", "", cell_niche))

merged_g_all$Manual_toplevel_pred <- anno_matched$Manual_toplevel_pred
merged_g_all$niches <- anno_matched$niches

# All matched successfully 
sum(is.na(merged_g_all$niches))

slide_short <- data.frame(Slide= merged_g_all$Slide)%>%
  left_join(slide_alias)

merged_g_all$Slide_short <- slide_short$slide_names
merged_g_all$Sample_Slide <- paste0(merged_g_all$Sample, " ", merged_g_all$Slide_short)

qsave(merged_orig_all, "./Intermediate/All_cells_g_annotated.qs")

# Add on T/N and Donor annotation 
merged_g_all$Sample_type <- ifelse(grepl("T", merged_g_all$Sample), "T", "N")
merged_g_all$Donor <- gsub("T|N","", merged_g_all$Sample)
merged_g_all$Donor <- paste0("Donor", merged_g_all$Donor)

# I think a min cutoff of 20 mols is at least required to have some
# idea of cell type
merged_g_all$Manual_toplevel_pred <- replace(merged_g_all$Manual_toplevel_pred, merged_g_all$nCount_RNA <20, "QC_fail")
merged_g_all$niches <- replace(merged_g_all$niches, merged_g_all$nCount_RNA <20, "QC_fail")

nichetab <- niche_match%>%
  group_by(Sample, Slide_short, fov)%>%
  mutate(Total = n())%>%
  mutate(fov = as.character(fov))%>%
  mutate(niches = as.character(niches))%>%
  group_by(fov, Sample, niches, Total, Slide_short)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Niche_pct = Count/Total *100)%>%
  mutate(fov = as.numeric(fov))%>%
  mutate(Slide_short = as.numeric(Slide_short))%>%
  arrange(Slide_short, fov)%>%
  mutate(Slide_short = factor(Slide_short, levels= unique(Slide_short)))%>%
  mutate(fov = factor(fov, levels= unique(fov)))%>%
  mutate(Sample = factor(Sample, levels= unique(Sample)))

Figure_s5 <- ggplot(data = nichetab, aes(x = fov, y = Niche_pct, fill = niches))+
  geom_bar(stat = "identity")+
  facet_wrap(~Slide_short + Sample, scales  = "free", ncol = 3)+
  scale_fill_manual(values = nichecols)+
  blank_theme+
  theme(legend.key.size = unit(3, 'mm'),
        legend.position = "top")+
  labs(x = "FOV", y = "Niche %", fill = "Niche")+
  ggtitle("Slide/Sample")

ggsave(plot = Figure_s5,"./Results_paper/Plots/Figure S5.pdf", width = 170, height = 180, units = "mm")

md_all <- bind_rows(merged_orig_all@meta.data, merged_g_all@meta.data)%>%
  mutate(Quality = ifelse(Barcode %in% niche_match$Barcode, "High", "Low"))%>%
  write_csv("./Results_paper/Tables/All cells metadata.csv")

# Plot the niches per FOV
nichetab <- md_all%>%
  mutate(fov_numeric = fov)%>%
  mutate(fov = paste0(fov, " ", Sample))%>%
  group_by(Slide, fov)%>%
  mutate(Total = n())%>%
  mutate(fov = as.character(fov))%>%
  mutate(niches = as.character(niches))%>%
  group_by(fov, fov_numeric, niches, Total, Slide,Sample)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Niche_pct = Count/Total *100)%>%
  left_join(slide_alias)%>%
  mutate(slide_numeric = as.numeric(slide_names))%>%
  arrange(slide_numeric,Sample, fov_numeric)%>%
  mutate(fov = factor(fov, levels = unique(fov)))%>%
  mutate(slide_names_long = factor(slide_names_long, levels = unique(slide_names_long)))

# Save the metadata in a readable format for the paper
md_save <- md_all%>%
  mutate(fov = as.character(fov))%>%
  group_by(fov, Slide, Sex, Sample, Sample_type, Donor, Sidedness, `Age at Dx`, `Primary tumour location`, `MMR status`)%>%
  summarise(Count = n(), Mean_genes_detected = mean(nFeature_RNA), Mean_probes_detected = mean(nCount_RNA))%>%
  ungroup()%>%
  dplyr::select(Slide, FOV = fov, Sample, Sample_type, Donor, Sex, 
                Age_at_diagnosis = `Age at Dx`, 
                Primary_tumour_location = `Primary tumour location`,
                MMR_status = `MMR status`, Sidedness, Total_cells = Count,
                Mean_genes_detected, Mean_probes_detected)%>%
  write_csv("./Results_paper/Tables/Slide-FOV sample metadata.csv")
  
plt <- ggplot(data = nichetab, aes(x = fov, y = Niche_pct, fill = niches))+
  geom_bar(stat = "identity")+
  facet_wrap(~slide_names_long, scales  = "free", ncol = 2)+
  scale_fill_manual(values = nichecols)+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'), legend.position = "bottom right")+
  labs(x = "FOV", y = "FOV niche %", fill = "Niche")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5))

#plt

save <- "./Results_paper/Extra_plots/FOV_niches.pdf"
ggsave(filename = save,plot = plt, width = 170, height = 170, units = "mm") 

plt <- ggplot(data = nichetab, aes(x = fov, y = Niche_pct, fill = niches))+
  geom_bar(stat = "identity")+
  facet_wrap(~Slide, scales  = "free", ncol = 1)+
  scale_fill_manual(values = nichecols)+
  blank_theme+
  labs(x = "FOV", y = "FOV niche %", fill = "Niche")+
  ggtitle("FOV/Slide")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  coord_flip()

plt

save <- save <- "./Results_paper/Extra_plots/FOV_niches_long.pdf"
ggsave(filename = save,plot = plt, width = 9, height = 60, limitsize = FALSE) 

# Plot the per-sample niche composition
nichetab_per_sample <- md_all%>%
  group_by(Sample, Sample_type)%>%
  mutate(Total = n())%>%
  mutate(niches = as.character(niches))%>%
  group_by(Sample, niches, Total,Sample_type)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Niche_pct = Count/Total *100)

Figure_2_b <- ggplot(data = nichetab_per_sample, aes(x = Sample, y = Niche_pct, fill = niches))+
  geom_bar(stat = "identity")+
  facet_wrap(~Sample_type, scales = "free_x")+
  scale_fill_manual(values = nichecols)+
  labs(x = "Sample", y = "Niche %", fill = "Niche")+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Read in the image PDFs
Tumour_cells <- "./Results_paper/Extra_plots/Figure_2_image_inputs/iii.png"
Tumour_niche <- "./Results_paper/Extra_plots/Figure_2_image_inputs/iv.png"
Normal_cells <- "./Results_paper/Extra_plots/Figure_2_image_inputs/i.png"
Normal_niche <- "./Results_paper/Extra_plots/Figure_2_image_inputs/ii.png"

gc()
tumour_c <- ggdraw() + 
  draw_image(Tumour_cells, scale = 1, width = 1.5,halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

tumour_n <- ggdraw() + 
  draw_image(Tumour_niche, scale = 1, width = 1.5,halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

gc()
Normal_c <- ggdraw() + 
  draw_image(Normal_cells, scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

Normal_n <- ggdraw() + 
  draw_image(Normal_niche, scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

merged$sample <- paste0("Sample", merged$Sample )

# Run sccomp to see differential niches
sc_result <- niche_match |>
  sccomp_glm( 
    formula_composition = ~ Sample_type, 
    .sample = Sample,
    .cell_group = niches, 
    bimodal_mean_variability_association = T,
    cores = 10
  )

plots <- plot_summary(sc_result) 

# Seems to be big differences in B cells and neutrophils
plots$boxplot

Figure_2_e <- plots$boxplot[[1]]
Figure_2_e <- Figure_2_e+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        legend.position = "top",
        axis.text.y =element_text(size=4))

# Make figure 2
top_left <- plot_grid(Figure_2_a, labels = c("a"), label_size = 8)
top_right <- plot_grid(Figure_2_b, labels = c("b"), label_size = 8, nrow = 1)
top_row <- plot_grid(top_left, top_right, rel_widths = c(1,1.2))
bottom_left <- plot_grid(Normal_c,Normal_n,tumour_c, tumour_n, labels = c("i", "ii", "iii", "iv"),vjust =0.3,  label_size = 6, nrow = 2,greedy =T)
bottom_right <- plot_grid(Figure_2_e, labels = "d", label_size = 8)
bottom_row <- plot_grid(bottom_left, bottom_right, rel_widths = c(1,1.2), labels = c("c"), label_size = 8, hjust = 0, vjust = -1)
Figure_2 <- plot_grid(top_row,bottom_row, nrow = 2, rel_heights = c(1,1.15))
ggsave(plot = Figure_2,"./Results_paper/Plots/Figure 2.png", width = 170, 
       height = 140, units = "mm", dpi = 600)
ggsave(plot = Figure_2,"./Results_paper/Plots/Figure 2.pdf", width = 170, 
       height = 140, units = "mm")


# Try a chi-square approach to describe enrichment of cell type in niches

# Make a matrix of cell types in niches
nichetab <- table( merged$Manual_toplevel_pred, merged$niches)%>%
  as.matrix()

chisq <- chisq.test(nichetab)
chisq

# http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
# Visualize Pearson residuals to see the highest contributers to each niche
save <- "./Results_paper/Extra_plots/Niche_celltype_pearsons_residuals.pdf"
pdf(save, width = 5, height = 6)
corrplot(chisq$residuals, is.cor = FALSE)
dev.off()

# Find markers of each cell type
Idents(merged) <- merged$Manual_toplevel_pred
clustmarks <- FindAllMarkers(merged)
write_csv(clustmarks, "./Results_paper/Tables/V1 cell type markers.csv")

Idents(merged_g) <- merged_g$Manual_toplevel_pred
clustmarks_g <- FindAllMarkers(merged_g)
write_csv(clustmarks_g, "./Results_paper/Tables/V2 cell type markers.csv")

# Get the top 3 markers of each cluster
top_3 <- clustmarks %>%
  arrange(cluster, -avg_log2FC)%>%
  group_by(cluster)%>%
  mutate(Count = 1:n())%>%
  ungroup()%>%
  filter(Count <=3)

marker_genes <- factor(unique(top_3$gene), levels = unique(top_3$gene))

figure_s_0_c <- DotPlot(object = merged, features = marker_genes, group.by = "Manual_toplevel_pred")+
  blank_theme+
  labs(y = "Cell type", x = "Probe")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(4, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))+
  scale_colour_gradient(low = "blue", high = "orange")

# Add on some extra metadata before saving
md_all <- merged@meta.data%>%
  left_join(donor_md)

rownames(md_all) <-md_all$Barcode
merged@meta.data <- md_all

# Make a plot of read counts per cell type
ct_frame <- merged@meta.data%>%
  group_by(Manual_toplevel_pred, Sample,Slide_short, Sample_type)%>%
  summarise(Mean_detection = mean(nCount_RNA))%>%
  group_by(Manual_toplevel_pred)%>%
  mutate(med_mean = mean(Mean_detection))%>%
  ungroup()%>%
  arrange(-med_mean)%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))

figure_s_0_d <- ggplot(data = ct_frame, aes(x = Manual_toplevel_pred, y = Mean_detection, fill = Manual_toplevel_pred))+
  facet_wrap(~Slide_short)+
  blank_theme+
  geom_jitter(aes(colour = Sample), width = 0.3, height = 0, size = 0.3)+
  guides(fill = "none")+
  scale_fill_manual(values = cell_type_colors)+
  scale_colour_manual(values = Sample_cols)+
  labs(x = "Cell type", y = "Mean probe detection per cell")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(4, 'mm'))+
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.7, linewidth = 0.2)

# Make misc figure S0
top_left <- plot_grid(Figure_2_a, labels = c("a"), label_size = 8)
top_right <- plot_grid(Figure_2_b, Figure_2_c, labels = c("b", "c"), label_size = 8, nrow = 2)
top_row <- plot_grid(top_left, top_right, rel_widths = c(1,1))
bottom_left <- plot_grid(tumour_c, tumour_n, Normal_c,Normal_n, labels = c("d"), label_size = 8, align = "h", axis = "bt", nrow = 2)
bottom_right <- plot_grid(Figure_2_e, labels = c("e"), label_size = 8)
bottom_row <- plot_grid(bottom_left, bottom_right, rel_widths = c(1,1))
Figure_s0 <- plot_grid(top_row,bottom_row, nrow = 2, rel_heights = c(1,1.15))
ggsave(plot = Figure_2,"./Results_paper/Plots/Figure 2.png", width = 170, 
       height = 140, units = "mm", dpi = 600)

# Neutrophil story
mmr_atlas$cell_sampletype <- paste0(mmr_atlas$clMidwayPr, " ", mmr_atlas$SPECIMEN_TYPE)

# Make a dotplot of neutrphil chemoattractant genes in the MMR atlas
genes = c("CXCR2", "IL1B",  "CCL3L1", "CCL4", "CCL3L3", "CXCL1","CXCL2", "CXCL3", "CXCL8", "CXCL5","CXCL6", "CXCL12","CXCR4", "IL6", "SPP1", "CCL20")

rownames(mmr_atlas)[grepl("^IL1", rownames(mmr_atlas))]

genes = c("DUOX2", "DUOXA2", "OLFM4", "CXCR1", "CXCR2", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "IL1RN", "IL1R2", "IL1R1", "IL1A", "IL1B", "TNF")

# Make a dot plot of the neutrophil related genes
DotPlot(object = mmr_atlas, features = genes, group.by = "cell_sampletype", dot.min = 0.01)+
  blank_theme+
  labs(y = "Sample + slide", x = "Probes")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.key.size = unit(4, 'mm'))+
  guides(size = guide_legend(title = "%\nexpressed"))+
  guides(color = guide_colorbar(title = "Expression"))

# Get gene expression on each slide to find out where genes are expressed
SLC2A1 <-gene_fov_finder("SLC2A1", merged)
MMP7 <-gene_fov_finder(c("MMP7"), merged)
TIMP1 <-gene_fov_finder("TIMP1", merged)
EPHA2 <-gene_fov_finder("EPHA2", merged)
OLFM4 <-gene_fov_finder("OLFM4", merged)
CD33 <-gene_fov_finder("CD33", merged)
CCL19 <-gene_fov_finder("CCL19", merged)
CXCL8 <-gene_fov_finder("CXCL8", merged)
CXCL10 <- gene_fov_finder("CXCL10", merged)
CXCL14 <- gene_fov_finder("CXCL14", merged)
CXCL13 <- gene_fov_finder("CXCL13", merged_g)
LAMP3 <- gene_fov_finder("LAMP3", merged)
CCR7 <- gene_fov_finder("CCR7", merged)
IFNG <- gene_fov_finder("IFNG", merged)
PDGFRA <- gene_fov_finder("PDGFRA", merged)
TXTN1 <- gene_fov_finder("TXTN1", merged)
TNXAB <- gene_fov_finder("TNXA/B", merged_g)
NRG1 <- gene_fov_finder("NRG1", merged)
RSPO3 <- gene_fov_finder("RSPO3", merged)
LGR5 <- gene_fov_finder("LGR5", merged_g)
PIGR <-gene_fov_finder("PIGR", merged)

rowSums(merged_g@assays$RNA$counts)[grepl("TNXA/B", rownames(merged_g))]

# CXCL8 looks to still be a good marker, even in the new runs
rownames(merged_g_all)[grepl("^TNX", rownames(merged_g_all))]
rownames(merged)[grepl("^TNX", rownames(merged))]

FeaturePlot(merged, features = "COL1A1", reduction = "umap_harmony")

unique(merged$Slide)

neus <- slide_list(merged, "Tumour_A", 8)
CXCL8 <- slide_list(merged, "Run5850.2186.2", 23)

TLS <- slide_list(merged, "Run5850.2186.2", 37)
GLUT1 <- slide_list(merged, "Run5850.2186.1", 35)
TIMP1 <- slide_list(merged, "Run5850.2186.1", 35)
neu <- slide_list(merged, "Tumour_A", 7)
t_active <- slide_list(merged, "Tumour_B", 9)
cxcl14 <- slide_list(merged, "Run5850.2186.2", 30)


# Plot cell size vs gene expression and gene expression in general
# Plot cell size for a supp figure
size_per_sample <- merged@meta.data%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Medain_area_cell = median(Area))%>%
  ungroup()%>%
  arrange(-Medain_area_cell)%>%
  mutate(Slide_short = factor(Slide_short, levels = unique(Slide_short)))

# Looks like there is a conversion from pixels to area
# Probably depends on the imager
merged_g$Area[1:5]
merged_g$Area.um2[1:5]

cor.test(merged_g$Area, merged_g$Area.um2)

size_per_sample_g <- merged_g@meta.data%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Medain_area_cell = median(Area))%>%
  ungroup()%>%
  arrange(-Medain_area_cell)%>%
  mutate(Slide_short = factor(Slide_short, levels = unique(Slide_short)))

size_per_sample <- bind_rows(size_per_sample, size_per_sample_g)%>%
  mutate(Slide_short = paste0("Slide", Slide_short))%>%
  arrange(Slide_short)%>%
  mutate(Slide_short = factor(Slide_short, levels = unique(Slide_short)))%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))

Figure_s6_a <- ggplot(data = size_per_sample, aes(x = Slide_short, y = Medain_area_cell, 
                                                  colour = Manual_toplevel_pred, group = Manual_toplevel_pred))+
  geom_point(size=0.5)+
  geom_line(linewidth=0.5)+
  labs( y = "Median pixel area per cell", colour = "Cell type", x = "Slide")+
  blank_theme+
  scale_colour_manual(values = cell_type_colors)+
  guides(colour = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Figure_s6_a

# Description of the 9 slides
slide_df <- data.frame(Slide = unique(size_per_sample$Slide_short))%>%
  mutate(Protocol = "TAP (V1)")%>%
  mutate(Protocol = replace(Protocol, Slide %in% c("Slide4", "Slide5"), "TAP (V2)"))%>%
  mutate(Protocol = replace(Protocol, Slide %in% c("Slide6", "Slide7","Slide8", "Slide9"), "Griffith (V3)"))%>%
  mutate(Probeset = "V1")%>%
  mutate(Probeset = replace(Probeset, Slide %in% c("Slide6", "Slide7","Slide8", "Slide9"), "V2"))%>%
  mutate(Slide = gsub("Slide", "Slide ", Slide))

rownames(slide_df) <- NULL

table_grob <- tableGrob(slide_df, rows = NULL, theme =  ttheme_minimal(base_size = 5,  padding = unit(c(5, 2), "mm")))

grid.newpage()
grid.draw(table_grob)

# Make a plot of cell size and transcript counts
size_vs_count <- merged@meta.data%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Median_area_cell = median(Area), Median_transcript_counts =  median(nCount_RNA))%>%
  ungroup()%>%
  group_by(Slide_short)%>%
  mutate(cor = cor.test(Median_area_cell, Median_transcript_counts)$estimate)%>%
  ungroup()%>%
  mutate(cor = round(cor,2))

size_vs_count_g <- merged_g@meta.data%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Median_area_cell = median(Area), Median_transcript_counts =  median(nCount_RNA))%>%
  ungroup()%>%
  group_by(Slide_short)%>%
  mutate(cor = cor.test(Median_area_cell, Median_transcript_counts)$estimate)%>%
  ungroup()%>%
  mutate(cor = round(cor,2))

size_vs_count <- bind_rows(size_vs_count, size_vs_count_g)%>%
  mutate(Slide_short = paste0("Slide ", Slide_short, ", r = ", cor))

Figure_s6_b <- ggplot(data = size_vs_count, aes(x = Median_area_cell, y = Median_transcript_counts, 
                                                 colour = Manual_toplevel_pred))+
  geom_point(size = 0.5)+
  facet_wrap(~Slide_short, scales = "free")+
  scale_colour_manual(values = cell_type_colors)+
  blank_theme+
  labs(x = "Median cell area (pixels)", 
       y = "Median probes detected per cell",
       shape = "Slide",
       colour = "Cell type")+
   theme(legend.key.size = unit(2.5, 'mm'),
         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Figure_s6_b

# Make a plot of the most highly detected genes
merged_raw_unfiltered <- qread("./Intermediate/Merged_raw_unfiltered.qs")

md_all <- read_csv("./Results_paper/Tables/All cells metadata.csv")

counts_list <- list()

# Loop over each slide and calculate the probe expression
for(i in 1:5){
  
  slide_filter <- md_all%>%
    filter(Slide_short == i)
  
  slide_counts <- merged_raw_unfiltered[,slide_filter$Barcode]
  
  exprs <- data.frame(Expression = rowSums(slide_counts@assays$RNA$counts))%>%
    rownames_to_column("Gene")%>%
    mutate(Slide = i)
  
  counts_list[[i]] <- exprs
  
}

merged_raw_unfiltered_g <- qread("./Intermediate/Merged_raw_unfiltered_g.qs")
# Loop over each slide and calculate the probe expression
for(i in 6:9){
  
  print(i)
  
  slide_filter <- md_all%>%
    filter(Slide_short == i)
  
  slide_counts <- merged_raw_unfiltered_g[,slide_filter$Barcode]
  
  exprs <- data.frame(Expression = rowSums(slide_counts@assays$RNA$counts))%>%
    rownames_to_column("Gene")%>%
    mutate(Slide = i)
  
  counts_list[[i]] <- exprs
  
}

counts_per_slide <- bind_rows(counts_list)%>%
  group_by(Slide)%>%
  arrange(-Expression)%>%
  mutate(Index = 1:n())%>%
  # Keep the top 10 expressed genes in each slide
  filter(Index <11)%>%
  ungroup()%>%
  mutate(Gene = factor(Gene, levels = unique(Gene)))%>%
  mutate(Slide = paste0("Slide ", Slide))

Figure_s6_c <- ggplot(data = counts_per_slide, aes(x = Gene, y =Expression))+
  geom_bar(stat = "identity")+
  facet_wrap(~Slide, scales = "free")+
  blank_theme+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic", size = 4))+
  labs(x = "Probe", y = "Total cellular probes detected")

Figure_s6_c

# Make a plot of cell size vs count
count_vs_size <- md_all%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Median_area_cell = median(Area), Total_cells_called = n())%>%
  arrange(Total_cells_called)%>%
  ungroup()%>%
  filter(Manual_toplevel_pred != "QC_fail")%>%
  group_by(Slide_short)%>%
  mutate(cor = cor.test(Median_area_cell, Total_cells_called)$estimate)%>%
  ungroup()%>%
  mutate(cor = round(cor,2))%>%
  mutate(Slide_short = paste0("Slide ", Slide_short, ", r = ", cor))

Figure_s6_d <- ggplot(data = count_vs_size, 
                            aes(x = Median_area_cell, y = log2(Total_cells_called), 
                                colour = Manual_toplevel_pred,
                                label = Manual_toplevel_pred))+
  geom_point(size = 0.5)+
  facet_wrap(~Slide_short, scales = "free")+
  scale_colour_manual(values = cell_type_colors)+
  blank_theme+
  labs(x = "Median cell area (pixels)", y = expression('Log'[2]*' total cells annotated'),
       colour = "Cell type")+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(colour = "none")

Figure_s6_d

# Make the figure
top_left <- plot_grid(table_grob, Figure_s6_a, ncol = 1, label_size = 8, labels = c("a", "b"), rel_heights = c(0.7,1))
top <- plot_grid(top_left, Figure_s6_b, labels = c("", "c"), label_size = 8, rel_widths = c(0.8,1.1))
bottom <- plot_grid(Figure_s6_c, Figure_s6_d, labels = c("d", "e"), label_size = 8, rel_widths = c(1.2,0.9))
figure_s_6 <- plot_grid(top, bottom, nrow = 2)
ggsave(plot = figure_s_6,"./Results_paper/Plots/Figure S6.pdf", width = 170, height = 170, units = "mm")
