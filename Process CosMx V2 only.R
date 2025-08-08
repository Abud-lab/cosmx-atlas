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
  library(FlowSOM)
  library(ggrepel)
  library(progressr)
})

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_v2/functions.R") 

# Set a homedir for the analysis so as to not use hardcoded file paths
setwd("/homevol/apattison/Data/Spatial")

# Set output folders
plots <- "./Results_paper/Plots_G/"

system(paste0("mkdir -p ", plots))

md_old <- read_csv("./raw/md_liam.csv")%>%
  mutate(fov = as.numeric(fov))%>%
  group_by(Slide)%>%
  mutate(fov = (fov+1)- min(fov))%>%
  ungroup()%>%
  mutate(Slide = gsub(" ", "_", Slide))%>%
  mutate(Barcode = paste0(Slide, "_", cell_ID, "_", fov))

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
  ungroup()%>%
  mutate(Slide = gsub(" ", "_", Slide))%>%
  mutate(Barcode = paste0(Slide, "_", cell_ID, "_", fov))

# Read in the Griffith metadata
g1_md <- read_csv("./raw/10501A/10501A_metadata_file.csv.gz")%>%
  mutate(Slide = "10501A")%>%
  mutate(Run_Tissue_name = as.character(Run_Tissue_name))

g2_md <- read_csv("./raw/120131/120131_metadata_file.csv.gz")%>%
  mutate(Slide = "120131")%>%
  mutate(Run_Tissue_name = as.character(Run_Tissue_name))

g3_md <- read_csv("./raw/148191A/148191A_metadata_file.csv.gz")%>%
  mutate(Slide = "148191A")%>%
  mutate(Run_Tissue_name = as.character(Run_Tissue_name))

g4_md <- read_csv("./raw/3461G/3461G_metadata_file.csv.gz")%>%
  mutate(Slide = "3461G")%>%
  mutate(Run_Tissue_name = as.character(Run_Tissue_name))

fov_line_up_g <- read_csv("./raw/fov_line_up_g.csv")

# Join all the Griffith samples together 
md_new_g <- bind_rows(g1_md, g2_md, g3_md, g4_md)%>%
  mutate(fov = as.numeric(fov))%>%
  left_join(fov_line_up_g)%>%
  group_by(Slide)%>%
  ungroup()%>%
  mutate(Slide = gsub(" ", "_", Slide))%>%
  mutate(Barcode = paste0(Slide, "_", cell_ID, "_", fov))

# Different numbers of FOVs each time!
unique(md_old$fov)
unique(md_new_all$fov)
unique(md_new_g$fov)

# Combine all the metadata
md_combined <- bind_rows(md_old, md_new_all, md_new_g)%>%
  select(-nCount_RNA, -nFeature_RNA)

g1 <- get_nano_counts_FIX(data.dir = "./raw/10501A/",sample_name = "10501A")
g2 <- get_nano_counts_FIX(data.dir = "./raw/120131/",sample_name = "120131")
g3 <- get_nano_counts_FIX(data.dir = "./raw/148191A/",sample_name = "148191A")
g4 <- get_nano_counts_FIX(data.dir = "./raw/3461G/",sample_name = "3461G")

merged <-  merge(x = g1, y = c(g2, g3, g4))%>%
  JoinLayers()

dim(merged@assays$RNA$counts)

cs <- colSums(merged@assays$RNA$counts)

hist(cs, breaks = 1000, xlim = c(0,100))

dim(merged)

rm(g1, g2, g3, g4)

qsave(merged, "./Intermediate/Merged_raw_unfiltered_g.qs")

#merged <- qread("./Intermediate/Merged_raw_unfiltered_g.qs")

# Make a plot of the top expressed genes
gene_exprs <- rowMeans(merged@assays$RNA$counts)

hist(gene_exprs, breaks = 100)

# Make a plot that compares single probes to single cell gene expression
mmr_atlas <- qread("/homevol/apattison/Data/Reference/GSE178341_lognorm_annotated.qs")

# Get the genes from the MMR atlas
gene_exprs_mmr <- rowMeans(mmr_atlas@assays$RNA@counts)

gene_exprs_df_mmr <- data.frame(gene_exprs_mmr)%>%
  rownames_to_column("Gene")

gene_exprs_df <- data.frame(gene_exprs)%>%
  rownames_to_column("Gene")%>%
  arrange(gene_exprs)%>%
  mutate(nums = 1:n())%>%
  left_join(gene_exprs_df_mmr)%>%
  mutate(labs = replace(Gene, nums < 1150, NA))

# Plot just the most expressed genes
gene_exprs_df_top <- gene_exprs_df%>%
  filter(nums > 1150)

# Plot the most highly expressed probes
figure_s_0_a <- ggplot(data = gene_exprs_df_top, aes(x = nums, y= gene_exprs, label = labs))+
  geom_point()+
  geom_text_repel(max.overlaps = 20,force = 10, colour = "blue", fontface = "italic")+
  labs(x = "Probe index", y = "Mean probes per cell")+
  blank_theme

figure_s_0_a

cor <- cor.test(gene_exprs_df$gene_exprs, gene_exprs_df$gene_exprs_mmr)

figure_s_0_b <- ggplot(data = gene_exprs_df, aes(x = log2(gene_exprs), y= log2(gene_exprs_mmr), label = labs))+
  geom_point(size =0.3)+
  annotate("text", x = 3, y =-2, label = paste0("Pearson's r = ", round(cor$estimate,2)), size = 2.5)+
  labs(x =  expression('Log'[2]*' mean CosMx probe count'), y = expression('Log'[2]*' mean 10x read count'))+
  blank_theme+
  theme(aspect.ratio = 1)

figure_s_0_b

# Plot the expression of the negative probes
negprb_plt <- merged@assays$RNA$counts[grepl("Neg", rownames(merged)),]
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
negprb <- rownames(merged)[grepl("Neg", rownames(merged))]
ctrl <- rownames(merged)[grepl("SystemControl", rownames(merged))]
merged <- merged@assays$RNA$counts
# https://kb.10xgenomics.com/hc/en-us/articles/360004729092-Why-do-I-see-high-levels-of-Malat1-in-my-gene-expression-data-
# 10x says dead and dying cells might have high MALAT1
# These genes do not help with celltype ID
merged <- merged[!rownames(merged) %in% c(negprb, ctrl, "NEAT1", "MALAT1"),]

# Make a seurat object, only keeping genes detected in > 50 cells
merged <- CreateSeuratObject(merged, min.cells = 50)
dim(merged)
tail(colnames(merged))
merged$Barcode <- colnames(merged)
unique(merged$orig.ident)
# Join the datasets on 'Barcode'
merged@meta.data <- left_join(merged@meta.data, md_combined)
# Reset the rownames
rownames(merged@meta.data) <- merged@meta.data$Barcode

# Save the object containing all cells
qsave(merged, "./Intermediate/All_cells_g.qs")

# merged <- qread("./Intermediate/All_cells_g.qs")

# Make a higher quality object to 
# Get the total number of cells
total_start <- length(merged$Barcode)

# Loop over each sample and calculate a silhouette score
files <- list.files("./Results_paper/analyses/single_sample_analyses/", pattern = "*.qs", recursive = T, full.names = T)

# A list for sil QC stats
sil_list <- list()
# Make a list of the cells that passed QC
good_cells_list <- list()
for (i in 1:length(files)){
  
  print(i)
  
  file <- files[i]
  
  sample_name <- gsub("_seurat_processed.qs", "", basename(file))
  
  sample <- qread(file)
  
  sil_df <- data.frame(Sample = sample_name, median_sil = median(sample$sil), 
                       medcount = median(sample$nCount_RNA),
                       clusters = max(as.numeric(sample$seurat_clusters)))
  
  sil_list[[i]] <- sil_df
  
  good_cells_list[[i]] <- sample@meta.data
  
}

# Read in all the clinical data
clinical <- read_csv("./raw/Patients for spatial clinical info.csv")%>%
  filter(!is.na(`ORG #`))%>%
  select(Sample = `ORG #`, `Age at Dx`, `Primary tumour location`,
         Sex, Sidedness, `MMR status`, `Surgery date`)%>%
  mutate(Sample = gsub("ORG", "", Sample))

all_sil <- bind_rows(sil_list)%>%
  mutate(Slide_sample = Sample)%>%
  mutate(Sample = gsub(".*_", "", Slide_sample))%>%
  left_join(clinical)%>%
  mutate(Slide = gsub("_[^_]+$", "", Slide_sample))

sil_plot <- all_sil%>%
  filter(!is.na(`Surgery date`))%>%
  mutate(`Surgery date` = as.Date(`Surgery date`,format = "%d/%m/%Y"))

date_num <- as.numeric(sil_plot$`Surgery date`)
date_cor <- cor.test(date_num, sil_plot$median_sil)

r2_date <- date_cor$estimate^2

good_samples <- all_sil%>%
  filter(median_sil > -0.025)

figure_s_temp <- ggplot(data = sil_plot, aes(x = `Surgery date`, y = median_sil, label = Sample, colour = Slide))+
  geom_point()+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  geom_text_repel()

figure_s_temp

count_cor <- cor.test(sil_plot$medcount, sil_plot$median_sil)
count_cor

r2_count <- count_cor$estimate^2
r2_count

figure_s_temp_2 <- ggplot(data = sil_plot, aes(x = medcount, y = median_sil, label = Sample, colour = Slide))+
  geom_point()+
  geom_text_repel()+
  geom_hline(yintercept = -0.025, linetype = 2)+
  blank_theme

figure_s_temp_2

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

# Drop any cells with low counts 
merged <- merged[, merged$nCount_RNA >=50]
dim(merged)
rownames(merged)

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

# Use scDblFinder on the data
# The tool by defualt expects a 10x like doublet rate. 
# I'm not sure this will be the case for nanostring but there might 
# be cells with overlapping cytoplasm, so worth trying to clean up

# Run without cluster info
# dbs <- scDblFinder(merged@assays$RNA$counts, 
#                    samples = merged$Sample_Slide, 
#                    BPPARAM=MulticoreParam(20))
# 
# # Look at the doublet rate
# table(dbs$scDblFinder.class)
# 
# # Drop doublets
# merged$scDblFinder.class <- dbs$scDblFinder.class
# merged <- merged[,merged$scDblFinder.class == "singlet"]

# Track the total after cell count filtering
#total_doublet_cutoff <- length(merged$Barcode)

qc_df <- data.frame(`Starting cells` = total_start,
                    `Cells in good samples` = total_good_samples,
                    `Cells passing RobZ` = total_robz_cutoff,
                    check.names = F)%>%
  write_csv("./Results_paper/Tables/QC_filtering_g.csv")

# Make a plot of QC filtering steps
qc_filt <- qc_df%>%
  gather(Step, Cells)%>%
  mutate(Step = factor(Step, levels = Step))

qc1_b <- ggplot(data  = qc_filt, aes(x = Step, y = Cells))+
  geom_bar(stat = "identity")+
  blank_theme+
  labs(x = "Filtering step", y = "Total cells")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

qc1_b

# Make a QC figure for the supp 
qc1_c <- VlnPlot(merged, features = c("nCount_RNA"),
        group.by = "Sample_Slide", pt.size = -1)+
  blank_theme+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))+
  labs(x = "Sample + slide", y = "Total probes", title = NULL)

# Set narrower line widths
qc1_c$layers[[1]]$aes_params$size = 0.25

qc1_d <- VlnPlot(merged, features = c("nFeature_RNA"),
        group.by = "Sample_Slide", pt.size = -1)+
  blank_theme+
  NoLegend()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4))+
  labs(x = "Sample + slide", y = "Unique probes (genes detected)",title = NULL)

qc1_d$layers[[1]]$aes_params$size = 0.25

# Log normalsie
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

qc1_e

top_row <- plot_grid(qc1_a, qc1_b, labels = c("a", "b"), label_size = 8, rel_widths = c(2.3,1), ncol = 2)
bottom_row <- plot_grid(qc1_c,qc1_d,qc1_e, labels = c("c", "d", "e"), 
                        label_size = 8, align = "h", axis = "bt", nrow = 1,
                        rel_widths = c(1.2, 1.2, 1))

Figure_s1 <- plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1.3,1))
ggsave(plot = Figure_s1,"./Results_paper/Plots_G/Figure S1.png", width = 170, 
       height = 170, units = "mm", dpi = 600)
ggsave(plot = Figure_s1,"./Results_paper/Plots_G/Figure S1.pdf", width = 170, 
       height = 170, units = "mm")

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
DimPlot(merged, reduction = "umap", label = T, group.by = "Sample")
DimPlot(merged, reduction = "umap", label = T, group.by = "Slide")

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

qc2_d <- DimPlot(merged, reduction="umap_harmony", group.by="Slide_short")+
  blank_theme+
  theme()+
  labs(colour = "Slide", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc2_e <- DimPlot(merged, reduction="umap_harmony", group.by="Sample", label = F, label.size = 2)+
  blank_theme+
  theme(
    legend.key.size = unit(4, 'mm'))+
  labs(colour = "Donor", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc2_f <- DimPlot(merged, reduction="umap_harmony", group.by="Sample_type", label = F, label.size = 2)+
  blank_theme+
  labs(colour = "Sample type", title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

top_row <- plot_grid(qc2_a, qc2_b, labels = c("a", "b"), label_size = 8, rel_widths = c(1,1.2), ncol = 2,align = "h", axis = "bt")
middle_row <- plot_grid( qc2_c, qc2_d, labels = c("c", "d"), label_size = 8, rel_widths = c(1,1.2), align = "h", axis = "bt", nrow = 1)
bottom_row <- plot_grid( qc2_e,qc2_f, labels = c("e", "f"), label_size = 8, rel_widths = c(1,1), align = "h", axis = "bt", nrow = 1)

# Save the figure
Figure_s2 <- plot_grid(top_row, middle_row, bottom_row, nrow = 3, rel_heights = c(1,1,1),align = "hv", axis = "bt")
ggsave(plot = Figure_s2,"./Results_paper/Plots_G/Figure S2.pdf", width = 170, height = 180, units = "mm")

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
predictions <- SingleR(test=merged@assays$RNA$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                       num.threads = 30, aggr.ref = T)

merged$MMR_pred_midlevel <- predictions$labels 

# Set the annotation to be just the MMR atlas preds
merged$Manual_toplevel_pred <- merged$MMR_pred_midlevel

# Run the MMR atlas at the toplevel
predictions <- SingleR(test=merged@assays$RNA$counts, 
                       ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clTopLevel,
                       num.threads = 30, aggr.ref = T)

merged$MMR_pred_toplevel <- predictions$labels

# Run the MMR atlas at the toplevel
predictions <- SingleR(test=merged@assays$RNA$data, 
                       ref=mmr_atlas@assays$RNA@data, labels=mmr_atlas$cl295v11SubFull,
                       num.threads = 30, aggr.ref = T)

merged$MMR_pred_lowlevel <- predictions$labels

# Intermediate save
#qsave(merged, "./Intermediate/lognorm_merged_integrated_annotated_g.qs")

#merged <- qread("./Intermediate/lognorm_merged_integrated_annotated_g.qs")

qc_3_a <- DimPlot(merged, reduction="umap_harmony", group.by="PrimaryCellAtlas_pred_main", label = T, label.size = 2, repel = T)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "Main primary cell atllas pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc_3_b <- DimPlot(merged, reduction="umap_harmony", group.by="MMR_pred_toplevel", label = T, label.size = 2)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "MMR atlas toplevel pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

qc_3_c <- DimPlot(merged, reduction="umap_harmony", group.by="MMR_pred_midlevel", label = T, label.size = 2, cols = cell_type_colors)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "MMR atlas midlevel pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

# I don't think resolution is there for the lowlevel preds
qc_3_d <- DimPlot(merged, reduction="umap_harmony", group.by="MMR_pred_lowlevel", label = T, label.size = 2, repel = T)+
  blank_theme+
  labs(colour = "Harmony\ncluster", title = "MMR atlas lowlevel pred", x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")+
  NoLegend()

# Save the figure
Figure_s3 <- plot_grid(qc_3_a, qc_3_b, qc_3_c, qc_3_d, labels = c("a", "b", "c", "d"), label_size = 8, align = "h", axis = "bt", nrow = 2)
ggsave(plot = Figure_s3,"./Results_paper/Plots_G/Figure S3.pdf", width = 170, height = 170, units = "mm")

# Get the celltype counts
ct_counts <- table(merged$MMR_pred_midlevel)%>%
  data.frame()%>%
  arrange(-Freq)

# Take a look at some clusters in fine deatil
# Try with the fibroblasts
#mac_neu <- merged[, merged$seurat_clusters %in% c(5, 11)]
# mac_neu <- merged[, merged$MMR_pred_midlevel %in% c("Mono", "Macro", "Granulo", "DC")]
# unique(merged$PrimaryCellAtlas_pred_main)

# Take a look at the object
# mac_neu
# # Run UMAP clustering on the dataset
# mac_neu <- RunPCA(mac_neu, verbose = FALSE)
# ElbowPlot(mac_neu, ndims = 30)+
#   blank_theme
# mac_neu <- RunHarmony(mac_neu, c("Sample", "Slide"), reduction="pca", reduction.save="harmony")
# mac_neu <- FindNeighbors(mac_neu, reduction="harmony", dims=1:15)
# mac_neu <- FindClusters(mac_neu, cluster.name = "harmony_clusters")
# mac_neu <- RunUMAP(mac_neu, reduction="harmony", dims=1:15, reduction.name="umap_harmony")

# Cluster 6 looks pretty strongly like monocytes here
# DimPlot(mac_neu,label = T,reduction = "umap_harmony")
# c5 <- FindMarkers(mac_neu, ident.1 = 5, ident.2 = 1)
# DimPlot(mac_neu,label = T, group.by = "MMR_pred_midlevel",reduction = "umap_harmony")
# DimPlot(mac_neu,label = T, group.by = "Sample",reduction = "umap_harmony")
# DimPlot(mac_neu,label = T, group.by = "Sample_type",reduction = "umap_harmony")
# DimPlot(mac_neu,label = T, group.by = "Donor",reduction = "umap_harmony")
# DimPlot(mac_neu,label = T, group.by = "PrimaryCellAtlas_pred_main",reduction = "umap_harmony")+
#   NoLegend()
# FeaturePlot(mac_neu, features = c("S100A8") ,reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("S100A9") ,reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("LYZ") ,reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("HCAR2"),reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("CD3E"),reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("MS4A4A"),reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("IGHG1"),reduction = "umap_harmony")
# FeaturePlot(mac_neu, features = c("CXCL8"),reduction = "umap_harmony")

# Try and look at the immune cells
#immune <- merged[, merged$MMR_pred_toplevel %in% c("TNKILC")]
# immune <- merged[, merged$seurat_clusters %in% c("6")]
# unique(merged$PrimaryCellAtlas_pred_main)
# 
# # Take a look at the object
# immune
# # Run UMAP clustering on the dataset
# immune <- RunPCA(immune, verbose = FALSE)
# ElbowPlot(mac_neu, ndims = 30)+
#   blank_theme
# immune <- RunHarmony(immune, c("Sample", "Slide"), reduction="pca", reduction.save="harmony")
# immune <- FindNeighbors(immune, reduction="harmony", dims=1:15)
# immune <- FindClusters(immune, cluster.name = "harmony_clusters")
# immune <- RunUMAP(immune, reduction="harmony", dims=1:15, reduction.name="umap_harmony")

# DimPlot(immune, label = T,reduction = "umap_harmony")
# DimPlot(immune, label = T,reduction = "umap_harmony", group.by = "MMR_pred_midlevel")
# FeaturePlot(immune, features = c("CD3E"),reduction = "umap_harmony")
# FeaturePlot(immune, features = c("CD8A"),reduction = "umap_harmony")
# FeaturePlot(immune, features = c("CD4"),reduction = "umap_harmony")

# Plot some marker genes
qc4_a <- FeaturePlot(merged, features = "EPCAM", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("EPCAM")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_b <- FeaturePlot(merged, features = "FN1", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("FN1")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_c <- FeaturePlot(merged, features = "IGHM", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("IGHM")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_d <- FeaturePlot(merged, features = "CD3E", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("CD3E")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_e <- FeaturePlot(merged, features = "SPP1", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("SPP1")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_f <- FeaturePlot(merged, features = "CD4", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("CD4")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_g <- FeaturePlot(merged, features = "ACTA2", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("ACTA2")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_h <- FeaturePlot(merged, features = "IGHG1", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("IGHG1")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

qc4_i <- FeaturePlot(merged, features = "VEGFA", max.cutoff = "q95", reduction = "umap_harmony")+
  blank_theme+
  theme(
    legend.key.size = unit(3, 'mm'))+
  labs(colour = expression(italic("VEGFA")), title = NULL, x = "UMAP\n harmony 1", y = "UMAP\n harmony 2")

Figure_s4 <- plot_grid( qc4_a, qc4_b, qc4_c, qc4_d, qc4_e, qc4_f, qc4_g, qc4_h,qc4_i, labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i"), label_size = 8, 
                        align = "hv", axis = "bt", nrow = 3)
ggsave(plot = Figure_s4,"./Results_paper/Plots_G/Figure S4.pdf", width = 170, height = 130, units = "mm")
ggsave(plot = Figure_s4,"./Results_paper/Plots_G/Figure S4.png", width = 170, height = 130, units = "mm", dpi = 600)

#merged <- qread("./Intermediate/lognorm_merged_integrated_annotated_g.qs")
rownames(merged@meta.data) <- merged$Barcode
Figure_1_a <- DimPlot(merged, reduction="umap_harmony", group.by="Manual_toplevel_pred", label = T, 
                      label.size = 2, repel = T, cols = cell_type_colors)+
  blank_theme+
  theme(legend.key.size = unit(2, 'mm'))+
  guides(color = guide_legend(override.aes = list(size = 0.5))) 

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

Figure_1_c_4 <- FeaturePlot(merged, features = c("Mean.CD68"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "CD68", title = NULL, x = NULL, y = NULL)

Figure_1_c_5 <- FeaturePlot(merged, features = c("Mean.Membrane"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "Membrane", title = NULL, x = NULL, y = NULL)

Figure_1_c_6 <- FeaturePlot(merged, features = c("Area.um2"), reduction = "umap_harmony",max.cutoff = "q95")+
  blank_theme+
  theme(
    legend.key.size = unit(2, 'mm'),
    legend.margin = margin(c(-1, -1, -1, -1)),
    plot.margin = unit(c(-1, -1, -1, -1), "cm"))+
  labs(colour = "Area.um2", title = NULL, x = NULL, y = NULL)

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
  #guides(fill = "none")+
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

plots <- plot_summary(sc_result) 

# Seems to be big differences in B cells and neutrophils
plots$boxplot

Figure_1_e <- plots$boxplot[[1]]
Figure_1_e <- Figure_1_e+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4))


# Read in some image PDFs
#Tumour_composite <- "~/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Tumour.jpg"
# normal_composite <- "~/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Normal.jpg"
# tumour_cells_f <- "~/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Figure_1_cells_tumour.png"
# normal_cells_f <- "~/Data/Spatial/Results_paper/Extra_plots/Image_inputs/Figure_1_cells_normal.png"

# tumour_comp <- ggdraw() + 
#   draw_image(
#     Tumour_composite, scale = 1, width = 1.5,halign = 0, valign = 0)+
#   blank_theme+
#   theme(panel.border = element_blank())
# 
# # Run GC to clear the rsession image cache
# gc()
# normal_comp <- ggdraw() + 
#   draw_image(
#     normal_composite,scale = 1, width = 1.5, halign = 0, valign = 0)+
#   blank_theme+
#   theme(panel.border = element_blank())
# 
# normal_cells <- ggdraw() + 
#   draw_image(
#     normal_cells_f,scale = 1,width = 1.5,halign = 0.2, valign = 0)+
#   blank_theme+
#   theme(panel.border = element_blank())
# 
# tumour_cells <- ggdraw() + 
#   draw_image(
#     tumour_cells_f,scale = 1,width = 1.5,halign = 0.2, valign = 0)+
#   blank_theme+
#   theme(panel.border = element_blank())
# 
left <- plot_grid(Figure_1_a,Figure_1_b, labels = c("a", "b"), label_size = 8,
                  rel_heights = c(1.6,1), nrow = 2)
right <- plot_grid(Figure_1_c_1, Figure_1_c_2, Figure_1_c_3, Figure_1_c_4, Figure_1_c_5, Figure_1_c_6, labels = c("c"), label_size = 8,
                   ncol = 2, scale = 1, align = "hv",axis = "bt")
top <- plot_grid(left, right, rel_widths = c(1,1.1))
#bottom_left <- plot_grid(normal_cells, normal_comp, tumour_cells, tumour_comp, nrow = 2, labels = "c", label_size = 8)
bottom_right <- plot_grid(Figure_1_e, nrow = 1, labels = c("e"),label_size = 8)
bottom <- plot_grid(bottom_right)
Figure_1 <- plot_grid(top, bottom, nrow = 2, rel_heights = c(1.34,1))
ggsave(plot = Figure_1,"./Results_paper/Plots_G/Figure 1.png", width = 170,
       height = 170, units = "mm", dpi = 600)

# Save the final object
qsave(merged, "./Intermediate/lognorm_merged_integrated_annotated_g.qs")

# Run singleR to annotate the cells and then group them into niches
# Load in the nanostring data again
g1 <- LoadNanostring.FIX (data.dir = "./raw/10501A/",fov = "fov_all")
g2 <- LoadNanostring.FIX (data.dir = "./raw/120131/",fov = "fov_all")
g3 <- LoadNanostring.FIX (data.dir = "./raw/148191A/",fov = "fov_all")
g4 <- LoadNanostring.FIX (data.dir = "./raw/3461G/",fov = "fov_all")

# Grab the metadata from the integrated object
md_all <- merged@meta.data

# Check the slide names for the next steps
unique(md_all$Slide)

slidelist <- list(g1, g2, g3, g4)
slidenames <- c("10501A_", "120131_", "148191A_", "3461G_")

slidelist_knn <- list()
for(i in 1:length(slidelist)){
  
  slide <- slidelist[[i]]
  
  slide_name <- slidenames[i]
  
  # Grab the knn matrix from each slide
  slidecounts <- grab_nichecounts(slide, slide_name, md_all, k_param = 20)
  
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
qsave(combined_niche, "./Intermediate/Combined_niche_knn_g.qs")
#combined_niche <- qread("./Intermediate/Combined_niche_knn_g.qs")

combined_niche[1:5,1:5]

sampled <- sample(ncol(combined_niche), 10000)

combined_niche_samp <- combined_niche[,sampled]

set.seed(42)

# Try flowSOM to see if I can see neighbourhoods
# fSOM <- FlowSOM(t(combined_niche),
#                   # Input options:
#                   compensate = FALSE,
#                   transform = FALSE,
#                   scale = F,
#                   # Metaclustering options:
#                   nClus = 8,
#                   seed = 42)
# 
# # Plot a summary of the data
# FlowSOMmary(fSOM, plotFile = "./Results_paper/Extra_plots/FlowSOMmary.pdf")

# Each cell commes out mapped to its 20 nearest neighbors (in terms of actual location)
combined_niche[1:5,1:5]
colSums(combined_niche[,1:5])

combined_niche_seu <- CreateSeuratObject(combined_niche)
combined_niche_seu <- NormalizeData(combined_niche_seu)

combined_niche_seu <- FindVariableFeatures(combined_niche_seu)

top10 <- head(VariableFeatures(combined_niche_seu), 19)

plot1 <- VariableFeaturePlot(combined_niche_seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot2
combined_niche_seu <- ScaleData(combined_niche_seu)
combined_niche_seu <- RunPCA(combined_niche_seu, features = VariableFeatures(object = combined_niche_seu))
DimPlot(combined_niche_seu, reduction = "pca") + NoLegend()

merged_md<- merged@meta.data

md_all_niche <- combined_niche_seu@meta.data%>%
  rownames_to_column("Barcode")%>%
  select(-nCount_RNA, -nFeature_RNA)%>%
  left_join(merged_md)

rownames(md_all_niche) <- md_all_niche$Barcode
combined_niche_seu@meta.data <- md_all_niche

DimPlot(combined_niche_seu, reduction = "pca", group.by = "Sample", label =T)

DimPlot(combined_niche_seu, reduction = "pca", group.by = "Slide", label =T)

# Is there a slide level batch effect in cell types that are annotated in the data?
DimPlot(combined_niche_seu, reduction = "pca", group.by = "Slide", label =T) + NoLegend()
# There is a clear 
Figure_2_a <- DimPlot(combined_niche_seu, reduction = "pca", group.by = "Sample_type", label =T) + 
  blank_theme+
  labs(title = "")+
  NoLegend()

Figure_2_a

# Scale the data before KNN
scaled <- ScaleData(combined_niche)
# Sample the object for the knn estimation
set.seed(42)

# Sample about 10% of cells
sampled <- sample(ncol(scaled), 20000)
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

# Set the seed again
set.seed(42)

# Kmeans is fast
# Set 20 iterations to hopfeully converge
results <- kmeans(x = t(scaled), centers = niches.k, iter.max = 20,
                  nstart = 30)

scaled[1:5,1:5]

#results$centers
niche_df <- data.frame(Barcode = colnames(scaled), niches = results$cluster)

md_all <- merged@meta.data

# Map the niches back onto the merged integrated object
# Rename the niches and match them to the metadata
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
  mutate(niches = replace(niches, niches == 1, "Macro mono rich"))%>%
  mutate(niches = replace(niches, niches == 2, "TLS"))%>%
  mutate(niches = replace(niches, niches == 3, "Normal like"))%>%
  mutate(niches = replace(niches, niches == 4, "Stroma"))%>%
  mutate(niches = replace(niches, niches == 5, "Border"))%>%
  mutate(niches = replace(niches, niches == 6, "Vasculature"))%>%
  mutate(niches = replace(niches, niches == 7, "Epithelial mass"))%>%
  mutate(niches = replace(niches, niches == 8, "Granulocyte rich"))%>%
  mutate(niches = replace(niches, niches == 9, "Normal immune infil"))

merged$niches <- niche_match$niches

unique(merged$niches)

nichetab <- merged@meta.data%>%
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
#all_cells <- qread("./Intermediate/All_cells.qs")

# Create a cell/niche label 
merged$cell_niche <- paste0(merged$Manual_toplevel_pred, "_", merged$niches)

low_qual <- all_cells[,!all_cells$Barcode %in% merged$Barcode]

# Run singleR for the low quailty cells to give them an annotaion to visualise
predictions <- SingleR(test=low_qual@assays$RNA$counts,aggr.ref = T,
                       ref=merged@assays$RNA$counts, labels=merged$cell_niche, 
                       num.threads = 30)

# Add a prediction to the low quailty cells
low_qual$cell_niche <- predictions$labels

# Read the annotation of the high quailty cells
hq_labs <- merged@meta.data%>%
  select(Barcode, cell_niche)

lq_labs <- low_qual@meta.data%>%
  select(Barcode,cell_niche_match= cell_niche)

anno_matched <- all_cells@meta.data%>%
  left_join(hq_labs)%>%
  left_join(lq_labs)%>%
  mutate(cell_niche = replace(cell_niche, is.na(cell_niche), cell_niche_match[is.na(cell_niche)]))%>%
  mutate(Manual_toplevel_pred = gsub("_.*", "", cell_niche))%>%
  # Non greedy match
  mutate(niches = gsub("^.*?_", "", cell_niche))

all_cells$Manual_toplevel_pred <- anno_matched$Manual_toplevel_pred
all_cells$niches <- anno_matched$niches

# All matched successfully 
sum(is.na(all_cells$niches))

slide_short <- data.frame(Slide= all_cells$Slide)%>%
  left_join(slide_alias)

all_cells$Slide_short <- slide_short$slide_names
all_cells$Sample_Slide <- paste0(all_cells$Sample, " ", all_cells$Slide_short)

# Add on T/N and Donor annotation 
all_cells$Sample_type <- ifelse(grepl("T", all_cells$Sample), "T", "N")
all_cells$Donor <- gsub("T|N","", all_cells$Sample)
all_cells$Donor <- paste0("Donor", all_cells$Donor)

# I think a min cutoff of 20 mols is at least required to have some
# idea of cell type
all_cells$Manual_toplevel_pred <- replace(all_cells$Manual_toplevel_pred, all_cells$nCount_RNA <20, "QC_fail")
all_cells$niches <- replace(all_cells$niches, all_cells$nCount_RNA <20, "QC_fail")

donor_md <- read_csv("./raw/Donor slide FOV info donor level.csv")

md_all <- all_cells@meta.data%>%
  left_join(donor_md)

rownames(md_all) <-md_all$Barcode
all_cells@meta.data <- md_all

# Save the annotated all cells object
qsave(all_cells, "./Intermediate/All_cells_annotated_g.qs")

nichetab <- merged@meta.data%>%
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

ggsave(plot = Figure_s5,"./Results_paper/Plots_G/Figure S5.pdf", width = 170, height = 180, units = "mm")

# Plot the niches per FOV
nichetab <- merged@meta.data%>%
  mutate(fov_numeric = fov)%>%
  mutate(fov = paste0(fov, " ", Sample))%>%
  group_by(Slide, fov)%>%
  mutate(Total = n())%>%
  mutate(fov = as.character(fov))%>%
  mutate(niches = as.character(niches))%>%
  group_by(fov, fov_numeric, niches, Total, Slide)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Niche_pct = Count/Total *100)%>%
  arrange(fov_numeric)%>%
  mutate(fov = factor(fov, levels = unique(fov)))

plt <- ggplot(data = nichetab, aes(x = fov, y = Niche_pct, fill = niches))+
  geom_bar(stat = "identity")+
  facet_wrap(~Slide, scales  = "free", ncol = 3)+
  scale_fill_manual(values = nichecols)+
  blank_theme+
  labs(x = "FOV", y = "FOV niche %", fill = "Niche")+
  ggtitle("FOV/Slide")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plt

save <- "./Results_paper/Extra_plots/FOV_niches.pdf"
ggsave(filename = save,plot = plt, width = 20, height = 8) 

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
ggsave(filename = save,plot = plt, width = 9, height = 30) 

# Plot the per-sample niche composition
nichetab_per_sample <- merged@meta.data%>%
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
Tumour_cells <- "./Results_paper/Extra_plots/Figure_2_image_inputs/Tumour_celltype.png"
Tumour_niche <- "./Results_paper/Extra_plots/Figure_2_image_inputs/Tumour_niche.png"
Normal_cells <- "./Results_paper/Extra_plots/Figure_2_image_inputs/Normal_celltype.png"
Normal_niche <- "./Results_paper/Extra_plots/Figure_2_image_inputs/Normal_niches.png"

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
sc_result <- merged |>
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
top_right <- plot_grid(Figure_2_b, Figure_2_c, labels = c("b", "c"), label_size = 8, nrow = 2)
top_row <- plot_grid(top_left, top_right, rel_widths = c(1,1))
bottom_left <- plot_grid(tumour_c, tumour_n, Normal_c,Normal_n, labels = c("d"), label_size = 8, align = "h", axis = "bt", nrow = 2)
bottom_right <- plot_grid(Figure_2_e, labels = c("e"), label_size = 8)
bottom_row <- plot_grid(bottom_left, bottom_right, rel_widths = c(1,1))
Figure_2 <- plot_grid(top_row,bottom_row, nrow = 2, rel_heights = c(1,1.15))
ggsave(plot = Figure_2,"./Results_paper/Plots_G/Figure 2.png", width = 170, 
       height = 140, units = "mm", dpi = 600)

# Plot cell size for a supp figure
size_per_sample <- merged@meta.data%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Medain_area_cell = median(Area))%>%
  ungroup()%>%
  arrange(-Medain_area_cell)%>%
  mutate(Slide_short = factor(Slide_short, levels = unique(Slide_short)))%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))

Figure_s6_a <- ggplot(data = size_per_sample, aes(x = Slide_short, y = Medain_area_cell, colour = Manual_toplevel_pred, group = Manual_toplevel_pred))+
  geom_point()+
  geom_line()+
  labs( y = "Median pixel area per sample", colour = "Cell type", x = "Slide")+
  blank_theme+
  scale_colour_manual(values = cell_type_colors)

Figure_s6_a

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

clustmarks <- FindAllMarkers(merged)%>%
  rownames_to_column("Gene")

write_csv(clustmarks, "./Results_paper/Tables/Object cell type markers.csv")

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
ggsave(plot = Figure_2,"./Results_paper/Plots_G/Figure 2.png", width = 170, 
       height = 140, units = "mm", dpi = 600)


# Save the final object
qsave(merged, "./Intermediate/lognorm_merged_integrated_annotated_g.qs")

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

# Make a plot of cell size and transcript counts
size_vs_count <- merged@meta.data%>%
  group_by(Manual_toplevel_pred, Slide_short)%>%
  summarise(Median_area_cell = median(Area), Median_transcript_counts =  median(nCount_RNA))%>%
  ungroup()

size_cor <- cor.test(size_vs_count$Median_area_cell, size_vs_count$Median_transcript_counts)

size_cor <-round(size_cor$estimate,digits = 2)

size_cor_fig <- ggplot(data = size_vs_count, aes(x = Median_area_cell, y = Median_transcript_counts, 
                                   colour = Manual_toplevel_pred, shape = Slide_short))+
  geom_point()+
  scale_colour_manual(values = cell_type_colors)+
  blank_theme+
  labs(x = "Median cell area (pixels)", 
       y = "Median probe detection per cell",
       shape = "Slide",
       colour = "Cell type")+
  guides(shape = "none")+
  theme(legend.key.size = unit(2.5, 'mm'))+
  annotate("text", x = 5000, y =350, label = paste0("Pearson's r = ", size_cor), size = 2.5)

# Make a plot of cell size and transcript counts
size_per_cell <- merged@meta.data%>%
  group_by(Manual_toplevel_pred, Sample, Slide_short)%>%
  summarise(Median_transcript_counts =  median(nCount_RNA))%>%
  group_by(Manual_toplevel_pred)%>%
  mutate(med_med = median(Median_transcript_counts))%>%
  ungroup()%>%
  arrange(-med_med)%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))

ct_size_fig <- ggplot(data = size_per_cell, aes(x = Manual_toplevel_pred, y = Median_transcript_counts))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(height = 0, width = 0.2, aes(colour = Sample, shape  = Slide_short))+
  scale_colour_manual(values = Sample_cols)+
  blank_theme+
  labs(x = "Cell type", 
       y = "Median probe detection per cell",
       colour = "Sample",
       shape = "Slide")+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Make a plot of cell size vs count
count_vs_size <- merged@meta.data%>%
  group_by(Manual_toplevel_pred)%>%
  summarise(Median_area_cell = median(Area), Total_cells_called = n())%>%
  arrange(Total_cells_called)%>%
  mutate(label = replace(Manual_toplevel_pred,6:n(), NA))

area_cells_cor <- cor.test(count_vs_size$Median_area_cell, 
                           count_vs_size$Total_cells_called)

area_cells_cor <- round(area_cells_cor$estimate,2)

count_vs_size_fig <- ggplot(data = count_vs_size, 
                            aes(x = Median_area_cell, y = log2(Total_cells_called), 
                                colour = Manual_toplevel_pred,
                                label = Manual_toplevel_pred))+
  geom_point()+
  scale_colour_manual(values = cell_type_colors)+
  blank_theme+
  labs(x = "Median cell area (pixels)", y = expression('Log'[2]*' total cells in dataset'),
       colour = "Cell type")+
  theme(legend.key.size = unit(2, 'mm'))+
  guides(colour = "none")+
  geom_text_repel(size = 2,max.overlaps = 20)+
 annotate("text", x = 4000, y =16, label = paste0("Pearson's r = ", area_cells_cor), size = 2.5)

top <- plot_grid(ct_size_fig, labels = c("a"), label_size = 8)
bottom <- plot_grid(size_cor_fig, count_vs_size_fig, labels = c("b", "c"), 
                    label_size = 8, ncol = 2, rel_widths = c(1.3,1))
figure_s_size <- plot_grid(top, bottom, nrow = 2)

ggsave(plot = figure_s_size,"./Results_paper/Plots_G/Figure S_size.pdf", width = 170, 
       height = 140, units = "mm")

# Make a table summarising donor characteristics 
donor_tab <- merged@meta.data%>%
  select(Sample, Slide, Donor, 53:69)%>%
  group_by_all()%>%
  summarise(Retained_cell_count = n())%>%
  ungroup()%>%
  write_csv("./Results_paper/Tables/Table S2: Sample information.csv")

