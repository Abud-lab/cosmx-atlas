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
  library(cluster)
})

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3//functions.R") 

# Set a homedir for the analysis so as to not use hardcoded file paths
setwd("/homevol/apattison/Data/Spatial")

# Set output folders
plots <- "./Results_paper/Plots/"

# Load up the MMR atlas as a singleR reference
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

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

# Add clinical data to the metadata
clinical <- read_csv("./raw/Patients for spatial clinical info.csv")%>%
  filter(!is.na(`ORG #`))%>%
  select(Sample = `ORG #`, `Age at Dx`, `Primary tumour location`, 
         Sex, Sidedness, `MMR status`, `Surgery date`)%>%
  mutate(Sample = gsub("ORG", "", Sample))

md_combined <- md_combined %>%
  left_join(clinical)

# Loop over each sample and cluster/annotate it 

# Read the counts only out of the nanostring objects
Tumor_A <- get_nano_counts(data.dir = "./raw/Run5654_Tumor_A/",sample_name = "Tumour_A")
Tumor_B <- get_nano_counts(data.dir = "./raw/Run5629_Tumor_B/",sample_name = "Tumour_B")
Normal_A <- get_nano_counts(data.dir = "./raw/Run5654_Normal_A/",sample_name = "Normal_A")
new_1 <- get_nano_counts(data.dir = "./raw/Run5850_2186_1/",sample_name = "Run5850.2186.1")
new_2 <- get_nano_counts(data.dir = "./raw/Run5850_2186_2/",sample_name = "Run5850.2186.2")
g1 <- get_nano_counts_FIX(data.dir = "./raw/10501A/",sample_name = "10501A")
g2 <- get_nano_counts_FIX(data.dir = "./raw/120131/",sample_name = "120131")
g3 <- get_nano_counts_FIX(data.dir = "./raw/148191A/",sample_name = "148191A")
g4 <- get_nano_counts_FIX(data.dir = "./raw/3461G/",sample_name = "3461G")

slide_list <- list(Tumor_A, Tumor_B, Normal_A, new_1, new_2, g1, g2, g3, g4)
rm(Tumor_A, Tumor_B, Normal_A, new_1, new_2, g1, g2, g3, g4)

outdir <- "./Results_paper/analyses/single_sample_analyses/"

for(i in 1:length(slide_list)){
  
  merged <- slide_list[[i]]
  
  merged <- merged@assays$RNA$counts
  
  negprb <- rownames(merged)[grepl("Neg|System", rownames(merged))]
  
  merged <- merged[!rownames(merged) %in% c(negprb, "NEAT1", "MALAT1"),]
  
  merged <- CreateSeuratObject(merged, min.cells = 50)
  merged$Barcode <- colnames(merged)
  merged@meta.data <- left_join(merged@meta.data, md_combined)
  rownames(merged@meta.data) <- merged@meta.data$Barcode
  
  for(sample in unique(merged$Sample)){
    
    slide <- unique(merged$Slide)
    outdir_sample <- paste0(outdir, slide, "_", sample, "/")
    system(paste0("mkdir -p ", outdir_sample))
    
    merged_s <- merged [,merged$Sample == sample]
    
    merged_s <- merged_s[, merged_s$nCount_RNA >=50]
    
    # Calculate the MAD values for counts features
    # Do this for each individual sample
    md <- merged_s@meta.data%>%
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
    merged_s@meta.data <- md
    
    min_QC_robz <- 3
    
    # Subset down based on robust Z score cutoffs for each sample
    merged_s <- subset(merged_s, subset = robzscore_nFeature_RNA < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)
    
    # Add on T/N and Donor annotation 
    merged_s$Sample_type <- ifelse(grepl("T", merged_s$Sample), "T", "N")
    merged_s$Donor <- gsub("T|N","", merged_s$Sample)
    
    # Run without cluster info
    dbs <- scDblFinder(merged_s@assays$RNA$counts, 
                       BPPARAM=MulticoreParam(20))
    
    # Look at the doublet rate
    table(dbs$scDblFinder.class)
    
    # Drop doublets
    merged_s$scDblFinder.class <- dbs$scDblFinder.class
    merged_s <- merged_s[,merged_s$scDblFinder.class == "singlet"]
    
    # Log normalsie
    merged_s <- NormalizeData(merged_s)
    merged_s <- FindVariableFeatures(merged_s, selection.method = "vst")
    
    # Identify the 20 most highly variable genes
    t20 <- head(VariableFeatures(merged_s), 20) 
    
    # Scale data and run PCA
    merged_s <- ScaleData(merged_s)
    merged_s <- RunPCA(merged_s, features = VariableFeatures(object = merged_s))
    
    print(merged_s[["pca"]], dims = 1:5, nfeatures = 5)
    
    # Find neighbour cells (in PCA space, not real space)
    merged_s <- FindNeighbors(merged_s, dims = 1:30)
    
    # Run with default cluster params
    merged_s <- FindClusters(merged_s)
    
    # Run UMAP
    merged_s <- RunUMAP(merged_s, dims = 1:30)
    
    # Add on the annotations for the MMR atlas
    predictions <- SingleR(test=merged_s@assays$RNA$counts, 
                           ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                           num.threads = 30, aggr.ref = T)
    
    merged_s$Manual_toplevel_pred <- predictions$labels 
    
    reduction <- "pca"
    dims <- 1:30
    dist.matrix <- dist(x = Embeddings(object = merged_s[[reduction]])[, dims]) 
    clusters <- merged_s$Manual_toplevel_pred
    sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
    merged_s$sil <- sil[, 3]
    
    plt <- DimPlot(merged_s, reduction = "umap", label = T, group.by = "Manual_toplevel_pred", cols = cell_type_colors)+
      ggtitle(paste0(slide, " ", sample))
    
    ggsave(plot = plt, filename = paste0(outdir_sample, slide, "_", sample,"_UMAP.pdf"), width = 10, height = 8)
    
    # Save the object
    qsave(merged_s, paste0(outdir_sample, slide, "_", sample,"_seurat_processed.qs"))
    
  }
  
}

# Try and improve clustering with less genes
good <- qread("./Results_paper/analyses/single_sample_analyses/120131_280T/120131_280T_seurat_processed.qs")
bad <- qread("./Results_paper/analyses/single_sample_analyses/3461G_137T/3461G_137T_seurat_processed.qs")
# 
# ElbowPlot(good)
# ElbowPlot(bad)
# 
# # https://github.com/satijalab/seurat/issues/1985
# # silhouette metric
# https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36
# 
# Loop over each sample and calculate a silhouette score
files <- list.files("./Results_paper/analyses/single_sample_analyses/", pattern = "*.qs", recursive = T, full.names = T)

# reduction <- "umap"
# dims <- 1:2

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
  good_cells_list[[i]] <- sample
  
}

all_sil <- bind_rows(sil_list)%>%
  mutate(Sample_slide = Sample)%>%
  mutate(Sample = gsub(".*_", "", Sample_slide))%>%
  left_join(clinical)%>%
  mutate(Slide = gsub("_[^_]+$", "", Sample_slide))

sil_plot <- all_sil%>%
  filter(!is.na(`Surgery date`))%>%
  mutate(`Surgery date` = as.Date(`Surgery date`,format = "%d/%m/%Y"))

date_num <- as.numeric(sil_plot$`Surgery date`)
cor.test(date_num, sil_plot$median_sil)

ggplot(data = sil_plot, aes(x = `Surgery date`, y = median_sil, label = Sample, colour = Slide))+
  geom_point()+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  geom_text_repel()

ggplot(data = sil_plot, aes(x = medcount, y = median_sil, label = Sample, colour = Slide))+
  geom_point()+
  geom_text_repel()