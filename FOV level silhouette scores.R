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
  library(FlowSOM)
  library(ggrepel)
  library(progressr)
  library(cluster)
})

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R") 

# Set a homedir for the analysis so as to not use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

# Set output folders
plots <- "./Results_paper/Plots/"

# Load up the MMR atlas as a singleR reference
mmr_atlas <- qread("/oldvol/apattison/Data/Reference/GSE178341_lognorm_annotated.qs")

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

outdir <- "./Results_paper/analyses/single_sample_analyses_fov/"

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
    
    merged_s <- merged [,merged$Sample == sample]
    
    # Get all the FOVs in the sample
    fovs <- unique(merged_s$fov)
    
    for(fov in fovs){
      
      sample_fov <- merged_s[,merged_s$fov == fov]
      
      dim(sample_fov)
      
      outdir_sample <- paste0(outdir, slide, "_", sample, "_",fov, "/")
      system(paste0("mkdir -p ", outdir_sample))
      
      sample_fov <- sample_fov[, sample_fov$nCount_RNA >=50]
      
      # Calculate the MAD values for counts features
      # Do this for each individual sample
      md <- sample_fov@meta.data%>%
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
      sample_fov@meta.data <- md
      
      min_QC_robz <- 3
      
      # Subset down based on robust Z score cutoffs for each sample
      sample_fov <- subset(sample_fov, subset = robzscore_nFeature_RNA < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)
      
      # Add on T/N and Donor annotation 
      sample_fov$Sample_type <- ifelse(grepl("T", sample_fov$Sample), "T", "N")
      sample_fov$Donor <- gsub("T|N","", sample_fov$Sample)
      
      # Run without cluster info
      dbs <- scDblFinder(sample_fov@assays$RNA$counts, 
                         BPPARAM=MulticoreParam(20))
      
      # Look at the doublet rate
      table(dbs$scDblFinder.class)
      
      # Drop doublets
      sample_fov$scDblFinder.class <- dbs$scDblFinder.class
      sample_fov <- sample_fov[,sample_fov$scDblFinder.class == "singlet"]
      
      # Log normalsie
      sample_fov <- NormalizeData(sample_fov)
      sample_fov <- FindVariableFeatures(sample_fov, selection.method = "vst")
      
      # Identify the 20 most highly variable genes
      t20 <- head(VariableFeatures(sample_fov), 20) 
      
      # Scale data and run PCA
      sample_fov <- ScaleData(sample_fov)
      sample_fov <- RunPCA(sample_fov, features = VariableFeatures(object = sample_fov))
      
      print(sample_fov[["pca"]], dims = 1:5, nfeatures = 5)
      
      # Find neighbour cells (in PCA space, not real space)
      sample_fov <- FindNeighbors(sample_fov, dims = 1:30)
      
      # Run with default cluster params
      sample_fov <- FindClusters(sample_fov)
      
      # Run UMAP
      sample_fov <- RunUMAP(sample_fov, dims = 1:30)
      
      # Add on the annotations for the MMR atlas
      predictions <- SingleR(test=sample_fov@assays$RNA$counts, 
                             ref=mmr_atlas@assays$RNA@counts, labels=mmr_atlas$clMidwayPr,
                             num.threads = 30, aggr.ref = T)
      
      sample_fov$Manual_toplevel_pred <- predictions$labels 
      
      reduction <- "pca"
      dims <- 1:30
      dist.matrix <- dist(x = Embeddings(object = sample_fov[[reduction]])[, dims]) 
      clusters <- sample_fov$Manual_toplevel_pred
      sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
      sample_fov$sil <- sil[, 3]
      
      plt <- DimPlot(sample_fov, reduction = "umap", label = T, group.by = "Manual_toplevel_pred", cols = cell_type_colors)+
        ggtitle(paste0(slide, " ", sample, " ", fov))
      
      ggsave(plot = plt, filename = paste0(outdir_sample, slide, "_", sample, "_", fov, "_UMAP.pdf"), width = 10, height = 8)
      
      # Save the object
      qsave(sample_fov, paste0(outdir_sample, slide, "_", sample, "_", fov, "_seurat_processed.qs"))
      
    }
  }
}

# Loop over the results and log sil value of each FOV
# 
# ElbowPlot(good)
# ElbowPlot(bad)
# 
# # https://github.com/satijalab/seurat/issues/1985
# # silhouette metric
# https://github.com/satijalab/Integration2019/blob/master/analysis_code/integration/integration_metrics.R#L36
# 
# Loop over each sample and calculate a silhouette score
files <- list.files("./Results_paper/analyses/single_sample_analyses_fov/", pattern = "*.qs", recursive = T, full.names = T)

# reduction <- "umap"
# dims <- 1:2

# List of silhouette values
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


clinical <- read_csv("./Results_paper/Tables/Table S1. Sample information.csv")

all_sil <- bind_rows(sil_list)%>%
  mutate(Sample_slide_fov = Sample)%>%
  mutate(Sample_slide_fov = gsub("Tumour_A", "Tumour-A", Sample_slide_fov))%>%
  mutate(Sample_slide_fov = gsub("Normal_A", "Normal-A", Sample_slide_fov))%>%
  mutate(Sample_slide_fov = gsub("Tumour_B", "Tumour-B", Sample_slide_fov))%>%
  separate(Sample_slide_fov, into = c('Slide', 'Sample', 'fov'), sep = "_")%>%
  mutate(Slide = gsub("-", "_",Slide))%>%
  left_join(slide_alias)%>%
  mutate(slide_names = paste0("Slide ",slide_names))%>%
  write_csv("./Results_paper/Tables/Per FOV silhouette scores.csv")

ggplot(data = all_sil, aes(x = Sample, y = median_sil, colour = Slide))+
  geom_point()+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cor.test(all_sil$medcount, all_sil$median_sil)

ggplot(data = all_sil, aes(x = medcount, y = median_sil, colour = Sample, label = fov))+
  geom_point()+
  facet_wrap(~slide_names, scales = "free")+
  geom_hline(yintercept = -0.025, linetype =2)+
  blank_theme+
  geom_text_repel(size =2, max.overlaps =50)+
  scale_colour_manual(values = Sample_cols)+
  labs(x = "Median probe count", y = "Median silhouette score", 
       title = "Silhouette score vs median probe count (r = 0.55)")

ggsave("./Results_paper/Plots/Figure S 00.pdf", width = 170, height = 150, units = "mm")
