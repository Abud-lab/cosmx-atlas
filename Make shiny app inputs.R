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
  library(RColorBrewer)
  library(Polychrome)
  library(parallel)
  library(progressr)
})

options(future.globals.maxSize = 8000 * 1024^2)

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R")

# Set a homedir for the analysis so as to not use hard coded file paths
setwd("/oldvol/apattison/Data/Spatial")

outdir <- "./Results/cached_fovs/"

md_all <- read_csv("./Results_paper/Tables/All cells metadata.csv")

# Make an object with barcodes and niches to match against
niche_match <- md_all%>%
  select(Barcode, niches)

# Load up the slides
Tumor_B <- LoadNanostring(data.dir = "./raw/Run5629_Tumor_B/", fov = "fov_all")
Tumor_A <- LoadNanostring(data.dir = "./raw/Run5654_Tumor_A/", fov = "fov_all")
Normal_A <- LoadNanostring(data.dir = "./raw/Run5654_Normal_A/", fov = "fov_all")
new_1 <- LoadNanostring(data.dir = "./raw/Run5850_2186_1/",fov = "fov_all")
new_2 <- LoadNanostring(data.dir = "./raw/Run5850_2186_2/",fov = "fov_all")
g1 <- LoadNanostring.FIX (data.dir = "./raw/10501A/",fov = "fov_all")
g2 <- LoadNanostring.FIX (data.dir = "./raw/120131/",fov = "fov_all")
g3 <- LoadNanostring.FIX (data.dir = "./raw/148191A/",fov = "fov_all")
g4 <- LoadNanostring.FIX (data.dir = "./raw/3461G/",fov = "fov_all")

# Annotate each slide
Tumor_B <- annotate_obj(Tumor_B, "Tumour_B_", md_all)
Tumor_A <- annotate_obj(Tumor_A, "Tumour_A_", md_all)
Normal_A <- annotate_obj(Normal_A, "Normal_A_", md_all)
Run5850.2186.1 <- annotate_obj(new_1, "Run5850.2186.1_", md_all)
Run5850.2186.2 <- annotate_obj(new_2, "Run5850.2186.2_", md_all)
g1 <- annotate_obj(g1, "10501A_", md_all)
g2 <- annotate_obj(g2, "120131_", md_all)
g3 <- annotate_obj(g3, "148191A_", md_all)
g4 <- annotate_obj(g4, "3461G_", md_all)

slides_list <- list("Tumor_B"= Tumor_B, 
                    "Tumor_A" = Tumor_A, 
                    "Normal_A" = Normal_A, 
                    "Run5850.2186.1" = Run5850.2186.1,
                    "Run5850.2186.2" = Run5850.2186.2,
                    "10501A_" = g1,
                    "120131_" = g2,
                    "148191A_" = g3,
                    "3461G_" = g4)

rm(Tumor_B, Tumor_A, Normal_A, Run5850.2186.1, Run5850.2186.2, g1, g2, g3, g4)

# Loop through each slide and cache the FOV
for(i in 1:length(slides_list)){
  
  slide <- slides_list[[i]]
  slide_name <- names(slides_list)[i]
  
  fovs <- unique(slide$fov)
  
  # Loop through and cache each FOV
  mclapply(FUN = cache_FOV,X = fovs, mc.cores = 3, niche_match, slide, slide_name, outdir)
  
}






