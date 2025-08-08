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
  library(slingshot)
  library(SpatialExperiment)
  library(scider)
})

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R") 

# Try clustering the highest quality sample and see how that goes
# Read in the slide FOV info
setwd("/oldvol/apattison/Data/Spatial")

slide_FOV <- read_csv("./Results_paper/Tables/Slide-FOV sample metadata.csv")

# Read in the metadata for everything
md_all <- read_csv("./Results_paper/Tables/All cells metadata.csv")

# Make each slide into a SpatialExperiment object and plot it in a simple format
read_slide_original <- function(data.dir, slide_name, md_all, outdir, fix= F){
  
  if(fix ==T){
    
    g2 <- LoadNanostring.FIX (data.dir = data.dir,fov = "fov_all")
  }else{
    g2 <- LoadNanostring (data.dir = data.dir,fov = "fov_all")
  }
  
  g2_md <- g2@meta.data%>%
    rownames_to_column("Cell")%>%
    mutate(Barcode = paste0(slide_name, Cell))%>%
    dplyr::select(Barcode, Cell)%>%
    left_join(md_all)
  
  rownames(g2_md) <- g2_md$Cell
  g2@meta.data <- g2_md
  
  # Make a plot of the full slide
  g2_plot <- g2
  
  # Set the boudary to centroids to get an XY per cell
  DefaultBoundary(g2_plot[["fov_all"]]) <- "centroids"
  
  # Get the XY coords
  fov_coords <- GetTissueCoordinates(g2_plot[["fov_all"]])
  
  fov_coords[1:2,1:2]
  
  cts <- g2_plot@assays$Nanostring$counts
  
  spe1 <- SpatialExperiment(
    assay = cts, 
    colData = fov_coords, 
    spatialCoordsNames = c("x", "y"))
  
  # Add cell type info to the object
  spe1$cell_type <- g2_plot$Manual_toplevel_pred
  spe1$FOV <- g2_plot$fov
  spe1$Sample <- g2_plot$Sample
  
  qsave(spe1, paste0(outdir, "/", slide_name, ".qs"))
  
  return(spe1)
  
}

get_plot <- function(spatial_obj, slide, rev_y =T){
  
  # Get the required info for a ggplot out of the SPE object
  toplot <- as.data.frame(SpatialExperiment::spatialCoords(spatial_obj))
  colnames(toplot) <- c("x", "y")
  cdata <- as.data.frame(SummarizedExperiment::colData(spatial_obj))
  toplot <- bind_cols(cdata, toplot)
  
  toplot$FOV <- as.character(toplot$FOV)
  
  df_text <- toplot %>% 
    mutate(FOV_sample = paste0(Sample," (", FOV, ")"))%>%
    group_by(FOV_sample) %>% 
    summarise(avg_x = mean(x),
              avg_y = mean(y),
              n = n()) %>%
    ungroup()
  
  yoff <- 0
  xoff <- 0
  
  # Plot the data
  plt <- ggplot(data= toplot, aes(x = x, y = y, colour = cell_type))+
    geom_point(size =0.1, alpha = 0.3)+
    scale_colour_manual(values = cell_type_colors)+
    labs(colour = "Cell type", title = slide)+
    guides(colour = "none")+
    geom_text(
      data = df_text, colour ="black", alpha = 1, point.padding = 0,
      nudge_x = 0,
      aes(label = FOV_sample, x = avg_x, y = avg_y)
    )
  
  if(rev_y == T){
    plt <- plt+scale_y_reverse()
  }
  
  return(plt)
}

unique(md_all$Slide)
outdir <- "/oldvol/apattison/Data/Spatial/raw_to_publish/SPE_objects/"

t_a <- read_slide_original(data.dir = "./raw/Run5654_Tumor_A/", slide_name = "Tumour_A_", 
                           md_all = md_all, outdir = outdir)
plt_TA <- get_plot(spatial_obj = t_a, slide = "Tumour A")
ggsave(plot = plt_TA, filename = paste0(outdir, "/Plots/Tumor_A.pdf"), width = 7, height = 7)
ggsave(plot = plt_TA, filename = paste0(outdir, "/Plots/Tumor_A.png"), width = 7, height = 7)

# Tumour B
t_b <- read_slide_original(data.dir = "./raw/Run5629_Tumor_B/", slide_name = "Tumour_B_", 
                           md_all = md_all, outdir = outdir)
plt_TB <- get_plot(spatial_obj = t_b, slide = "Tumour B")+
  coord_flip()
plt_TB
ggsave(plot = plt_TB, filename = paste0(outdir, "/Plots/Tumor_B.pdf"), width = 7, height = 7)
ggsave(plot = plt_TB, filename = paste0(outdir, "/Plots/Tumor_B.png"), width = 7, height = 7)

# Normal A
n_a <- read_slide_original(data.dir = "./raw/Run5654_Normal_A/", slide_name = "Normal_A_", 
                           md_all = md_all, outdir = outdir)
plt_n_a <- get_plot(spatial_obj = n_a, slide = "Normal A", rev_y = F)
ggsave(plot = plt_n_a, filename = paste0(outdir, "/Plots/Normal_A.pdf"), width = 7, height = 7)
ggsave(plot = plt_n_a, filename = paste0(outdir, "/Plots/Normal_A.png"), width = 7, height = 7)

# Slide 4
s4 <- read_slide_original(data.dir = "./raw/Run5850_2186_1/", slide_name = "Run5850.2186.1_", 
                           md_all = md_all, outdir = outdir)
plt_s4 <- get_plot(spatial_obj = s4, slide = "Run5850.2186.1", rev_y = F)
ggsave(plot = plt_s4, filename = paste0(outdir, "/Plots/Run5850.2186.1.pdf"), width = 7, height = 7)
ggsave(plot = plt_s4, filename = paste0(outdir, "/Plots/Run5850.2186.1.png"), width = 7, height = 7)
# Slide 5
s5 <- read_slide_original(data.dir = "./raw/Run5850_2186_2/", slide_name = "Run5850.2186.2_", 
                          md_all = md_all, outdir = outdir)
plt_s5 <- get_plot(spatial_obj = s5, slide = "Run5850.2186.2", rev_y = F)
ggsave(plot = plt_s5, filename = paste0(outdir, "/Plots/Run5850.2186.2.pdf"), width = 7, height = 7)
ggsave(plot = plt_s5, filename = paste0(outdir, "/Plots/Run5850.2186.2.png"), width = 7, height = 7)
# G slides
s_10501A <- read_slide_original(data.dir = "./raw/10501A/", slide_name = "10501A_", 
                          md_all = md_all, outdir = outdir, fix =T)
plt_s_10501A <- get_plot(spatial_obj = s_10501A, slide = "10501A", rev_y = F)
ggsave(plot = plt_s_10501A, filename = paste0(outdir, "/Plots/10501A.pdf"), 
       width = 7, height = 7)
ggsave(plot = plt_s_10501A, filename = paste0(outdir, "/Plots/10501A.png"), 
       width = 7, height = 7)

s_120131 <- read_slide_original(data.dir = "./raw/120131/", slide_name = "120131_", 
                                md_all = md_all, outdir = outdir, fix =T)
plt_s_120131 <- get_plot(spatial_obj = s_120131, slide = "120131", rev_y = F)
ggsave(plot = plt_s_120131, filename = paste0(outdir, "/Plots/120131.pdf"), 
       width = 7, height = 7)
ggsave(plot = plt_s_120131, filename = paste0(outdir, "/Plots/120131.png"), 
       width = 7, height = 7)

s_148191A <- read_slide_original(data.dir = "./raw/148191A/", slide_name = "148191A_", 
                                md_all = md_all, outdir = outdir, fix =T)
plt_s_148191A <- get_plot(spatial_obj = s_148191A, slide = "148191A", rev_y = F)
ggsave(plot = plt_s_148191A, filename = paste0(outdir, "/Plots/148191A.pdf"), 
       width = 7, height = 7)
ggsave(plot = plt_s_148191A, filename = paste0(outdir, "/Plots/148191A.png"), 
       width = 7, height = 7)

s_3461G <- read_slide_original(data.dir = "./raw/3461G/", slide_name = "3461G_", 
                                 md_all = md_all, outdir = outdir, fix =T)
plt_s_3461G <- get_plot(spatial_obj = s_3461G, slide = "3461G", rev_y = F)
ggsave(plot = plt_s_3461G, filename = paste0(outdir, "/Plots/3461G.pdf"), 
       width = 7, height = 7)
ggsave(plot = plt_s_3461G, filename = paste0(outdir, "/Plots/3461G.png"), 
       width = 7, height = 7)

