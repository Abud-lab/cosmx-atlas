suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(SeuratDisk)
  library(parallel)
})

# Set a homedir for the analysis so as to not use hard coded file paths
setwd("/oldvol/apattison/Data/Spatial")

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R") 

# Make a metadata file for GEO submission (will be copy-pasted as required)
geo_meta <- read_csv("./Results_paper/Tables/Slide-FOV sample metadata.csv")

# Get a list of all the files
files <- list.files("/oldvol/apattison/Data/Spatial/raw_to_publish/", recursive = T, full.names = F, pattern = "*.csv.gz")

files_slide <- data.frame(file = files)%>%
  mutate(Slide= gsub("\\/.*","", files))%>%
  mutate(Slide = gsub("_2186_1", ".2186.1",Slide))%>%
  mutate(Slide = gsub("5850_2186_2", "5850.2186.2",Slide))%>%
  mutate(Slide = gsub("Run5.*_N", "N",Slide))%>%
  mutate(Slide = gsub("Run5.*_T", "T",Slide))%>%
  mutate(Slide = gsub("Tumor_B", "Tumour_B",Slide))%>%
  mutate(Slide = gsub("Tumor_A", "Tumour_A",Slide))%>%
  mutate(file = gsub(".*/", "", file))%>%
  mutate(Type = sub(".*?[-_]", "", file))%>%
  mutate(Type = gsub("2186_1-|2186_2-|2186_1_|2186_2_|Tumor_A_|Tumor_A-|Normal_A-|Normal_A_|Tumor_B_|Tumor_B-", "", Type))%>%
  pivot_wider(names_from = Type, values_from = file)

to_submit <- geo_meta%>%
  mutate(Sample_name = paste0(Sample, "-", Slide, "-", "FOV_", FOV))%>%
  mutate(Title = paste0(Sample, "-", Slide))%>%
  mutate(`Source name` =ifelse(grepl("T", Sample), "Colon,tumour,FFPE", "Colon,normal,FFPE"))%>%
  mutate(processed_data_file  = paste0(Sample_name, ".rds"))%>%
  mutate(Donor = replace(Donor, Donor == "DonorC", NA))%>%
  left_join(files_slide)%>%
  write_csv("/oldvol/apattison/Data/Spatial/Results_paper/Tables/GEO_metadata_table.csv")

# Read in the qs objects, name them and save them as RDS objects for GEO submission
qs_files <- list.files("/pvol/shiny_data/cached_fovs_smaller/", recursive = T, pattern = "*.qs", full.names = T)

# test <- qread("/pvol/shiny_data/cached_fovs/10501A_/FOV_1.qs")
# 
# test[['fov_all']] <- NULL
# 
# ImageDimPlot(test, fov = "zoom1", 
#              group.by = "Manual_toplevel_pred", 
#              boundaries = "segmentation",
#              #cells = cells_plot, 
#              molecules = "EPCAM",
#              mols.size = 0.2, 
#              border.color = "black",
#              dark.background = T,
#              coord.fixed = FALSE,
#              nmols = 60000)+
#   coord_flip()+
#   theme(aspect.ratio = 1)+
#   ggtitle("test")+
#   labs(fill = "Cell type")
# 
# test

FOVs <- to_submit$Sample_name

convert_qs <- function(FOv ,to_submit,qs_files){
  
  line <- to_submit%>%
    data.frame()%>%
    filter(Sample_name == FOv)
  
  line$Slide <- gsub("Tumour", "Tumor",line$Slide)
  
  to_read <- qs_files[grepl(line$Slide, qs_files)]
  to_read <- to_read[grepl(paste0("FOV_", line$FOV, ".qs"), to_read)]
  
  qfile <- qread(to_read)
  
  saveRDS(object = qfile, file = paste0("/oldvol/apattison/Data/Spatial/Results_paper/rds_for_submission/",
                               line$Sample_name, ".rds"))
  
}

mclapply(FOVs, convert_qs,to_submit,qs_files, mc.cores = 10)
