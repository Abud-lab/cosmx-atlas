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
library(Polychrome)
library(ggrepel)
library(GSA)
library(parallel)
library(cowplot)

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R")

# Set a homedir for the analysis so as to no use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

merged.integrated <- qread("./Intermediate/lognorm_merged_integrated_annotated.qs")

# Make a table of donor and FOV information
donor_fov <- merged.integrated@meta.data%>%
  select(Sample, Slide, Sample_type, fov)%>%
  distinct()%>%
  write_csv("./Results_paper/Tables/Table S1: Donor slide FOV info.csv")

location_files <- list.files("/oldvol/apattison/Data/Spatial/raw/", 
                             pattern = "_tx_file.csv$", recursive = T,
                             full.names = T)

nuc_test <-  read_csv(location_files[1])

head(nuc_test)

test_image <- head(nuc_test)
  
unique(nuc_test$CellComp)
niche_matrix <- qread("./Intermediate/Combined_niche_knn_all.qs")
nm <- niche_matrix[1:19,1:5]
# # Get the counts for all slides at once
# nuclear_counts <- lapply(location_files, get_nuclear_reads)
# 
# # Order the genes the same way
# for(i in 1:length(nuclear_counts)){
#   
#   nuclear_counts[[i]] <- nuclear_counts[[i]][rownames(nuclear_counts[[1]]),]
#   
# }
# 
# # Check the rows are aligned
# sum(rownames(nuclear_counts[[1]]) == rownames(nuclear_counts[[4]]))

# # Bind the results together
# big_mat <- do.call(cbind, nuclear_counts)
# 
# # Remove everything before the first 
# colnames(big_mat) <- sub(".*?_", "", colnames(big_mat))
# head(colnames(big_mat))
# tail(colnames(big_mat))
# # Save the result
# qsave(big_mat, "./Intermediate/Nuclear_cell_table.qs")
# #big_mat <- qread("./Intermediate/Nuclear_cell_table.qs")
# dim(big_mat)

# Figure out differential niche usage between conditions
# Make a pseudobulk object ----
# Start by only looking at the tumour samples
pseudobulk_dir <- "./Results_paper/analyses/Pseudobulk_niches_HQ/"

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots/"))

# Remove spaces that break limma
merged.integrated$niches <- gsub(" ", "_", merged.integrated$niches)

to_bulk_anno <- merged.integrated@meta.data%>%
  dplyr::select(Barcode, Manual_toplevel_pred, Sample, fov, Slide, niches, fov)%>%
  mutate(Sample_type =ifelse(grepl("T", Sample), "T", "N"))%>%
  mutate(Sample_name = Sample)%>%
  mutate(Manual_toplevel_pred = gsub(" |-", "_", Manual_toplevel_pred))%>%
  # Include the FOV for this DE analysis
  mutate(Niche_sample_cell = paste0(niches,"_", Sample, "_", Manual_toplevel_pred, "_", Slide))

# Keep only conditions where we have a decent number of cells
to_drop <- to_bulk_anno %>%
  group_by(Niche_sample_cell)%>%
  summarise(count = n())%>%
  # Keep cell types with >= 20 cells per sample
  mutate(drop = count <= 20)

table(to_drop$drop)

keep <- to_drop%>%
  filter(count >= 20)

to_bulk_anno <- to_bulk_anno %>%
  filter(Niche_sample_cell %in% keep$Niche_sample_cell)

# Make sure the objects have the same cells
to_bulk_anno$Barcode[1:5]
pb_counts <- merged.integrated@assays$RNA$counts[, to_bulk_anno$Barcode] 

# Get the number of samples to test
nums <- 1:length(unique(to_bulk_anno$Niche_sample_cell))

# Parallelise making pseudobulk counts
pb_list <- mclapply(X = nums, FUN = get_pb, to_bulk_anno, aggr_col = "Niche_sample_cell", mc.cores = 15)

# Bind the PB counts
bound <- bind_cols(pb_list)

saveRDS(bound, paste0(pseudobulk_dir, "Pseudobulk_counts.rds"))

#bound <- readRDS(paste0(pseudobulk_dir, "Pseudobulk_counts.rds"))

# Convert back to a matrix
bound_counts <- bound%>%
  dplyr::select(-Gene)%>%
  as.matrix()
rownames(bound_counts) <- bound$Gene

# See how many cell types x samples we've ended up with
dim(bound_counts)

# Make a condensed annotation for the PB
condensed_SC_anno <- to_bulk_anno%>%
  dplyr::select(-Barcode)%>%
  distinct(Niche_sample_cell,Slide, .keep_all = T)

# Read in the clinical data as well
clinical <- read_csv("./Results_paper/Tables/Slide-FOV sample metadata.csv")

# Set up the annotation and keep tumour only
bulk_anno <- data.frame(Niche_sample_cell = colnames(bound_counts), check.names = F)%>%
  left_join(condensed_SC_anno)%>%
  filter(Manual_toplevel_pred != "Doublet")%>%
  # Look at tumour only
  #filter(Sample_type == "T")%>%
  mutate(Donor = gsub("T|N","", Sample_name))%>%
  mutate(Sample_type_cell = paste0(Sample_type,"_", niches, "_", Manual_toplevel_pred))%>%
  # Save the Pseudobulk annotation
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

#bulk_anno <- read_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# Make a DGElist
rnaseq <- DGEList(bound_counts[,bulk_anno$Niche_sample_cell])

# Remove weird characters
rownames(rnaseq) <- gsub("\\+", "_pos_", rownames(rnaseq))
colnames(rnaseq) <- gsub("\\+", "_pos_", colnames(rnaseq))
colnames(rnaseq) <- gsub("\\'|-|\\(|\\)", "_", colnames(rnaseq))
bulk_anno$Niche_sample_cell <- gsub("\\+", "_pos_", bulk_anno$Niche_sample_cell)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Niche_sample_cell) == length(bulk_anno$Niche_sample_cell)

design <- model.matrix(~0 + Sample_type_cell + Slide , data = bulk_anno)
# Neaten up design row and colnames
colnames(design) <- gsub("Sample_type_cell|\\+", "", colnames(design))
colnames(design) <- gsub("\\'|-|\\(|\\)", "_", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

# Get the negative probes
negprb <- rownames(rnaseq)[grepl("NegPrb", rownames(rnaseq))]
# Drop negative probes and uninformative genes
rnaseq <- rnaseq[!rownames(rnaseq) %in% c(negprb,"NEAT1", "MALAT1"),]

hist(rnaseq$counts, breaks = 100000, xlim = c(0,50))

hist(cpm(rnaseq,log = F), breaks = 100000, xlim = c(0,1000))

keep <- rowSums(cpm(rnaseq)>100)>(0.3 * nrow(design))
#keep <- filterByExpr(rnaseq, design = design, min.count = 10)
table(keep)
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]

rnaseq <- calcNormFactors(rnaseq)

# Save the counts as a reference for single cell
cpm <- cpm(rnaseq, log =T)
qsave(cpm, paste0(pseudobulk_dir, "Pseudobulk_log2_CPMs.qs"))

# Do a glimma MDS of batch removed counts
system(paste0("mkdir -p ", pseudobulk_dir,"glimma/mds/"))

# Plot the MDS on the batch removed dataset
cpm_raw <- cpm(rnaseq, log=T)

slide_removed <- removeBatchEffect(cpm_raw, batch = bulk_anno$Slide)

rownames(bulk_anno) <- bulk_anno$Niche_sample_cell

# Save and MDS plot per cell type
for(celltype in unique(bulk_anno$Manual_toplevel_pred)){
  
  print(celltype)
  
  mds_save <- paste0(paste0(pseudobulk_dir,"glimma/mds/", celltype, "_MDS.html"))
  
  bulk_filt <- bulk_anno%>%
    filter(Manual_toplevel_pred == celltype)
  
  rseq_filt <- slide_removed[,bulk_filt$Niche_sample_cell]

  if(ncol(rseq_filt) <3){
    next
  }
  
  htmlwidgets::saveWidget(glimmaMDS(rseq_filt, groups = bulk_filt,
                                    labels = bulk_filt$Sample_type_cell), mds_save)
}

v <- voom(rnaseq, design = design, plot = T)
corfit <- duplicateCorrelation(v, design, block=bulk_anno$Sample)
v <- voom(rnaseq, design, plot=TRUE,block=bulk_anno$Sample, correlation=corfit$consensus)
fit <- lmFit(v , design, block=bulk_anno$Sample, correlation=corfit$consensus)

saveRDS(fit, paste0(pseudobulk_dir, "fit_object.rds"))

#fit <- readRDS(paste0(pseudobulk_dir, "fit_object.rds"))

# Automate a contrast matrix
# Make sure baseline is the first group
niches <- paste0 ("T_",unique(bulk_anno$niches))

cell_n_tab <- table(bulk_anno$Sample_type_cell)%>%
  data.frame()

niches

# Compare each niche to all other niches

# Remove non-required groups
contrast_lines <- character()
i2 <- 0
for(i in 1:length(unique(bulk_anno$Manual_toplevel_pred))){
  
  cell_type <- unique(bulk_anno$Manual_toplevel_pred)[i]
  
  # Look over all the contrasts assuming the first is the baseline
  for(i3 in 1:length(niches)){
    
    contrast_1 <- niches[i3]
    contrast_2 <- niches[niches != contrast_1]
    contrast_2 <- paste0(contrast_2, "_", cell_type)
    contrast_1 <- paste0(contrast_1, "_", cell_type)
    
    # Make sure the contrast/cell combo is present
    contrast_2 <- contrast_2[contrast_2 %in% colnames(design)]
    
    if(contrast_1 %in% colnames(design) & (length(contrast_2) > 0)){
      
      i2 <- i2+1
      
      clen <- length(contrast_2)
      # Make a combined contrast2 (all conditions/length)
      contrast_2 <- paste0(contrast_2, collapse = "+")
      contrast_2 <- paste0("(",contrast_2, ")/",clen)
      
      # Test condition - baseline
      contrast_line <- paste0(contrast_1, "-", contrast_2)
      
      print(contrast_line)
      
      contrast_lines[i2] <- c(contrast_line)
      
    }
    
  }
  
}

cont.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                              (contrast_lines),levels=list(design))))

colnames(cont.matrix) <- gsub("-.*", "", colnames(cont.matrix))

unique(bulk_anno$Sample_type_cell)

contMat_2 <-  makeContrasts(TLS_B_T_vs_N =T_TLS_B - N_TLS_B, 
                            TLS_TCD4_T_vs_N =T_TLS_TCD4 - N_TLS_TCD4, 
                            levels = design)

cont.matrix <- cbind(cont.matrix, contMat_2)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

de_summary <- data.frame(summary(summa.fit), check.names = F)%>%
  dplyr::select(Contrast = 2, `Direction if significant` = 1, `Number of genes` = 3)%>%
  mutate(`Direction if significant` = factor(`Direction if significant`, levels = c("Up", "Down", "NotSig")))%>%
  arrange(`Direction if significant`, `Number of genes`)%>%
  write_csv(paste0(pseudobulk_dir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  filter(`Number of genes` >0)%>%
  group_by(Contrast)%>%
  mutate(Total_genes = sum(`Number of genes`))%>%
  ungroup()%>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] *-1))%>%
  arrange(-Total_genes)%>%
  mutate(Contrast = gsub("_", " ", Contrast))%>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`))+
  geom_bar(stat = "identity")+
  blank_theme

ggsave(filename =  paste0(pseudobulk_dir, "/Plots/Niche_DE_genes.pdf"), width = 15, height = 10)

# Make the output directory
system(paste0("mkdir -p ", pseudobulk_dir, "/toptables/"))
#system(paste0("rm ", pseudobulk_dir, "/toptables/*.csv"))

VOL <- paste0(pseudobulk_dir, "Volcano_plots/")
system(paste0("mkdir -p ", VOL))

# Make a vector of all contrast types to remove
# This works as long as the contrast names don't overlap the cell type names
all_cty <- paste0("T_",unique(bulk_anno$niches), collapse = "_|")
# Add on the last _
all_cty <- paste0(all_cty, "_")

VOL <- paste0(pseudobulk_dir, "glimma/volcano/")
MA <- paste0(pseudobulk_dir, "glimma/MA/")
system(paste0("mkdir -p ", VOL))
system(paste0("mkdir -p ", MA))

# Get all the toptables
for(contrast in colnames(cont.matrix)){
  
  output <- paste0(pseudobulk_dir, "toptables/", contrast, ".csv")
  
  toptable <- topTable(fit.cont,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("SYMBOL")%>%
    write_csv(output)
  
  conts <- gsub(all_cty, "", contrast)
  
  if(contrast == "TLS_B_T_vs_N"){
    conts <- "B"
  }else if(contrast == "TLS_TCD4_T_vs_N"){
    conts <- "TCD4"
  }
  
  bulk_filt <- bulk_anno%>%
    filter(Manual_toplevel_pred == conts)
  
  rnaseq_filt <- rnaseq[,bulk_filt$Niche_sample_cell]

  MA_fac <- factor(bulk_filt$Sample_type_cell, levels = unique(bulk_filt$Sample_type_cell))
  
  vol_save <- paste0(pseudobulk_dir, "glimma/volcano/",contrast,"_", "glimma_volcano.html")
  htmlwidgets::saveWidget(glimmaVolcano(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = round(rnaseq_filt$counts),
                                        dge = rnaseq_filt, groups = MA_fac), vol_save)
  
  MA_save <- paste0(pseudobulk_dir, "glimma/MA/",contrast,"_", "glimma_MA.html")
  htmlwidgets::saveWidget(glimmaMA(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = round(rnaseq_filt$counts),
                                        dge = rnaseq_filt, groups = MA_fac), MA_save)
  
}

# Drop the glimma files
system(paste0("rm -r ",pseudobulk_dir, "/glimma/*/*_files"))

# Compile the toptables to compare with CPDB
all_toptables <- list.files(paste0(pseudobulk_dir, "toptables/"), full.names = T)

tt_list <- list()
for(i in 1:length(all_toptables)){
  
  contrast <- gsub(".csv", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)%>%
    filter(adj.P.Val < 0.05)
  
  tt_list[[i]] <- tt
  
  # Do a volcano plot
  # tt <- tt %>%
  #   mutate(FC = ifelse(logFC > 0, "Up (FDR < 0.05)", "Down (FDR < 0.05)"))%>%
  #   mutate(FC = replace(FC, adj.P.Val >= 0.05, "FDR > 0.05"))%>%
  #   mutate(FC = factor(FC, levels = c("FDR > 0.05", "Down (FDR < 0.05)","Up (FDR < 0.05)")))%>%
  #   arrange(FC)%>%
  #   mutate(label = replace(SYMBOL,! (abs(logFC) > 2 | adj.P.Val < 0.05), NA))
  # 
  # custom_volcano <- ggplot(data = tt, aes(x = logFC, y = -log10(adj.P.Val), label = label,  colour= FC))+
  #   geom_point(aes())+
  #   blank_theme+
  #   geom_text_repel(max.overlaps =40,size=4, colour = "black", aes(fontface=3))+
  #   geom_vline(xintercept = 0,linetype = "dashed")+
  #   #xlim(-1*max(abs(toptable$logFC)), max(abs(toptable$logFC)))+
  #   scale_colour_manual(values = c("grey", "blue", "red"), drop = F)+
  #   guides(alpha = "none", colour = "none", size = "none")+
  #   labs(x = expression('Log'[2]*' fold change'), y = expression('-Log'[10]*' FDR'), colour = "Significance")+
  #   ggtitle(gsub("_", " ", contrast))
  # 
  # ggsave(filename =  paste0(pseudobulk_dir, "/Volcano_plots/", contrast, ".pdf"),plot = custom_volcano, width = 7, height = 7)
  # 
}

# Compile toptables and save the significant results
toptables_signif <- bind_rows(tt_list)%>%
  mutate(cell_type = gsub(".*Epithelium_", "", contrast))%>%
  group_by(SYMBOL)%>%
  mutate(Gene_DE_N_times = n())%>%
  ungroup()%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

# Read in the hallmark gene set
# Gene set collections
gsea_gmt_dir <- "~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
keep <- c(20, 31)
#keep <- c(31)
collections <- collections[keep]
# Add on my fetal signature
collections <- c(collections,"~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt")

# Make a directory
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)
  
  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix, mc.cores = 4)
  
}

# Save some gene sets to to plot/to use as DE in other datasets
# Compile the camera results
all_camera <- list.files(paste0(pseudobulk_dir,"gsea/camera/"), full.names = T)

clist <- list()
for(i in 1:length(all_camera)){
  
  contrast <- gsub("\\.csv", "", basename(all_camera[i]))
  
  tt <- read_csv(all_camera[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    select(-Contrast)
  
  clist[[i]] <- tt
  
}

# Compile toptables and save the significant results
camera_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

# Look at the TLSs a bit more
camera_GO_tls <- camera_compiled%>%
  filter(grepl("^GO", `Gene set`))%>%
  filter(grepl("TLS_", contrast))

# Look at the TLSs a bit more
camera_Hallmark_granulo <- camera_compiled%>%
  filter(grepl("^HALLMARK", `Gene set`))%>%
  filter(grepl("_Granulocyte_rich", contrast))

camera_GO_granulo <- camera_compiled%>%
  filter(grepl("^GO", `Gene set`))%>%
  filter(grepl("_Granulocyte_rich", contrast))

# Make some plots of TLS/neutrophil gene sets
granulo_macro_camera <- read_csv(paste0(pseudobulk_dir,"gsea/camera/c5.all.v2023.1_Granulocyte_rich_Macro.csv"))%>%
  filter(grepl("^GOBP_", `Gene set`))

plt <- gsea_plot (granulo_macro_camera, 20, "Granulo rich macrophage")
plt
  
# Make some plots of TLS/neutrophil gene sets
tls_mast_camera <- read_csv(paste0(pseudobulk_dir,"gsea/camera/c5.all.v2023.1_TLS_Mast.csv"))%>%
  filter(grepl("^GOBP_", `Gene set`))

plt <- gsea_plot (tls_mast_camera, 20, "TLS mast cell")
plt

# Make some plots of TLS/neutrophil gene sets
tls_macro_camera <- read_csv(paste0(pseudobulk_dir,"gsea/camera/c5.all.v2023.1_TLS_Macro.csv"))%>%
  filter(grepl("^GOBP_", `Gene set`))

plt <- gsea_plot (tls_macro_camera, 20, "TLS Macrophage")
plt

# Make some plots of TLS/neutrophil gene sets
tls_DC_camera <- read_csv(paste0(pseudobulk_dir,"gsea/camera/c5.all.v2023.1_TLS_DC.csv"))%>%
  filter(grepl("^GOBP_", `Gene set`))

plt <- gsea_plot (tls_macro_camera, 20, "TLS DC")
plt

# Make some plots of TLS/neutrophil gene sets
tls_DC_camera <- read_csv(paste0(pseudobulk_dir,"gsea/camera/c5.all.v2023.1_TLS_Mono.csv"))%>%
  filter(grepl("^GOBP_", `Gene set`))

plt <- gsea_plot (tls_macro_camera, 20, "TLS Mono")
plt

# Get the log2CPM to remove the slide level batch effect
batch_removed <- removeBatchEffect(v, batch = bulk_anno$Sample)

batch_removed[1:5,1:5]

# See the options for plotting
unique(bulk_anno$Manual_toplevel_pred)
unique(bulk_anno$niches)

bulk_filt <- bulk_anno%>%
  filter(Manual_toplevel_pred == "DC")%>%
  filter(niches %in% c("Immune_rich", "Stroma"))%>%
  arrange(niches, Sample)

# Get a gene set to plot in the correct order using the batch removed counts
gene_set_heat <- batch_removed[indexed$GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN,bulk_filt$Niche_sample_cell]%>%
  t()%>%
  scale()%>%
  t()

ha = HeatmapAnnotation(Niche = bulk_filt$niches,
  Celltype = bulk_filt$Manual_toplevel_pred,
                       Sample = bulk_filt$Sample,
                       Sample_type = bulk_filt$Sample_type,
                       Slide = bulk_filt$Slide,
                       col = list(Celltype = mycols,
                                  Niche = nichecols,
                         Sample_type = c("T" = "red", "N" = "blue")))

pdf(paste0(pseudobulk_dir, "gsea/DC genes heatmap.pdf"), width = 12, height = 7)
Heatmap(gene_set_heat, top_annotation = ha, 
        show_column_names = F, 
        cluster_columns = F)
dev.off()

# Do an RLE plot
# cpm_mat <- edgeR::cpm(rnaseq, log =T)
# mn  <- apply(cpm_mat, 1, median)
# rle <- data.frame(sweep(cpm_mat, MARGIN=1, STATS=mn, FUN='-'))
# List the cell types

unique(bulk_anno$Manual_toplevel_pred)
# Make a heatmap of the TLS niche
genes <- c("CCL19", "CCL21", "CCL2", "FGF7", "VCAM1", "ACKR1", "TXK", "CXCL9", "CXCL10")

genes <- genes[genes %in% rownames(batch_removed)]

TLS <- plot_hm(bulk_anno = bulk_anno, celltypes = unique(bulk_anno$Manual_toplevel_pred), 
               genes = genes,
               cpm = v$E, nichecols = nichecols)

pdf(paste0(pseudobulk_dir, "Plots/TLS genes heatmap.pdf"), width = 12, height = 5)
print(TLS)
dev.off()

# Plot some B cell TLS genes as a dotplot
unique(merged.integrated$niches)
unique(merged.integrated$Manual_toplevel_pred)

# Keep the B cells
B <- merged.integrated[,merged.integrated$Manual_toplevel_pred %in% c("B", "Plasma")]
B_niche <- B[,B$niches %in% c("TLS", "Epithelial_mass", "Normal_like")]

B_niche$Sample_cell_niche <- paste(B_niche$Sample_type, B_niche$Manual_toplevel_pred, B_niche$niches)

# Make a dotplot of some of the TLS genes
tls_genes <- c("IGHD", "PTHLH", "CCL19", "LAMP3", "CD19")

unique(B_niche$Sample_cell_niche)

ordering <- c("T B Epithelial_mass", "N B Epithelial_mass", "T Plasma Epithelial_mass", 
              "N Plasma Epithelial_mass", "T B Normal_like", "N B Normal_like", 
              "T Plasma Normal_like", "N Plasma Normal_like", "T B TLS", 
              "N B TLS", "T Plasma TLS", "N Plasma TLS")

B_niche$Sample_cell_niche <- factor(B_niche$Sample_cell_niche, levels = ordering)

DotPlot(B_niche, features = tls_genes, group.by = "Sample_cell_niche")

# Make a plot of some key T cell genes
T_cells <- merged.integrated[,merged.integrated$Manual_toplevel_pred %in% c("TCD4", "TCD8", "Tgd")]
T_niche <- T_cells[,T_cells$niches %in% c("TLS", "Epithelial_mass", "Normal_like")]

T_niche$Sample_cell_niche <- paste(T_niche$Sample_type, T_niche$Manual_toplevel_pred, T_niche$niches)

unique(T_niche$Sample_cell_niche)

ordering <-c("T TCD4 Epithelial_mass", "N TCD4 Epithelial_mass", "T TCD8 Epithelial_mass", 
             "N TCD8 Epithelial_mass", "T Tgd Epithelial_mass", "N Tgd Epithelial_mass", 
             "T TCD4 Normal_like", "N TCD4 Normal_like", "T TCD8 Normal_like", 
             "N TCD8 Normal_like", "T Tgd Normal_like", "N Tgd Normal_like", 
             "T TCD4 TLS", "N TCD4 TLS", "T TCD8 TLS", "N TCD8 TLS", 
             "T Tgd TLS", "N Tgd TLS")


T_niche$Sample_cell_niche <- factor(T_niche$Sample_cell_niche, levels = ordering)

# Make a dotplot of some of the TLS genes
T_genes <- c("TOX", "EOMES", "CD28", "CCL19")

DotPlot(T_niche, features = T_genes, group.by = "Sample_cell_niche")

bulk_anno_filt <- bulk_anno%>%
  filter(Manual_toplevel_pred %in% c("TCD4", "TCD8", "Tgd"))

plt <- plot_hm_niches(bulk_anno = bulk_anno_filt, celltypes = unique(bulk_anno$Manual_toplevel_pred), 
        genes = T_genes,
        cpm = v$E, nichecols = nichecols_underscore, niches_keep = unique(bulk_anno$niches), 
        donor_keep = unique(bulk_anno$Donor), cellcols = cell_type_colors)

print(plt)

# Read in the MMR atlas and plot some of the DC genes
mmr_atlas <- qread("/oldvol/apattison/Data/Reference/GSE178341_lognorm_annotated.qs")

# Read in the TLS toptable for fibroblasts
signif <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/Compiled_toptables_significant_genes.csv")
macro_tls <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/toptables/TLS_Macro.csv")
mono_tls <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/toptables/TLS_Mono.csv")
endo_tls <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/toptables/TLS_Endo.csv")
fibro_tls <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/toptables/TLS_Fibro.csv")
cd4_tls <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/toptables/TLS_TCD4.csv")
DC_tls <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/toptables/TLS_DC.csv")

DC_tls_signif <- DC_tls%>%
  filter(adj.P.Val < 0.05)

# Find the CD4 TLS genes that are not in the DC DE contrast
TLS_cd4_specific <- cd4_tls %>%
  filter(!SYMBOL %in% DC_tls_signif$SYMBOL)%>%
  filter(adj.P.Val < 0.05)

signif_TLS <- signif%>%
  filter(grepl("TLS_", contrast))

p1_a <- volcano_plot(fibro_tls, title = "TLS Fibro", top_n = 10)
p1_b <- volcano_plot(endo_tls, title = "TLS Endo", top_n = 10)
p1_c <- volcano_plot(cd4_tls, title = "TLS CD4 T", top_n = 10)
p1_d <- volcano_plot(DC_tls, title = "TLS DC", top_n = 10)

topleft <- plot_grid( p1_a, p1_b,p1_c,p1_d,   labels = c("i", "ii", "iii", "iv"), label_size = 8, 
                        nrow = 2)

DC <- mmr_atlas[,mmr_atlas$clMidwayPr %in% c("DC","Fibro", "Endo", "TCD4")]

DC <- DC[,DC$MMRStatus %in% c("MMRd", "MMRp")]

genes <- c("CD4", "IL7R","CCL5","CCR7", "SELL", "CCL19","LAMP3","CD38", "PLAC8", "CCL21", "ACKR1", "FGF7", "VCAM1", "TXK", "CLU")

bottomright <- DotPlot(object = DC, features = genes, group.by = "cl295v11SubFull",dot.scale = 3)+
  blank_theme+
  labs(y = "Cell type", x = "Probe")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'),legend.position = "top")+
  guides(size = guide_legend(title = "%\nexpressing"))+
  guides(color = guide_colorbar(title = "Expression"))

merged.integrated_TLS <- merged.integrated[,merged.integrated$Manual_toplevel_pred %in% c("DC","Fibro", "Endo", "TCD4")]
merged.integrated_TLS <- merged.integrated_TLS[,merged.integrated_TLS$Sample_type %in% c("T")]

merged.integrated_TLS$niche_cell <- paste0(merged.integrated_TLS$niches," ", merged.integrated_TLS$Manual_toplevel_pred)

topright <- DotPlot(object = merged.integrated_TLS, features = genes, group.by = "niche_cell", dot.scale = 3)+
  blank_theme+
  labs(y = "Cell type", x = "Probe")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "top")+
  guides(size = guide_legend(title = "%\nexpressing"))+
  guides(color = guide_colorbar(title = "Expression"))

topright

top <- plot_grid(topleft, topright,labels =  c("a", "b"), label_size = 8)
bottom <- plot_grid(bottomright, bottomright, rel_widths = c(1.5, 1), 
                    labels = c("c", "d"), label_size = 8)

Figure_TLS <- plot_grid(top,bottom, nrow = 2)
ggsave(plot = Figure_TLS,"./Results_paper/Plots/Figure TLS.pdf", width = 170, height = 200, units = "mm")

# Find some genes that are highly expressed in specific FOVs
CCL19 <- gene_fov_finder(gene = "CCL19", seurat_obj = merged.integrated)

# Read in the Griffith data
merged.integrated_g <- qread("./Intermediate/lognorm_merged_integrated_annotated_g.qs")

# Can't see too mcuh of this in the entire dataset
CCL19_g <- gene_fov_finder(gene = "CCL19", seurat_obj = merged.integrated_g)


