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
library(parallel)
library(GSA)
library(DescTools)
library(cowplot)

# What to do here:
#https://support.bioconductor.org/p/125489/#125602
# Set a homedir for the analysis so as to no use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

merged <- qread("./Intermediate/lognorm_merged_integrated_annotated.qs")

# Source the R functions I have created and plot themes/colours
source("~/Code/Spatial/Paper_code_v2/functions.R")

# Set the output dir
pseudobulk_dir <- "./Results_paper/analyses/Pseudobulk_T_vs_N_HQ/"

system(paste0("mkdir -p ", pseudobulk_dir, "Plots/"))
system(paste0("mkdir -p ",pseudobulk_dir, "toptables/"))
system(paste0("mkdir -p ",pseudobulk_dir, "tables/"))

# Make a pseudobulk object ----
to_bulk_anno <- merged@meta.data%>%
  dplyr::select(Barcode, Manual_toplevel_pred, Sample, Slide, fov)%>%
  mutate(Sample_type =ifelse(grepl("T", Sample), "T", "N"))%>%
  mutate(Sample_name = Sample)%>%
  mutate(Manual_toplevel_pred = gsub(" |-", "_", Manual_toplevel_pred))%>%
  mutate(donor_sampletype_cell = paste0(Sample, "_", Sample_type, "_", Manual_toplevel_pred, "_", Slide))

# Keep only conditions where we have a decent number of cells
to_drop <- to_bulk_anno %>%
  group_by(donor_sampletype_cell)%>%
  summarise(count = n())%>%
  # Keep cell types with >10 cells per sample
  mutate(drop = count <= 10)

keep <- to_drop%>%
  filter(!drop)

table(to_drop$drop)

to_bulk_anno <- to_bulk_anno %>%
  filter(donor_sampletype_cell %in% keep$donor_sampletype_cell)

cell_type_counts <- to_bulk_anno %>%
  group_by(Sample_name, Sample_type, Manual_toplevel_pred, Slide)%>%
  summarise(count = n())%>%
  group_by(Sample_name,Slide)%>%
  mutate(total = sum(count))%>%
  ungroup()%>%
  mutate(Pct = count/total*100)%>%
  arrange(-Pct)%>%
  mutate(Sample_name = factor(Sample_name, levels = unique(Sample_name)))%>%
  mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))%>%
  dplyr::rename(Count = count, Total_cells_per_sample = total, Percent_total_cells = Pct)%>%
  write_csv(paste0(pseudobulk_dir, "Sample_cell_type_compositions.csv"))

# Make sure the objects have the same cells
pb_counts <- merged@assays$RNA$counts

# Get the number of samples to test
nums <- 1:length(unique(to_bulk_anno$donor_sampletype_cell))

# Parallelise making pseudobulk counts
pb_list <- mclapply(X = nums, FUN = get_pb, to_bulk_anno,aggr_col = "donor_sampletype_cell", mc.cores = 15)

# Bind the PB counts
bound <- bind_cols(pb_list)

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
  distinct(donor_sampletype_cell ,.keep_all = T)

# Contrast only one sample replicate to avoid repeated contrasts from the same sample
reptab <-condensed_SC_anno %>%
  select(Sample, Slide)%>%
  distinct()%>%
  filter(!duplicated(Sample))%>%
  mutate(Sample_slide = paste0(Sample, "_", Slide))

bulk_anno <- data.frame(donor_sampletype_cell = colnames(bound_counts), check.names = F)%>%
  left_join(condensed_SC_anno)%>%
  # Drop any cell types that could not be identified from scmatch
  filter(!is.na(Manual_toplevel_pred))%>%
  mutate(Sample_slide = paste0(Sample, "_", Slide))%>%
  # Only contrast one rep of each sample
  mutate(Sample_type_old = Sample_type)%>%
  mutate(Sample_type = replace(Sample_type, !Sample_slide %in% reptab$Sample_slide, "Duplicate"))%>%
  mutate(Sample_type = replace(Sample_type, Sample_type == "Duplicate", paste0(Sample_type[Sample_type == "Duplicate"], Sample_type_old[Sample_type == "Duplicate"])))%>%
  mutate(Sample_type_cell = paste0(Sample_type, "_", Manual_toplevel_pred))%>%
  mutate(Donor = gsub("T|N","", Sample_name))%>%
  mutate(Slide_fov = paste0(Slide, "_", fov))%>%
  # Save the Pseudobulk annotation
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# Make a DGElist
rnaseq <- DGEList(bound_counts[,bulk_anno$donor_sampletype_cell])

# Remove weird characters
rownames(rnaseq) <- gsub("\\+", "_pos_", rownames(rnaseq))
colnames(rnaseq) <- gsub("\\+", "_pos_", colnames(rnaseq))
colnames(rnaseq) <- gsub("\\'|-|\\(|\\)", "_", colnames(rnaseq))
bulk_anno$donor_sampletype_cell <- gsub("\\+", "_pos_", bulk_anno$donor_sampletype_cell)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$donor_sampletype_cell) == length(bulk_anno$donor_sampletype_cell)

design <- model.matrix(~0 + Sample_type_cell, data = bulk_anno)
# Neaten up design row and colnames
colnames(design) <- gsub("Sample_type_cell|\\+", "", colnames(design))
colnames(design) <- gsub("\\'|-|\\(|\\)", "_", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

hist(edgeR::cpm(rnaseq), breaks = 20000, xlim = c(0,1000))
  

# Keep genes with CPM >200 the negative probe expression in >10% of samples (all of them)
keep <- rowSums(edgeR::cpm(rnaseq)>200)>(0.1 * nrow(design))
table(keep)
rownames(rnaseq)[!keep]
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]

rnaseq <- calcNormFactors(rnaseq)

# No need to filter as all probes are expressed from single-cell
# Get the TMM normalised CPMs in case anyone wants to work with it
cpm_save <- edgeR::cpm(rnaseq, log =T)%>%
  data.frame()%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_log2_CPMs.csv"))

# Do a glimma MDS of batch removed counts
system(paste0("mkdir -p ", pseudobulk_dir,"glimma/mds/"))

# Save and MDS plot per cell type
for(celltype in unique(bulk_anno$Manual_toplevel_pred)){
  
  print(celltype)
  
  mds_save <- paste0(paste0(pseudobulk_dir,"glimma/mds/", celltype, "_MDS.html"))
  
  htmlwidgets::saveWidget(glimmaMDS(rnaseq[,grepl(celltype, colnames(rnaseq))], groups = bulk_anno[grepl(celltype, colnames(rnaseq)),],
                                    labels = bulk_anno$donor_sampletype_cell[grepl(celltype, colnames(rnaseq))]), mds_save)
  
}

# Voomlmfit for single cell
fit <- voomLmFit(rnaseq, design = design, plot = T, block = bulk_anno$Slide)

# Automate a contrast matrix
# Make sure baseline is the first group
contrasts <- unique(bulk_anno$Sample_type)
contrasts_manual <- "T-N"
# Remove non-required groups
contrast_lines <- character()
i2 <- 0
for(i in 1:length(unique(bulk_anno$Manual_toplevel_pred))){
  
  cell_type <- unique(bulk_anno$Manual_toplevel_pred)[i]
  
  # Look over all the contrasts assuming the first is the baseline
  for(i3 in 1:length(contrasts_manual)){
    
    contrast_line <- contrasts_manual[i3]
    contrast_2 <- gsub("-.*", "",contrast_line)
    contrast_1 <- gsub(".*-", "",contrast_line)
    contrast_1 <- paste0(contrast_1, "_", cell_type)
    contrast_2 <- paste0(contrast_2, "_", cell_type)
    
    if(contrast_1 %in% colnames(design) & contrast_2 %in% colnames(design)){
      
      i2 <- i2+1
      # Test condition - baseline
      contrast_line <- paste0(contrast_2, "-", contrast_1)
      
      print(contrast_line)
      
      contrast_lines[i2] <- c(contrast_line)
      
    }
    
  }
  
}

cont.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                              (contrast_lines),levels=list(design))))

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

de_summary <- data.frame(summary(summa.fit), check.names = F)%>%
  dplyr::select(Contrast = 2, `Direction if significant` = 1, `Number of genes` = 3)%>%
  mutate(`Direction if significant` = factor(`Direction if significant`, levels = c("Up", "Down", "NotSig")))%>%
  arrange(`Direction if significant`, -`Number of genes`)%>%
  write_csv(paste0(pseudobulk_dir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  filter(`Number of genes` >10)%>%
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

VOL <- paste0(pseudobulk_dir, "glimma/volcano/")
system(paste0("mkdir -p ", VOL))

# This works as long as the contrast names don't overlap the cell type names
all_cty <- paste0(unique(bulk_anno$Sample_type), collapse = "_|")
# Add on the last _
all_cty <- paste0(all_cty, "_")

# Get all the toptables
for(contrast in colnames(cont.matrix)){
  
  output <- paste0(pseudobulk_dir, "toptables/", contrast, ".csv")
  
  toptable <- topTable(fit.cont,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("SYMBOL")%>%
    write_csv(output)
  
  conts <- gsub(all_cty, "", contrast)
  
  # Get back to just the cell type
  cont_first <- gsub("-.*", "", conts)
  cont_second <- gsub(".*-", "", conts)
  
  conts_both <- paste0(cont_first, "|", cont_second)
  
  rnaseq_filt <- rnaseq[,grepl(conts_both, colnames(rnaseq))]
  
  bulk_filt <- bulk_anno%>%
    filter(donor_sampletype_cell %in% colnames(rnaseq_filt))
  
  MA_fac <- factor(bulk_filt$Sample_type_cell, levels = unique(bulk_filt$Sample_type_cell))
  
  vol_save <- paste0(pseudobulk_dir, "glimma/volcano/",contrast,"_", "glimma_volcano.html")
  htmlwidgets::saveWidget(glimmaVolcano(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = round(rnaseq_filt$counts),
                                        dge = rnaseq_filt, groups = MA_fac), vol_save)
}

# Drop the glimma files
system(paste0("rm -r ",pseudobulk_dir, "/glimma/*/*_files"))

# Compile the toptables to compare with CPDB
all_toptables <- list.files(paste0(pseudobulk_dir, "toptables/"), full.names = T)

# Read in the high ambient RNA genes
#high_ambient <- read_csv("/homevol/apattison/Data/Spatial/Results/Tables/High_ambient_RNA_genes.csv")
de_ambient <- read_csv("/homevol/apattison/Data/Spatial/Results/Tables/T_vs_N_ambient_DE.csv")

VOL <- paste0(pseudobulk_dir, "Volcano_plots/")
system(paste0("mkdir -p ", VOL))

tt_list <- list()
for(i in 1:length(all_toptables)){
  
  contrast <- gsub(".csv", "", basename(all_toptables[i]))
  
  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
  # Do a volcano plot
  tt <- tt %>%
    # Drop the high ambient genes
    filter(!SYMBOL %in% de_ambient$SYMBOL)%>%
    filter(adj.P.Val < 0.05)%>%
    mutate(FC = ifelse(logFC > 0, "Up (FDR < 0.05)", "Down (FDR < 0.05)"))%>%
    mutate(FC = replace(FC, adj.P.Val >= 0.05, "FDR > 0.05"))%>%
    mutate(FC = factor(FC, levels = c("FDR > 0.05", "Down (FDR < 0.05)","Up (FDR < 0.05)")))%>%
    arrange(FC)%>%
    mutate(label = replace(SYMBOL,! (abs(logFC) > 2 | adj.P.Val < 0.05), NA))
  
  custom_volcano <- ggplot(data = tt, aes(x = logFC, y = -log10(adj.P.Val), label = label,  colour= FC))+
    geom_point(aes())+
    blank_theme+
    geom_text_repel(max.overlaps =40,size=4, colour = "black", aes(fontface=3))+
    geom_vline(xintercept = 0,linetype = "dashed")+
    xlim(-1*max(abs(tt$logFC)), max(abs(tt$logFC)))+
    scale_colour_manual(values = c("grey", "blue", "red"), drop = F)+
    guides(alpha = "none", colour = "none", size = "none")+
    labs(x = expression('Log'[2]*' fold change'), y = expression('-Log'[10]*' FDR'), colour = "Significance")+
    ggtitle(gsub("_", " ", contrast))
  
  ggsave(filename =  paste0(pseudobulk_dir, "/Volcano_plots/", contrast, ".pdf"),plot = custom_volcano, width = 7, height = 7)
  
}

# Compile toptables and save the significant results
toptables_compiled <- bind_rows(tt_list)%>%
  mutate(cell_type = gsub(".*-", "", contrast))%>%
  mutate(cell_type = gsub("N_", "", cell_type))

toptables_signif <- toptables_compiled %>%
  mutate(High_ambient = ifelse(SYMBOL %in% de_ambient$SYMBOL, "Possible ambient DE", "No ambient DE"))%>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  mutate(FC_direction = ifelse(logFC > 0, "Up", "Down"))%>%
  group_by(SYMBOL,FC_direction)%>%
  mutate(Gene_DE_N_times = n())%>%
  ungroup()%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

# GSEA analysis
# Make a voom object for camera analysis with blocking included
v <- voom(rnaseq, block = bulk_anno$Slide, correlation = fit$correlation)

# Define a function to shut up some other functions
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Read in the hallmark gene set
# Gene set collections
gsea_gmt_dir <- "~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
keep <- c(5,1,20, 28,31,9,16)
#keep <- c(31)
collections <- collections[keep]
# Add on my fetal signature
collections <- c(collections,"~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt")

# Make a directory
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)
  
  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix, mc.cores = 3)
  
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
  # Fix the cell type naming
  mutate(cell_type = gsub(".*-N_", "", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

# Summarise the CPM for each cell type for
cpm <- edgeR::cpm(rnaseq, log =T)

sum(colnames(cpm) == bulk_anno$donor_sampletype_cell) == ncol(cpm)

t_vs_n_gene (cpm, "SPP1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
t_vs_n_gene (cpm, "CRP", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
t_vs_n_gene (cpm, "IL6", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
t_vs_n_gene (cpm, "CXCL8", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
t_vs_n_gene (cpm, "KRT23", bulk_anno, paste0(pseudobulk_dir,"Plots/"))

# Do an outlier analysis per cell type
bulk_ordered <- bulk_anno%>%
  arrange(Sample_type,Manual_toplevel_pred,Slide)

outlier_list <- list()
for(i in 1:length(unique(bulk_ordered$Manual_toplevel_pred))){
  
  cell_type <- unique(bulk_ordered$Manual_toplevel_pred)[i]
  
  if(cell_type == "Hepatocytes"){
    
    next
    
  }
  
  bulk_ordered_celltype <- bulk_ordered %>%
    filter(Manual_toplevel_pred == cell_type)
  
  # Grab some of the outlier genes
  outliers <- cpm[,bulk_ordered_celltype$donor_sampletype_cell]%>%
    t()%>%
    # Use the robust z score
    scale()%>%
    t()%>%
    data.frame(check.rows = F, check.names = F)%>%
    rownames_to_column("Gene")%>%
    gather(donor_sampletype_cell, Z_score, -Gene)%>%
    filter(abs(Z_score) >2.5)%>%
    left_join(bulk_ordered_celltype)
  
  outlier_list[[i]] <- outliers
  
}

all_outs <- bind_rows(outlier_list)%>%
  arrange(-Z_score)

bulk_stromal <- bulk_anno%>%
  arrange(Sample_type_old,Sample, Manual_toplevel_pred,Slide)%>%
  filter(Manual_toplevel_pred %in% c("Epi", "Fibro", "Granulo", "Mono", "Macro", "Peri", "Endo"))

ha = HeatmapAnnotation(Celltype = bulk_stromal$Manual_toplevel_pred,
                       Sample = bulk_stromal$Sample,
                       "Sample type" = bulk_stromal$Sample_type_old,
                       col = list(Celltype = cell_type_colors, Sample = Sample_cols,
                                  "Sample type" = sample_type_cols))

# Look at some neutrophil genes
neu_genes <- c("CXCL1","CXCL2", "CXCL3", "CXCL8", "CXCL5", "CXCL12","CXCR4", "IL6", "SPP1", "HCAR1")
neu_genes <- neu_genes[neu_genes %in% rownames(cpm)]
scaled <- cpm%>%
  t()%>%
  scale()%>%
  t()

to_plot <- scaled[neu_genes,bulk_stromal$donor_sampletype_cell]

pdf(paste0(pseudobulk_dir, "Plots/Neutrophil chemoattractants.pdf"), width = 9, height = 5)
Heatmap(to_plot, top_annotation = ha, show_column_names = F, cluster_columns = F, cluster_rows = F, name  = "Z score")
dev.off()

# Make some volcano plots to compare T vs N for cancer cells
T_vs_N <- read_csv("./Results_paper/analyses/Pseudobulk_T_vs_N_HQ/toptables/T_Epi-N_Epi.csv")
T_vs_N_V2 <- read_csv("./Results_paper/analyses/Pseudobulk_T_vs_N_G_HQ/toptables/T_Epi-N_Epi.csv")

# Read in T vs N for the MMR atlas
T_vs_N_MMR <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/toptables/T_Epi-N_Epi.csv")

# Compare the p values and fold changes
T_vs_N_MMR_join <- T_vs_N_MMR%>%
  select(SYMBOL, mmr_logFC = logFC, mmr_adj.P.Val = adj.P.Val)

T_vs_N_MMR_plot <- T_vs_N%>%
  left_join(T_vs_N_MMR_join)%>%
  # Drop MMR NAs
  filter(!is.na(mmr_adj.P.Val))%>%
  mutate(SYMBOL = replace(SYMBOL, adj.P.Val > 0.05, NA))%>%
  mutate(FC = ifelse(logFC > 0, "Up", "Down"))%>%
  mutate(FC = replace(FC, adj.P.Val > 0.05, "NS"))

T_vs_N_MMR_plot_signif <- T_vs_N_MMR_plot%>%
  filter(adj.P.Val < 0.05)

cor <-cor.test(T_vs_N_MMR_plot_signif$logFC, T_vs_N_MMR_plot_signif$mmr_logFC)

p1 <- ggplot(data = T_vs_N_MMR_plot, aes(x = logFC, y = mmr_logFC, label = SYMBOL, colour = FC))+
  geom_point(size = 0.5)+
  geom_text_repel(max.overlaps = 30, force = 10, fontface = "italic", size =2.5)+
  scale_colour_manual(values = c("Blue", "Grey", "Red"))+
  guides(colour = "none")+
  blank_theme+
  annotate("text", x = -2, y =8, label = paste0("Pearson's r = ", round(cor$estimate,2), "\n (significant genes)"), size = 2.5)+
  labs(x =  expression('Log'[2]*' FC CosMx'), y = expression('Log'[2]*' FC MMR atlas'), colour = "FDR < 0.05\n (CosMx)")

# V1 volcano plot
T_vs_N_vol <- volcano_plot(T_vs_N, title = "T vs N (CosMx V1)", top_n = 10)
T_vs_N_vol

# V2 volcano plot
T_vs_N_V2_vol <- volcano_plot(T_vs_N_V2, title = "T vs N (CosMx V2)", top_n = 10)
T_vs_N_V2_vol

# MMR volcano plot
T_vs_N_MMR_vol <- volcano_plot(T_vs_N_MMR, title = "T vs N (MMR atlas)", top_n = 10)
T_vs_N_MMR_vol

# Make some GSEA plos of tumour stem cell expression
T_vs_N_MMR_GSEA <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Compiled_significant_gene_sets_camera.csv")

T_vs_N_MMR_fetal_GSEA <- T_vs_N_MMR_GSEA%>%
  filter(grepl("Fetal_gene_sigs", contrast))%>%
  filter(cell_type == "Epi")

T_vs_N_MMR_HALLMARK_GSEA <- T_vs_N_MMR_GSEA%>%
  filter(grepl("HALLMARK", `Gene set`))%>%
  filter(cell_type == "Epi")

gsea_plot_fetal <- gsea_plot(T_vs_N_MMR_fetal_GSEA, top_n = 4, title = "T vs N MMR atlas", xlab = "Stem gene set")
gsea_plot_hallmark <- gsea_plot(T_vs_N_MMR_HALLMARK_GSEA, top_n = 20, title = "T vs N MMR atlas")

mmr_p_vs_d <- read_csv("./Results_paper/analyses/MMR_atlas_T_split_vs_N/Compiled_toptables_significant_genes.csv")%>%
  filter(grepl("-MMRp_Epi",contrast))

mmr_p_vs_d_gsea <- read_csv("./Results_paper/analyses/MMR_atlas_T_split_vs_N/Compiled_significant_gene_sets_camera.csv")%>%
  filter(grepl("-MMRp",contrast))

mmr_p_vs_d_vol <- volcano_plot(mmr_p_vs_d, title = "MMRd vs MMRp (MMR atlas)", top_n = 20)

# Check the expression of some stem cell genes in the MMR atlas
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

stem_genes <- GSA.read.gmt("/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt")

stem_genes$genesets
stem_genes$geneset.names

T_vs_N_MMR_genes_stem <- unlist(stem_genes$genesets)%>%
  unique()

T_vs_N_MMR_stem <- T_vs_N_MMR%>%
  filter(SYMBOL %in% T_vs_N_MMR_genes_stem)%>%
  filter(adj.P.Val < 0.01)%>%
  top_n(30)%>%
  arrange(-t)

# Add some selected known stem cell genes that don't rank as highly but are still there

genes <- T_vs_N_MMR_stem$SYMBOL[1:30]

mmr_atlas_epi <- mmr_atlas[,mmr_atlas$clMidwayPr == "Epi"]
mmr_atlas_epi$MMRStatus <- replace(mmr_atlas_epi$MMRStatus, is.na(mmr_atlas_epi$MMRStatus), "Normal")

stem_dot <- DotPlot(object = mmr_atlas_epi, features = genes, group.by = "MMRStatus",dot.scale = 3)+
  blank_theme+
  labs(y = "MMR", x = "Gene")+
  theme(axis.text.x = element_text(angle = 45,face = "italic",vjust = 1, hjust=1),
        legend.key.size = unit(2, 'mm'), legend.position = "top")+
  guides(size = guide_legend(title = "%\nexpressing"))+
  guides(color = guide_colorbar(title = "Expression"))

# Make a plot of the numbers of genes in each cell type for MMRp/d
mmr_p_vs_d_all <- read_csv("./Results_paper/analyses/MMR_atlas_T_split_vs_N/Compiled_toptables_significant_genes.csv")
mmr_p_vs_d_DE <- mmr_p_vs_d_all%>%
  filter(grepl("-MMRp",contrast))%>%
  mutate(cell_type = gsub(".*MMRp_", "", cell_type))%>%
  mutate(Direction = ifelse(logFC > 0 , "Up", "Down"))%>%
  group_by(cell_type, Direction)%>%
  summarise(Count = n())%>%
  ungroup()%>%
  mutate(Count = replace(Count, Direction == "Down", Count[Direction == "Down"] *-1))%>%
  arrange(-Count)%>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

gene_numbers <- ggplot(data = mmr_p_vs_d_DE, aes(y = cell_type, x = Count, fill = Direction))+
  geom_bar(stat = "identity")+
  blank_theme+
  scale_fill_manual(values = c("blue", "red"), drop = F)

# Make a figure our of the volcano plots
top_row <- plot_grid(T_vs_N_vol, T_vs_N_V2_vol,T_vs_N_MMR_vol, p1, labels = c("a", "b", "c", "d"), 
                        label_size = 8, ncol = 4, align = "h", axis = "bt")

middle_middle <- plot_grid(gsea_plot_fetal, stem_dot,ncol = 1, labels = c("f", "g"), label_size = 8)

middle_left <- plot_grid(gsea_plot_hallmark, labels = c("e"), 
                        label_size = 8, nrow = 1, align = "h", axis = "bt")
middle_row <- plot_grid(middle_left, middle_middle, nrow = 1, rel_widths = c(0.9,1))
bottom_row <- plot_grid(gene_numbers, gene_numbers, nrow = 1, rel_widths = c(1.5,1), labels = c("h", "i"), 
                        label_size = 8)
all_t_vs_n <- plot_grid(top_row, middle_row, bottom_row, nrow = 3)

ggsave(plot = all_t_vs_n,"./Results_paper/Plots/Figure 3 T vs N.pdf", width = 170, height = 200, units = "mm")

# Make a plot of the DE genes for each cell type for T vs N for CosMx V1
all_t_vs_n <- list.files("./Results_paper/analyses/Pseudobulk_T_vs_N_HQ/toptables/", full.names = T)
all_t_vs_n <- all_t_vs_n[!grepl("Granulo",all_t_vs_n)]

plotlist <- list()
for(i in 1:length(all_t_vs_n)){
  
  DE <- read_csv(all_t_vs_n[i])
  
  title <- gsub("_", " ", basename(all_t_vs_n[i]))
  title <- gsub(".csv", "", title)
  
  plt <- volcano_plot(DE, title = title, top_n = 20, text_size = 2)
  
  plt <- rasterize(plt, layers='Point', dpi=100)
  
  plotlist[[i]]<- plt
}

grd <- plot_grid(plotlist = plotlist, ncol = 3,labels = letters[1:length(plotlist)], label_size = 8)
ggsave(plot = grd,"./Results_paper/Plots/Figure S14.pdf", width = 170, height = 400, units = "mm")

# Make a plot of the DE genes for each cell type for T vs N
all_mmr_t_vs_n <- list.files("./Results_paper/analyses/MMR_atlas_T_vs_N/toptables/", full.names = T)

plotlist <- list()
for(i in 1:length(all_mmr_t_vs_n)){
  
  DE <- read_csv(all_mmr_t_vs_n[i])
  
  title <- gsub("_", " ", basename(all_mmr_t_vs_n[i]))
  title <- gsub(".csv", "", title)
  
  plt <- volcano_plot(DE, title = title, top_n = 20, text_size = 2)
  
  plt <- rasterize(plt, layers='Point', dpi=100)
  
  plotlist[[i]]<- plt
}

grd <- plot_grid(plotlist = plotlist, ncol = 3,labels = letters[1:length(plotlist)], label_size = 8)
ggsave(plot = grd,"./Results_paper/Plots/Figure S15.pdf", width = 170, height = 400, units = "mm")

# Make a plot of the DE genes for each cell type for mmrp vs d
all_mmr_t_vs_n <- list.files("./Results_paper/analyses/MMR_atlas_T_split_vs_N/toptables/", full.names = T)
all_mmr_t_vs_n <-all_mmr_t_vs_n[grepl("-MMRp",all_mmr_t_vs_n)] 
plotlist <- list()
for(i in 1:length(all_mmr_t_vs_n)){
  
  DE <- read_csv(all_mmr_t_vs_n[i])
  
  title <- gsub("_", " ", basename(all_mmr_t_vs_n[i]))
  title <- gsub(".csv", "", title)
  
  plt <- volcano_plot(DE, title = title, top_n = 20, text_size = 2)
  
  plt <- rasterize(plt, layers='Point', dpi=100)
  
  plotlist[[i]]<- plt
}

grd <- plot_grid(plotlist = plotlist, ncol = 3,labels = letters[1:length(plotlist)], label_size = 8)
ggsave(plot = grd,"./Results_paper/Plots/Figure S16.pdf", width = 170, height = 400, units = "mm")

# Make a plot of the DE genes for each cell type for TLS high vs low
all_TLS <- list.files("./Results_paper/analyses/MMR_atlas_TLS/toptables/", full.names = T)
all_TLS <- all_TLS[!grepl("*ILC*",all_TLS)]

plotlist <- list()
for(i in 1:length(all_TLS)){
  
  DE <- read_csv(all_TLS[i])
  
  title <- gsub("_", " ", basename(all_TLS[i]))
  title <- gsub(".csv", "", title)
  title <- gsub("TLS", "LA", title)
  
  plt <- volcano_plot(DE, title = title, top_n = 20, text_size = 2)
  
  plt <- rasterize(plt, layers='Point', dpi=100)
  
  plotlist[[i]]<- plt
}

grd <- plot_grid(plotlist = plotlist, ncol = 3,labels = letters[1:length(plotlist)], label_size = 8)
ggsave(plot = grd,"./Results_paper/Plots/Figure S TLS.pdf", width = 170, height = 400, units = "mm")

# Make a plot of the DE genes for each cell type for granulo high vs low
all_neu <- list.files("./Results_paper/analyses/MMR_atlas_Neu/toptables/", full.names = T)

plotlist <- list()
for(i in 1:length(all_neu)){
  
  DE <- read_csv(all_neu[i])
  
  title <- gsub("_", " ", basename(all_neu[i]))
  title <- gsub(".csv", "", title)
  
  plt <- volcano_plot(DE, title = title, top_n = 20, text_size = 2)
  
  plt <- rasterize(plt, layers='Point', dpi=100)
  
  plotlist[[i]]<- plt
}

grd <- plot_grid(plotlist = plotlist, ncol = 3,labels = letters[1:length(plotlist)], label_size = 8)
ggsave(plot = grd,"./Results_paper/Plots/Figure S Neu.pdf", width = 170, height = 400, units = "mm")




