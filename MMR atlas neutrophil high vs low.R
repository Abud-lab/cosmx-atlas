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
library(GSA)
library(parallel)
library(ActivePathways)
library(rstatix)
library(ggrepel)
library(cowplot)
library(GSVA)

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_v2/functions.R")

# Set a homedir for the analysis so as to no use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

# Read in the MMR atlas with logtransformed counts
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

pseudobulk_dir <- "./Results_paper/analyses/MMR_atlas_Neu/"

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

# Take a loot at the MMR atlas and see if you can do any better
# with processing

mmr_atlas$Sample_name <- mmr_atlas$PatientTypeID

mmr_atlas@assays$RNA@data[1:20,1:20]

mmr_atlas$Manual_toplevel_pred <- mmr_atlas$clMidwayPr

# Study annotations
unique(mmr_atlas$clMidwayPr)

granulo_tab <- mmr_atlas@meta.data%>%
  select(Sample_name = PatientTypeID, Manual_toplevel_pred,MMRStatus)%>%
  mutate(T_N = gsub(".*_", "", Sample_name))%>%
  mutate(isgranulo = ifelse(Manual_toplevel_pred == "Granulo", 1, 0))%>%
  group_by(Sample_name,MMRStatus)%>%
  mutate(Total = n())%>%
  group_by(Sample_name, Total,MMRStatus)%>%
  summarise(Total_granulo = sum(isgranulo))%>%
  ungroup()%>%
  mutate(Pct_granulo = Total_granulo/Total*100)%>%
  mutate(MMRStatus = replace(MMRStatus, is.na(MMRStatus), "Normal sample"))%>%
  mutate(MMRStatus = factor(MMRStatus, levels = c("Normal sample", "MMRp", "MMRd")))%>%
  group_by(MMRStatus)%>%
  mutate(Mean_granulo_pct = mean(Pct_granulo))%>%
  ungroup()

# Save the granulocyte tab
write_csv(granulo_tab, "./Results_paper/Tables/MMR atlas granulo pcts.csv")

hist(granulo_tab$Pct_granulo)

stat.test <- granulo_tab %>% 
  wilcox_test(Pct_granulo ~ MMRStatus) %>%
  add_significance()
stat.test

Figure_4_a <- ggplot(data = granulo_tab, aes(x = MMRStatus, y = Pct_granulo))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(height = 0, width = 0.3)+
  labs(x = "MMR status", y = "Percent granulocytes")+
  blank_theme+
  geom_hline(yintercept = 2, linetype = 2)
  
Figure_4_a

# Samples with > 2% granulocyte infiltration
granulo_keep <- granulo_tab$Sample_name[granulo_tab$Pct_granulo >2]

granulo_pcts <- granulo_tab%>%
  select(Sample_name, Pct_granulo)

granulo_pcts_tumour <- granulo_pcts%>%
  filter(!grepl("_N", Sample_name))%>%
  mutate(High = Pct_granulo >2)

table(granulo_pcts_tumour$High)
  
to_bulk_anno <- mmr_atlas@meta.data%>%
  dplyr::select(Barcode = sampleID, Manual_toplevel_pred, Sample_name = PatientTypeID,
                Sample_type = SPECIMEN_TYPE, MMRStatus, Sex, TumorStage)%>%
  mutate(T_N = gsub(".*_", "", Sample_name))%>%
  mutate(Sample_cell = paste0(Sample_name, "_", Manual_toplevel_pred))%>%
  filter(T_N == "T")%>%
  mutate(Sample_type = ifelse(Sample_name %in% granulo_keep, "Granulo_high", "Granulo_low"))%>%
  left_join(granulo_pcts)

# Keep only conditions where we have a decent number of cells
to_drop <- to_bulk_anno %>%
  group_by(Sample_cell,Sample_type, T_N)%>%
  summarise(count = n())%>%
  ungroup()

keep <- to_drop%>%
  filter(count >= 10)%>%
  filter(!grepl("_NA", Sample_type))

to_bulk_anno <- to_bulk_anno %>%
  filter(Sample_cell %in% keep$Sample_cell)

# Plot the niche composition of each sample
cell_type_counts <- to_bulk_anno %>%
  group_by(Sample_name, Sample_type, Manual_toplevel_pred)%>%
  summarise(count = n())%>%
  group_by(Sample_name)%>%
  mutate(total = sum(count))%>%
  ungroup()%>%
  mutate(Pct = count/total*100)%>%
  arrange(-Pct)%>%
  dplyr::rename(Count = count, Total_cells_per_sample = total, Percent_total_cells = Pct)%>%
  write_csv(paste0(pseudobulk_dir, "Sample_compositions.csv"))

# Get a big vector of different colours
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Cell type composition plot
plt <- ggplot(data = cell_type_counts, aes(x = Sample_name, y = Percent_total_cells, fill =Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  facet_wrap(~Sample_type, scales = "free_x")+
  labs(x = "Sample", y = "Percent of total", fill = "Niche")+
  blank_theme+
  scale_fill_manual(values = col_vector)+
  guides(fill = guide_legend(ncol = 2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plt

save <- paste0(pseudobulk_dir, "/Plots/Per sample cell type composition.pdf")
ggsave(filename = save, plot = plt, width = 8, height = 5)

# Get the raw counts
pb_counts <- mmr_atlas@assays$RNA@counts

# Get the number of samples to test
nums <- 1:length(unique(to_bulk_anno$Sample_cell))

# Parallelise making pseudobulk counts
pb_list <- mclapply(X = nums, FUN = get_pb, to_bulk_anno, aggr_col = "Sample_cell", mc.cores = 15)

# Bind the PB counts
bound <- bind_cols(pb_list)

bound[1:5,1:5]

saveRDS(bound, paste0(pseudobulk_dir, "bound_counts.rds"))

#bound <- readRDS( paste0(pseudobulk_dir, "bound_counts.rds"))

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
  distinct(Sample_cell,.keep_all = T)

bulk_anno <- data.frame(Sample_cell = colnames(bound_counts), check.names = F)%>%
  left_join(condensed_SC_anno)%>%
  mutate(Donor = gsub("_T|_N","", Sample_name))%>%
  mutate(Sample_type_cell = paste0(Sample_type, "_", Manual_toplevel_pred))%>%
  # Save the Pseudobulk annotation
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# bulk_anno <- read_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# Make a DGElist
rnaseq <- DGEList(bound_counts[,bulk_anno$Sample_cell])

# Remove weird characters
rownames(rnaseq) <- gsub("\\+", "_pos_", rownames(rnaseq))
colnames(rnaseq) <- gsub("\\+", "_pos_", colnames(rnaseq))
colnames(rnaseq) <- gsub("\\'|-|\\(|\\)", "_", colnames(rnaseq))
bulk_anno$Sample_cell <- gsub("\\+", "_pos_", bulk_anno$Sample_cell)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Sample_cell) == length(bulk_anno$Sample_cell)

# The experimental setup is matched
design <- model.matrix(~0 + Sample_type_cell + MMRStatus + Sex, data = bulk_anno)

# Neaten up design row and colnames
colnames(design) <- gsub("Sample_type_cell|\\+", "", colnames(design))
colnames(design) <- gsub("\\'|-|\\(|\\)", "_", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

#keep <- filterByExpr(study_counts, design = design)
keep <- rowSums(cpm(rnaseq)>1)>(0.05 * nrow(design))
table(keep)
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]

rnaseq <- calcNormFactors(rnaseq)

# No need to filter as all probes are expressed from single-cell
# Get the TMM normalised CPMs in case anyone wants to work with it
cpm <- edgeR::cpm(rnaseq, log =T)%>%
  data.frame()%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_log2_CPMs.csv"))

# Do a glimma MDS of batch removed counts
system(paste0("mkdir -p ", pseudobulk_dir,"glimma/mds/"))
#mds_save <- paste0(paste0(pseudobulk_dir,"glimma/", "MDS.html"))
#htmlwidgets::saveWidget(glimmaMDS(rnaseq, groups = bulk_anno, labels = bulk_anno$Sample_cell), mds_save)

# Save and MDS plot per cell type
for(celltype in unique(bulk_anno$Manual_toplevel_pred)){

  print(celltype)

  mds_save <-paste0(pseudobulk_dir,"glimma/mds/", celltype, "_MDS.html")

  # Filter downn to just one cell type
  rseq_filt <- rnaseq[,grepl(celltype, colnames(rnaseq))]

  if(ncol(rseq_filt) <3){
    next
  }

  htmlwidgets::saveWidget(glimmaMDS(rseq_filt, groups = bulk_anno[grepl(celltype, colnames(rnaseq)),],
                                    labels = bulk_anno$donor_sampletype_cell[grepl(celltype, colnames(rnaseq))]), mds_save)
}

# Normalise and fit linear model
#fit <- voomLmFit(rnaseq, design = design, plot = T)

# Voom for gsea
v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

saveRDS(fit, paste0(pseudobulk_dir, "fit.rds"))
#fit <- readRDS(paste0(pseudobulk_dir, "fit.rds"))

# Automate a contrast matrix
# Make sure baseline is the first group
unique(bulk_anno$Sample_type)

contrasts_manual <- "Granulo_high-Granulo_low"

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
  arrange(`Direction if significant`, `Number of genes`)%>%
  write_csv(paste0(pseudobulk_dir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  filter(`Number of genes` >2)%>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] *-1))%>%
  arrange(`Direction if significant`,-`Number of genes`)%>%
  mutate(Contrast = gsub("_", " ", Contrast))%>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

de_genes <- ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`))+
  geom_bar(stat = "identity")+
  labs(fill = "Direction")+
  blank_theme+
  theme(legend.key.size = unit(4, 'mm'))+
  scale_fill_manual(values = c("red", "blue"))

de_genes

# Make the output directory
system(paste0("mkdir -p ", pseudobulk_dir, "/toptables/"))

#system(paste0("rm ", pseudobulk_dir, "/toptables/*.csv"))

VOL <- paste0(pseudobulk_dir, "glimma/volcano/")
system(paste0("mkdir -p ", VOL))

# Get all the toptables
for(contrast in colnames(cont.matrix)){

  output <- paste0(pseudobulk_dir, "toptables/", contrast, ".csv")

  toptable <- topTable(fit.cont,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("SYMBOL")%>%
    write_csv(output)

  # Get back to just the cell type
  cont_first <- gsub("-.*", "", contrast)
  cont_first <- gsub("Granulo_high_", "", cont_first)

  rnaseq_filt <- rnaseq[,grepl(cont_first, colnames(rnaseq))]

  bulk_filt <- bulk_anno%>%
    filter(Sample_cell %in% colnames(rnaseq_filt))

  MA_fac <- factor(bulk_filt$Sample_type_cell, levels = unique(bulk_filt$Sample_type_cell))

  vol_save <- paste0(pseudobulk_dir, "glimma/volcano/",contrast,"_", "glimma_volcano.html")
  htmlwidgets::saveWidget(glimmaVolcano(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = round(rnaseq_filt$counts),
                                        dge = rnaseq_filt, groups = MA_fac), vol_save)

}

# Delete extra glimma files
system(paste0("rm -r ", pseudobulk_dir, "glimma/*/*_files"))

# Compile the toptables to compare with CPDB
all_toptables <- list.files(paste0(pseudobulk_dir, "toptables/"), full.names = T)

tt_list <- list()
for(i in 1:length(all_toptables)){

  contrast <- gsub(".csv", "", basename(all_toptables[i]))

  tt <- read_csv(all_toptables[i])%>%
    mutate(contrast = contrast)

  tt_list[[i]] <- tt


}

# Compile toptables and save the significant results
toptables_compiled <- bind_rows(tt_list)%>%
  mutate(cell_type = gsub(".*Granulo_low_", "", contrast))

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

# Gene set collections
gsea_gmt_dir <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
keep <- c(5,1,20, 28,31,9,16)
#keep <- c(31)
collections <- collections[keep]
# Add on my fetal signature
collections <- c(collections,"/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt")

# Make a directory
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/fry/"))

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)

  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix, mc.cores = 12,  pseudobulk_dir = pseudobulk_dir)

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
  mutate(cell_type = gsub(".*Granulo_low_", "", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

# Save some gene sets to to plot/to use as DE in other datasets
# Compile the camera results
all_fry <- list.files(paste0(pseudobulk_dir,"gsea/fry/"), full.names = T)

flist <- list()
for(i in 1:length(all_fry)){
  
  contrast <- gsub("\\.csv", "", basename(all_fry[i]))
  
  tt <- read_csv(all_fry[i], col_types = cols(.default = "c"))%>%
    mutate(contrast = contrast)%>%
    select(-Contrast)
  
  flist[[i]] <- tt
  
}

# Compile toptables and save the significant results
fry_compiled <- bind_rows(flist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  # Fix the cell type naming
  mutate(cell_type = gsub(".*Granulo_low_", "", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

fry_top_mono <- fry_compiled%>%
  filter(cell_type == "Mono")

# Look at some key gene sets
go_fry <- fry_compiled%>%
  filter(grepl("^GO", `Gene set`))

# Look at some key gene sets
go_fry <- fry_compiled%>%
  filter(grepl("^GO", `Gene set`))

# Look at some key gene sets
HALLMARK_fry <- fry_compiled%>%
  filter(grepl("^HALLMARK", `Gene set`))

HALLMARK <- camera_compiled%>%
  filter(grepl("^HALLMARK", `Gene set`))

stem <- camera_compiled%>%
  filter(grepl("*Enterocyte_", `Gene set`))%>%
  filter(cell_type == "Epi")

# Read in the hallmark gene sets for the monocyte and granulocyte populations
hallmark_mono <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/gsea/camera//h.all.v2023.1_Granulo_high_Mono-Granulo_low_Mono.csv")
hallmark_granulo <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/gsea/camera//h.all.v2023.1_Granulo_high_Granulo-Granulo_low_Granulo.csv")
hallmark_fibro <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/gsea/camera/h.all.v2023.1_Granulo_high_Fibro-Granulo_low_Fibro.csv")
hallmark_epi <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/gsea/camera/h.all.v2023.1_Granulo_high_Epi-Granulo_low_Epi.csv")

tt_mono <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/toptables/Granulo_high_Mono-Granulo_low_Mono.csv")%>%
  filter(adj.P.Val < 0.05)

tt_granulo <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/toptables/Granulo_high_Granulo-Granulo_low_Granulo.csv")%>%
  filter(adj.P.Val < 0.05)

tt_epi <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/toptables/Granulo_high_Epi-Granulo_low_Epi.csv")%>%
  filter(adj.P.Val < 0.05)

tt_fibro <- read_csv("./Results_paper/analyses/MMR_atlas_Neu/toptables/Granulo_high_Fibro-Granulo_low_Fibro.csv")%>%
  filter(adj.P.Val < 0.05)

# Make a folder to store the figure 4 inputs
fig_4_folder <- "./Results_paper/Extra_plots/Figure_4_inputs/"

# Read in figure 4a
four_a_1 <- "./Results_paper/Extra_plots/Image_inputs/4_a_1.png"
gc()
four_a_1 <- ggdraw() + 
  draw_image(
    four_a_1,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

# Read in figure 4a
four_a_2 <- "./Results_paper/Extra_plots/Image_inputs/4_a_2.png"
gc()
four_a_2 <- ggdraw() + 
  draw_image(
    four_a_2,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

four_a_3 <- "./Results_paper/Extra_plots/Image_inputs/4a_3.jpg"
gc()
four_a_3 <- ggdraw() + 
  draw_image(
    four_a_3,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

four_a_4 <- "./Results_paper/Extra_plots/Image_inputs/4a_4.jpg"
gc()
four_a_4 <- ggdraw() + 
  draw_image(
    four_a_4,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

# Make volcano plots of key cell types
fig_4b_i <- volcano_plot(tt_mono, "Granulocyte high monocytes")
  
fig_4b_i

fig_4b_ii <-  volcano_plot(tt_granulo, "Granulocyte high granulocytes")
  
fig_4b_ii

fig_4b_iii <-  volcano_plot(tt_fibro, "Granulocyte high fibroblasts")

fig_4b_iii

fig_4b_iv <-  volcano_plot(tt_epi, "Granulocyte high epi")

fig_4b_iv

four_b <- plot_grid(fig_4b_i,fig_4b_ii,fig_4b_iii,fig_4b_iv, labels = c("b"), 
                    label_size = 8, nrow = 2)

ggsave(plot = four_b,filename =  "./Results_paper/Extra_plots/Figure_4_inputs/Volcano.pdf",
       width = 6, height = 6)

# Look at the fetal gene signature
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), 
                       remove.empty=TRUE)

# Get the stem/foetal cell genes
stem_down <- gene_set_formatted$Adult_stem_vs_Adult_Enterocyte_down
stem_up <- gene_set_formatted$Adult_stem_vs_Adult_Enterocyte_up
fetal_down <- gene_set_formatted$Second_trim_Enterocyte_vs_Adult_Enterocyte_down
fetal_up <- gene_set_formatted$Second_trim_Enterocyte_vs_Adult_Enterocyte_up

tt_epi_stem_down <- tt_epi%>%
  filter(SYMBOL %in% stem_down)%>%
  filter(logFC <0)%>%
  mutate(`Gene set` = "stem_down")

tt_epi_fetal_up <- tt_epi%>%
  filter(SYMBOL %in% fetal_up)%>%
  filter(logFC >0)%>%
  mutate(`Gene set` = "fetal_up")

tt_epi_fetal_down <- tt_epi%>%
  filter(SYMBOL %in% fetal_down)%>%
  filter(logFC <0)%>%
  mutate(`Gene set` = "fetal_down")

tt_epi_stem_up <- tt_epi%>%
  filter(SYMBOL %in% stem_up)%>%
  filter(logFC >0)%>%
  mutate(`Gene set` = "stem_up")

# Find genes that might also be in the spatial
merged.integrated <- qread("./Intermediate/lognorm_merged_integrated_annotated.qs")

combined_stem <- list(tt_epi_stem_down, tt_epi_fetal_up, tt_epi_fetal_down, tt_epi_stem_up)%>%
  bind_rows()%>%
  filter(!duplicated(SYMBOL))%>%
  arrange(adj.P.Val)%>%
  mutate(is_in_CosMx = SYMBOL %in% rownames(merged.integrated))%>%
  write_csv("./Results_paper/analyses/MMR_atlas_Neu/Genes_to_validate.csv")

granulo_tab <- read_csv("./Results_paper/Tables/MMR atlas granulo pcts.csv")%>%
  filter(Pct_granulo >2)

mmr_atlas$granulo_high <- ifelse(mmr_atlas$Sample_name %in% granulo_tab$Sample_name,"High", "Low")

mmr_atlas_epi <- mmr_atlas[,mmr_atlas$clMidwayPr == "Epi"]

DotPlot(mmr_atlas_epi,combined_stem$SYMBOL, group.by = "granulo_high")

fig_4_c <- function() {
  par(
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0),
    cex=0.25
  )
  barcodeplot(fit.cont$t[,"Granulo_high_Epi-Granulo_low_Epi"],
              index=indexed$Adult_stem_vs_Adult_Enterocyte_down,
              labels = c("Down","Up"),
              main="Stem vs enterocyte down", 
              xlab = "Limma t statistic (granulocyte high vs low tumour cells)")
  
}

fig_4_c <- ggdraw(fig_4_c)+
  theme_bw(base_size = 2)

# Read the T vs N results as well
t_vs_n_epi <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/toptables/T_Epi-N_Epi.csv")
indexed_2 <- ids2indices(gene.sets = gene_set_formatted, identifiers =t_vs_n_epi$SYMBOL, remove.empty=TRUE)

tt_epi_stem_t_n <- t_vs_n_epi%>%
  filter(SYMBOL %in% stem)

fig_4_c_2 <- function() {
  par(
    mar = c(3, 3, 1, 1),
    mgp = c(2, 1, 0),
    cex=0.25
  )
  barcodeplot(t_vs_n_epi$t,
                           index=indexed_2$Adult_stem_vs_Adult_Enterocyte_down,
                           labels = c("Down","Up"),
                           main="Stem vs enterocyte down", 
                           xlab = "Limma t statistic (tumour vs normal epithelial cells)")
  
}

pdf("./Results_paper/Extra_plots/Figure_4_inputs/Tumour normal barcode.pdf", width = 7, height = 3)
barcodeplot(t_vs_n_epi$t,
            index=indexed_2$Adult_stem_vs_Adult_Enterocyte_down,
            labels = c("Down","Up"),
            main="Stem vs enterocyte down", 
            xlab = "Limma t statistic (tumour vs normal epithelial cells)")
dev.off()

pdf("./Results_paper/Extra_plots/Figure_4_inputs/Granulo high low barcode.pdf", width = 7, height = 3)
barcodeplot(fit.cont$t[,"Granulo_high_Epi-Granulo_low_Epi"],
            index=indexed$Adult_stem_vs_Adult_Enterocyte_down,
            labels = c("Down","Up"),
            main="Stem vs enterocyte down", 
            xlab = "Limma t statistic (granulocyte high vs low tumour cells)")
dev.off()

fig_4_c_2 <- ggdraw(fig_4_c_2)+
  theme_bw(base_size = 2)

four_a <- plot_grid(four_a_1,four_a_3, four_a_2,four_a_4, labels = c("a"), 
                    label_size = 8, nrow = 2, rel_widths = c(1.2, 1))

four_c <- plot_grid(fig_4_c, fig_4_c_2, labels = c("c"), 
                    label_size = 8, nrow = 2)

four_e_1 <- "./Results_paper/Extra_plots/Image_inputs/4e1.png"
gc()
four_e_1 <- ggdraw() + 
  draw_image(
    four_e_1,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

four_e_2 <- "./Results_paper/Extra_plots/Image_inputs/4e_2.png"
gc()
four_e_2 <- ggdraw() + 
  draw_image(
    four_e_2,scale = 1, width = 1.5, halign = 0, valign = 0)+
  blank_theme+
  theme(panel.border = element_blank())

four_e <- plot_grid(four_e_1, four_e_2, labels = c("e"), 
                    label_size = 8, nrow = 2)

# Make the
four_d_1 <- gsea_plot(hallmark_mono, 10, "Granulocyte rich monocytes")
ggsave(filename = "./Results_paper/Extra_plots/Figure_4_inputs/Granulocyte rich monocytes hallmark.pdf", 
       plot = four_d_1,
       width = 7, height = 3)
four_d_2 <- gsea_plot(hallmark_granulo, 10, "Granulocyte rich granulocytes")
ggsave(filename = "./Results_paper/Extra_plots/Figure_4_inputs/Granulocyte rich granulocytes hallmark.pdf", 
       plot = four_d_2,
       width = 7, height = 3)
four_d_3 <- gsea_plot(hallmark_fibro, 10, "Granulocyte rich fibroblasts")
ggsave(filename = "./Results_paper/Extra_plots/Figure_4_inputs/Granulocyte rich fibroblasts hallmark.pdf", 
       plot = four_d_3,
       width = 7, height = 3)

four_d <- plot_grid(four_d_1, four_d_2, four_d_3, four_c)

top <- plot_grid(four_a, four_b)
bottom <- plot_grid( four_d, four_e, nrow = 1, rel_widths = c(2, 1))

# What stats are you putting on old vs young
Figure_4 <- plot_grid(top, bottom, nrow = 2)

# Plot CXCL14 in fibros and PPBP in neus
plot_genes <- cpm%>%
  filter(Gene %in% c("CXCL14", "PPBP"))%>%
  gather(Sample_cell, CPM, -Gene)%>%
  left_join(bulk_anno)%>%
  filter(Manual_toplevel_pred %in% c("Fibro", "Macro", "Mono", "Epi"))%>%
  mutate(Sample_type = gsub("_", " ", Sample_type))%>%
  mutate(Sample_type = factor(Sample_type, levels = c("Granulo low", "Granulo high")))

# https://ggplot2.tidyverse.org/reference/position_jitterdodge.html
ggplot(data = plot_genes, aes(x = Sample_type, y = CPM, fill = Manual_toplevel_pred))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0),
             aes(colour = Manual_toplevel_pred))+
  geom_boxplot(outlier.alpha = 0, alpha = 0.6)+
  facet_wrap(~Gene)+
  blank_theme+
  scale_color_manual(values = cell_type_colors)+
  scale_fill_manual(values = cell_type_colors)+
  guides(fill = "none")+
  labs(x = "Granulocyte infiltration", y = expression('Log'[2]*' CPM'),
       colour = "Cell type")

ggsave("./Results_paper/Extra_plots/Figure_4_inputs/Gene expression boxplots.pdf", width = 5, height = 2)

ggsave(plot = Figure_4,"./Results_paper/Plots/Figure 4.pdf", 
       width = 170, height = 170, units = "mm")

bulk_anno$Sample_type <- factor(bulk_anno$Sample_type, levels = c("Granulo_low", "Granulo_high"))

# Try looking at just the pct granulo in the Epi cells
design <- model.matrix(~0 + Manual_toplevel_pred + Manual_toplevel_pred:Sample_type + MMRStatus+ Sex + TumorStage, data = bulk_anno)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Sample_cell) == length(bulk_anno$Sample_cell)

v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

fit <- eBayes(fit)
summa.fit <- decideTests(fit)
summary(summa.fit)

tt <- topTable(fit,coef="Manual_toplevel_predEpi:Sample_typeGranulo_high",number = Inf)%>%
  rownames_to_column("SYMBOL")

# Try again using granulo % as a numeric rather than as a factor
design <- model.matrix(~0 + Manual_toplevel_pred + Manual_toplevel_pred:Pct_granulo + MMRStatus+ Sex + TumorStage, data = bulk_anno)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Sample_cell) == length(bulk_anno$Sample_cell)

v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
summa.fit <- decideTests(fit)
summary(summa.fit)

tt <- topTable(fit,coef="Manual_toplevel_predMono:Pct_granulo",number = Inf)%>%
  rownames_to_column("SYMBOL")

tt_epi <- topTable(fit,coef="Manual_toplevel_predEpi:Pct_granulo",number = Inf)%>%
  rownames_to_column("SYMBOL")

tt_fibro <- topTable(fit,coef="Manual_toplevel_predFibro:Pct_granulo",number = Inf)%>%
  rownames_to_column("SYMBOL")



# Look at the fetal gene signature
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), 
                       remove.empty=TRUE)

camera_result <- camera(y = v ,index = indexed, design = design, contrast = "Manual_toplevel_predEpi:Pct_granulo")%>%
  rownames_to_column("Gene set")%>%
  dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")

fry_result <- fry(y = v ,index = indexed, design = design, contrast = "Manual_toplevel_predEpi:Pct_granulo")%>%
  rownames_to_column("Gene set")%>%
  dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
  filter(FDR <= 0.05)

# Look at the fetal gene signature
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), 
                       remove.empty=TRUE)

camera_result <- camera(y = v ,index = indexed, design = design, contrast = "Manual_toplevel_predEpi:Pct_granulo")%>%
  rownames_to_column("Gene set")%>%
  dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")

fry_result <- fry(y = v ,index = indexed, design = design, contrast = "Manual_toplevel_predEpi:Pct_granulo")%>%
  rownames_to_column("Gene set")%>%
  dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")


# Look at the fetal gene signature
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/c5.all.v2023.1.Hs.symbols.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), 
                       remove.empty=TRUE)

camera_result <- camera(y = v ,index = indexed, design = design, contrast = "Manual_toplevel_predEpi:Pct_granulo")%>%
  rownames_to_column("Gene set")%>%
  dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")

# Plot some key genes
cpm <- edgeR::cpm(rnaseq, log =T)

# Plot the expression of some genes vs granulo %s
TGFB1 <- cpm["TGFB1",]%>%
  data.frame()%>%
  rownames_to_column("Sample_cell")%>%
  dplyr::rename(TGFB1 = 2)%>%
  left_join(bulk_anno)%>%
  filter(Manual_toplevel_pred == "Epi")%>%
  filter(Pct_granulo > 0.5)%>%
  mutate(Log_pct_granulo = log(Pct_granulo +1))

cor.test(TGFB1$Pct_granulo, TGFB1$TGFB1)

ggplot(data = TGFB1, aes(x = Log_pct_granulo, y  = TGFB1))+ geom_point()

PDE10A <- cpm["PDE10A",]%>%
  data.frame()%>%
  rownames_to_column("Sample_cell")%>%
  dplyr::rename(PDE10A = 2)%>%
  left_join(bulk_anno)%>%
  filter(Manual_toplevel_pred == "Epi")%>%
  filter(Pct_granulo > 0)%>%
  mutate(Log_pct_granulo = log(Pct_granulo +1))

cor.test(PDE10A$Pct_granulo, PDE10A$PDE10A)

ggplot(data = PDE10A, aes(x = Log_pct_granulo, y  = PDE10A))+ geom_point()

plot_gene_celltype (cpm, "PIGR", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL14", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "OLFM4", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "TGFB1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CD33", bulk_anno, paste0(pseudobulk_dir,"Plots/"))

# Look at the data from the perspective of high neutrophil T vs T vs N vs N ----
CPM <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Pseudobulk_log2_CPMs.csv")
IL24 <- CPM%>%
  filter(Gene == "IL24")%>%
  gather(Sample_cell, CPM, -Gene)
anno <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Pseudobulk_annotation.csv")%>%
  left_join(granulo_pcts)%>%
  left_join(IL24)%>%
  mutate(Cell_type = gsub(".*_", "", Sample_type_cell))%>%
  mutate(Sample_type_g = ifelse(Sample_name %in% granulo_keep, "Granulo_high", "Granulo_low"))%>%
  mutate(Sample_type_g = paste0(Sample_type_g, " ", Sample_type))

ggplot(data = anno, aes(x = log2(Pct_granulo), y = CPM, colour = Sample_type))+ geom_jitter()+
  facet_wrap(~Cell_type)

# Plot the fibroiblasts only
only_fibro <- anno%>%
  filter(Cell_type == "Fibro")%>%
  mutate(Sample_type_g = replace(Sample_type_g, Sample_type_g == "Granulo_low T","Granulo low tumour"))%>%
  mutate(Sample_type_g = replace(Sample_type_g, Sample_type_g == "Granulo_low N","Normal"))%>%
  mutate(Sample_type_g = replace(Sample_type_g, Sample_type_g == "Granulo_high T","Granulo high tumour"))%>%
  mutate(Sample_type_g = factor(Sample_type_g, levels = c("Normal", "Granulo low tumour", "Granulo high tumour")))

cor.test(only_fibro$Pct_granulo, only_fibro$CPM)

ggplot(data = only_fibro, aes(x = log2(Pct_granulo), y = CPM, colour = Sample_type))+ geom_point()+
  facet_wrap(~Sample_type)

p1 <- ggplot(data = only_fibro, aes(x = Sample_type_g, y = CPM))+ 
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter()+
  blank_theme+
  labs(x = "Granulocyte infiltration", y = substitute( paste(italic('IL24'), ' log'[2]*' CPM' )))

p1

ggsave(plot = p1, filename = "./Results_paper/analyses/MMR_atlas_Neu/Plots/IL24 granulo.pdf", )

# Make a DGElist
rnaseq <- DGEList(bound_counts[,bulk_anno$Sample_cell])

# Do a DE analysis of the fibros with granulo pct as a numeric
# The experimental setup is matched
design <- model.matrix(~Pct_granulo:Manual_toplevel_pred + Manual_toplevel_pred + MMRStatus + Sex, data = bulk_anno)

# Neaten up design row and colnames
rownames(design) <- rownames(rnaseq$samples)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Sample_cell) == length(bulk_anno$Sample_cell)

#keep <- filterByExpr(study_counts, design = design)
keep <- rowSums(cpm(rnaseq)>1)>(0.05 * nrow(design))
table(keep)
rnaseq <- rnaseq[keep,, keep.lib.sizes=FALSE]

rnaseq <- calcNormFactors(rnaseq)

v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

fit <- eBayes(fit)
summa.fit <- decideTests(fit)
summary(summa.fit)

fibro <- topTable(fit,coef = "Pct_granulo:Manual_toplevel_predFibro", n= Inf)%>%
 # filter(adj.P.Val < 0.05)%>%
  rownames_to_column("Gene")

macro <- topTable(fit,coef = "Pct_granulo:Manual_toplevel_predMacro", n= Inf)%>%
  rownames_to_column("Gene")

Epi <- topTable(fit,coef = "Pct_granulo:Manual_toplevel_predEpi", n= Inf)%>%
  rownames_to_column("Gene")

TCD8 <- topTable(fit,coef = "Pct_granulo:Manual_toplevel_predTCD8", n= Inf)%>%
  rownames_to_column("Gene")

granulo <- topTable(fit,coef = "Pct_granulo:Manual_toplevel_predPlasma", n= Inf)%>%
  rownames_to_column("Gene")

anno <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Pseudobulk_annotation.csv")%>%
  left_join(granulo_pcts)%>%
  mutate(Cell_type = gsub(".*_", "", Sample_type_cell))%>%
  mutate(Sample_type_g = ifelse(Sample_name %in% granulo_keep, "Granulo_high", "Granulo_low"))%>%
  mutate(Sample_type_g = paste0(Sample_type_g, " ", Sample_type))%>%
  mutate(Sample_type_g = replace(Sample_type_g, Sample_type_g == "Granulo_low T","Granulo low tumour"))%>%
  mutate(Sample_type_g = replace(Sample_type_g, Sample_type_g == "Granulo_low N","Normal"))%>%
  mutate(Sample_type_g = replace(Sample_type_g, Sample_type_g == "Granulo_high T","Granulo high tumour"))%>%
  mutate(Sample_type_g = factor(Sample_type_g, levels = c("Normal", "Granulo low tumour", "Granulo high tumour")))

bulk_anno_epi <- anno%>%
  filter(Manual_toplevel_pred == "Epi")

# Try a GSVA analysis for the MMR atlas
CPM_mat <- as.matrix(CPM[,2:ncol(CPM)])
rownames(CPM_mat) <- CPM$Gene
ordered <- CPM_mat[,bulk_anno_epi$Sample_cell]
ordered[1:5,1:5]

revsc_gene_set <- c(
  "C9orf116", "SMIM33", "C6orf141", "ABCB1", "ABI3", "ACSM3", "AKR1B10", "AKR1B1", "AKR1B10", 
  "ALCAM", "ANXA1", "ANXA3", "ANXA5", "AP1S3", "ARHGAP44", "ARHGEF4", "ARL14", "ARL4C", "AVPI1", 
  "BAIAP2L2", "BASP1", "BOK", "CAMKK1", "CAMSAP2", "CAP2", "CAPN5", "CA12", "CA2", "CAVIN3", 
  "CCDC92", "CCL2", "CCN2", "CD14", "CD44", "CD68", "CELF2", "CGNL1", "CKAP2", "CKB", "CLCA4", 
  "CLDN1", "CLDN4", "CLIC3", "CLU", "CNN2", "CRYBG1", "CTSH", "CWH43", "CXCL1", "CXCL10", 
  "CXCL16", "CXCL2", "CXCL5", "CXCL9", "DAPK1", "DDAH1", "DYRK3", "EMP1", "EMP2", "ENAH", 
  "ENDOD1", "EPHA2", "EPHA4", "EPN3", "EPOP", "EPS8L1", "ERCC5", "EREG", "ETV5", "EYA2", 
  "FABP5", "FFAR4", "FHL2", "FKBP9", "FLNA", "FMNL2", "FOSL1", "FUT2", "FXYD3", "GADD45B", 
  "GASK1B", "GATM", "GCNT1", "GCNT3", "GJB3", "GLIPR1", "GPNMB", "GPR137B", "GPX2", "GTSE1", 
  "HCAR2", "HMOX1", "HSPA1B", "HSPA2", "HSPB1", "HSPH1", "HYAL1", "ICA1", "ICAM1", "ID2", 
  "IFITM3", "IL17RE", "IL1RN", "IL33", "ILDR1", "INKA2", "IRF1", "ITGA2", "ITGAV", "ITGB4", 
  "ITPRIPL2", "KCNN4", "KIFC3", "KLF7", "KRT15", "KRT18", "KRT23", "LAMA3", "LAPTM5", "LCN2", 
  "LDLR", "LGALS1", "LPCAT4", "LRG1", "MDM2", "MFSD4A", "MISP3", "MRAS", "MSANTD3", "MSLN", 
  "MSN", "MST1R", "MTCL1", "MUC3A", "MUC4", "MYH9", "NFE2L3", "NFKBIA", "NFKBIZ", "NOS2", 
  "NOXO1", "OAS1", "OSBPL5", "P2RY2", "P2RY6", "PAQR8", "PDZK1IP1", "PHLDA1", "PHLDA2", 
  "PIMREG", "PLAU", "PLAUR", "PLCD3", "PLEKHS1", "PLK2", "PMEPA1", "PPIC", "PPL", "PSRC1", 
  "PTGS2", "RAB11FIP5", "RAB15", "RACK1", "RALGAPA2", "RBMS1", "RGS19", "RHOD", "RIPK2", 
  "RPPH1", "RUNX1", "S100A14", "S100A4", "S100A6", "SAMD5", "SCARNA17", "SCEL", "SEMA7A", 
  "SERPINB9", "SH3BGRL2", "SLC19A2", "SLC39A10", "SLC44A2", "SLC7A5", "SLCO4A1", "SLPI", 
  "SNORA20", "SNORA73A", "SNORA73B", "SNORD104", "SNORD118", "SNORD43", "SNORD92", "SOCS3", 
  "SORL1", "SPATA6", "SPP1", "STX11", "SULF2", "SYT8", "TCIM", "TCIRG1", "TGFBI", "TGM2", 
  "TIMP2", "TINAGL1", "TLR4", "TMEM184C", "TMEM43", "TNF", "TNFAIP3", "TNFAIP8", "TNFRSF1B", 
  "TNFSF9", "TNIP3", "TNNI2", "TNNT2", "TNS4", "TRAM2", "TRIM15", "TRIM40", "TSPAN4", "TSPAN6", 
  "TUBA1A", "VAMP5", "VIM", "VNN1", "WFDC2", "WWTR1", "ZNF296", "ZNF703", "ZMAT3", "HLA-DQB1", 
  "HLA-G"
)

# Look at the fetal gene signature
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names

gs <- list("revsc_gene_set" = revsc_gene_set)

gene_set_formatted <- c(gene_set_formatted, gs)
str(gene_set_formatted)

# Run GSVA across the Epis
param <- gsvaParam(exprData = ordered, gene_set_formatted)
gsva.es <- gsva(param, verbose=FALSE)

scores <- data.frame(gsva.es)%>%
  rownames_to_column("Gene_set")%>%
  gather(Sample_cell, GSVA_score, -Gene_set)%>%
  left_join(bulk_anno_epi)

ggplot(data = scores, aes(x = log2(Pct_granulo+0.1), y = GSVA_score, colour = Sample_type_g))+ 
  geom_point()+
  facet_wrap(~Gene_set)

ggplot(data = scores, aes(x = Sample_type_g, y = GSVA_score))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(height = 0, width = 0.3)+
  labs(x = "MMR status", y = "GSVA_score")+
  blank_theme+
  facet_wrap(~Gene_set)
