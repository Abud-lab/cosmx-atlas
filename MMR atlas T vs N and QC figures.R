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
library(cowplot)
library(grid)
library(gridExtra)
library(ggrepel)

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/functions.R")

# Set a homedir for the analysis so as to no use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

# Read in the MMR atlas with logtransformed counts
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

pseudobulk_dir <- "./Results_paper/analyses/MMR_atlas_T_vs_N/"

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

# Find the markers genes for each cell type
Idents(mmr_atlas) <- mmr_atlas$clMidwayPr
marks <- FindAllMarkers(mmr_atlas) 

write_csv(marks, paste0(pseudobulk_dir, "Celltype_markers.csv"))

granulo <- FindMarkers(mmr_atlas,ident.1 = "Granulo")%>%
  rownames_to_column("Gene")%>%
  write_csv(paste0(pseudobulk_dir, "Celltype_markers.csv"))

# Take a loot at the MMR atlas and see if you can do any better
# with processing
Idents(mmr_atlas) <- mmr_atlas$PatientTypeID
VlnPlot(mmr_atlas, features = c("nCount_RNA"),pt.size = -1)+NoLegend()

mmr_atlas$Sample_name <- mmr_atlas$PatientTypeID

mmr_atlas@assays$RNA@data[1:20,1:20]

# Figure out differential niche usage between conditions
table(mmr_atlas$clMidwayPr)
table(mmr_atlas$Sex)
table(mmr_atlas$MMRStatus)
table(mmr_atlas$SPECIMEN_TYPE)
mmr_atlas$Manual_toplevel_pred <- mmr_atlas$clMidwayPr

# Study annotations
unique(mmr_atlas$clMidwayPr)

to_bulk_anno <- mmr_atlas@meta.data%>%
  dplyr::select(Barcode = sampleID, Manual_toplevel_pred, Sample_name = PatientTypeID,
                Sample_type = SPECIMEN_TYPE, MMRStatus, Sex)%>%
  mutate(T_N = gsub(".*_", "", Sample_name))%>%
  mutate(Sample_cell = paste0(Sample_name, "_", Manual_toplevel_pred))

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

# Keep only the T samples for now
#mmr_atlas <- mmr_atlas[,mmr_atlas$PatientTypeID %in% to_bulk_anno$Sample_name]
mmr_atlas$type <-mmr_atlas$SPECIMEN_TYPE
table(mmr_atlas$type)

# Run sccomp with contrasts for tumour vs normal
sc_result <- mmr_atlas |>
  sccomp_glm(
    formula_composition = ~type,
    .sample = PatientTypeID,
    .cell_group = Manual_toplevel_pred,
    bimodal_mean_variability_association = T,
    cores = 10
  )

plots <- plot_summary(sc_result)

CT_DA <- plots$boxplot[[1]]
CT_DA <- CT_DA+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4))

pdf(paste0(pseudobulk_dir, "Plots/sccomp DA boxplot.pdf"),
    width = 7, height = 4.5)
CT_DA
dev.off()

# Get the raw counts
pb_counts <- mmr_atlas@assays$RNA@counts

# Remove the single cell object save RAM
rm(mmr_atlas)

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

# Figure out the samples I want to keep
sample_keep <- condensed_SC_anno%>%
  filter(Sample_type %in% c("T", "N"))%>%
  mutate(Donor = gsub("_T|_N","", Sample_name))%>%
  dplyr::select(Sample_type, Sample_name,Donor)%>%
  distinct()%>%
  group_by(Donor)%>%
  mutate(Count = n())%>%
  filter(Count == 2)

bulk_anno <- data.frame(Sample_cell = colnames(bound_counts), check.names = F)%>%
  left_join(condensed_SC_anno)%>%
  mutate(Donor = gsub("_T|_N","", Sample_name))%>%
  mutate(Sample_type_cell = paste0(Sample_type, "_", Manual_toplevel_pred))%>%
  # Drop samples that aren't matched
  filter(Sample_name %in% sample_keep$Sample_name)%>%
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
design <- model.matrix(~0 + Sample_type_cell + Donor, data = bulk_anno)

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

cpm[1:5,1:5]

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
  arrange(`Direction if significant`, `Number of genes`)%>%
  write_csv(paste0(pseudobulk_dir, "Significant_genes_summary.csv"))

plot_summary <- de_summary %>%
  filter(`Direction if significant`!= "NotSig")%>%
  filter(`Number of genes` >10)%>%
  mutate(`Number of genes` = replace(`Number of genes`, `Direction if significant` == "Down", `Number of genes`[`Direction if significant` == "Down"] *-1))%>%
  arrange(`Direction if significant`,-`Number of genes`)%>%
  mutate(Contrast = factor(Contrast, levels = unique(Contrast)))

ggplot(data = plot_summary, aes(y = Contrast, x = `Number of genes`, fill = `Direction if significant`))+
  geom_bar(stat = "identity")

# Make the output directory
system(paste0("mkdir -p ", pseudobulk_dir, "/toptables/"))

system(paste0("rm ", pseudobulk_dir, "/toptables/*.csv"))

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
  cont_first <- gsub("T_", "", cont_first)

  cpm_filt <- cpm[,grepl(cont_first, colnames(cpm))]

  bulk_filt <- bulk_anno%>%
    filter(Sample_cell %in% colnames(cpm_filt))

  MA_fac <- factor(bulk_filt$Sample_type_cell, levels = unique(bulk_filt$Sample_type_cell))

  vol_save <- paste0(pseudobulk_dir, "glimma/volcano/",contrast,"_", "glimma_volcano.html")
  htmlwidgets::saveWidget(glimmaVolcano(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = cpm_filt,transform.counts = "none",
                                        groups = MA_fac), vol_save)

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
  mutate(cell_type = gsub(".*N_", "", contrast))

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

KRT23 <- toptables_compiled%>%
  filter(SYMBOL == "KRT23")

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

# Run the epi gene set for the epi cells
all_epi_gene_sets <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Gut_cell_atlas_epithelial_cell_types.gmt"

epi_reslt <- run_GSEA(contrast = "T_Epi-N_Epi", collection = all_epi_gene_sets,rnaseq =  rnaseq,v= v, design = design, 
         cont.matrix,pseudobulk_dir = pseudobulk_dir)

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)

  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix,pseudobulk_dir = pseudobulk_dir, mc.cores = 7)

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

# Make some barcode plots for the top gut cell atlas gene sets
# Read in the fetal gene set
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Gut_cell_atlas_epithelial_cell_types.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), remove.empty=TRUE)

t_vs_n_epi <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/toptables/T_Epi-N_Epi.csv")

t_vs_n_epi_IBD_stem_ped <- t_vs_n_epi%>%
  filter(adj.P.Val < 0.05)%>%
  filter(SYMBOL %in% gene_set_formatted$Pediatric_IBD_Stem_cells_up)

# Make some barcode plots of some key gene sets
pdf(paste0(pseudobulk_dir, "Plots/", "T_Epi_Pediatric_IBD_Enterocyte_up.pdf"), width = 7, height = 6)
barcodeplot(fit.cont$t[,"T_Epi-N_Epi"],
            index=indexed$Pediatric_IBD_Enterocyte_up,
            index2 = indexed$Pediatric_IBD_Enterocyte_up,
            labels = c("Down in tumours","Up in tumours"),
            main="Pediatric IBD Enterocyte up")
dev.off()

# Make some barcode plots of some key gene sets
pdf(paste0(pseudobulk_dir, "Plots/", "Pediatric_IBD_Stem_cells_up.pdf"), width = 7, height = 6)
barcodeplot(fit.cont$t[,"T_Epi-N_Epi"],
            index=indexed$Pediatric_IBD_Stem_cells_up,
            index2 = indexed$Pediatric_IBD_Stem_cells_up,
            labels = c("Down in tumours","Up in tumours"),
            main="Pediatric IBD Stem cells up")
dev.off()

# Read in the fetal gene set
collection <- "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt"
gene_set <- quiet(GSA.read.gmt(collection))
gene_set_formatted <- gene_set$genesets
names(gene_set_formatted) <- gene_set$geneset.names
indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), remove.empty=TRUE)

# Make some barcode plots of some key gene sets
pdf(paste0(pseudobulk_dir, "Plots/", "Tumour vs normal epithelial cells fetal genes barcodeplot.pdf"), width = 7, height = 6)
barcodeplot(fit.cont$t[,"T_Epi-N_Epi"],
            index=indexed$Second_trim_Enterocyte_vs_Adult_Enterocyte_up,
            index2 = indexed$Second_trim_Enterocyte_vs_Adult_Enterocyte_down,
            labels = c("Down in fetal enterocyte","Up in fetal enterocyte"),
            main="Tumour vs normal epithelial cells")
dev.off()

# Make some barcode plots of some key gene sets
pdf(paste0(pseudobulk_dir, "Plots/", "Tumour vs normal epithelial cells fetal genes barcodeplot.pdf"), width = 7, height = 6)
barcodeplot(fit.cont$t[,"T_Epi-N_Epi"],
            index=indexed$Second_trim_Enterocyte_vs_Adult_Enterocyte_up,
            index2 = indexed$Second_trim_Enterocyte_vs_Adult_Enterocyte_down,
            labels = c("Down in fetal enterocyte","Up in fetal enterocyte"),
            main="Tumour vs normal epithelial cells")
dev.off()

pdf(paste0(pseudobulk_dir, "Plots/", "Tumour vs normal epithelial cells stem cell genes barcodeplot.pdf"), width = 7, height = 6)
barcodeplot(fit.cont$t[,"T_Epi-N_Epi"],
            index=indexed$Adult_stem_vs_Adult_Enterocyte_up,
            index2 = indexed$Adult_stem_vs_Adult_Enterocyte_down,
            labels = c("Down in stem","Up in stem"),
            main="Tumour vs normal epithelial cells")
dev.off()

# Plot the expression of a given gene in a given celltype
cpm <- edgeR::cpm(rnaseq, log =T)

unique(bulk_anno$Sample_type)
bulk_mmr <- bulk_anno%>%
  mutate(MMRStatus = replace(MMRStatus, is.na(MMRStatus), ""))%>%
  mutate(Sample_type = paste0(Sample_type, " ", MMRStatus))

# t_vs_n_gene (cpm, "SPP1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCL8", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCL5", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCL1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCL2", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCL3", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCR2", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "IL23A", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
# t_vs_n_gene (cpm, "CXCL13", bulk_anno, paste0(pseudobulk_dir,"Plots/"))

plot_gene_celltype (cpm, "CSPG4", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "FZD1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "NECTIN3", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))


plot_gene_celltype (cpm, "CLU", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CLU", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "FCGBP", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "FXYD5", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL14", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL13", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL11", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL6", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL24", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL11RA", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "BEST4", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "OTOP2", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "LEFTY1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL1RN", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "NRG1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "PDGFRA", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL24", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "WNT10B", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "KRT23", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CCL2", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "TIMP1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "DLL1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "GREM1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "BMP4", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "BMP7", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL23A", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "TNF", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL1B", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "IL1A", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCR2", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL8", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL2", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL3", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL5", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "SPP1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CEACAM8", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CRP", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "MMP7", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CCL18", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
splot_gene_celltype (cpm, "TIMP1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "PLAC8", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "ITGAX", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "SLC2A1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL13", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "REG1A", bulk_anno, paste0(pseudobulk_dir,"Plots/"))


plot_gene_celltype (cpm, "MRC1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCR2", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL1", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL2", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL3", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL5", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))

# Show the CXCL8/12 axis
c8 <- plot_gene_celltype (cpm, "CXCL8", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
c12 <- plot_gene_celltype (cpm, "CXCL12", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
c5 <- plot_gene_celltype (cpm, "CXCL5", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))
fap <- plot_gene_celltype (cpm, "FAP", bulk_mmr, paste0(pseudobulk_dir,"Plots/"))

genes <- c("CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CXCL12", 
           "CXCL14", "PPBP", "ACKR1", "SPP1", "MMP1", "MMP3", "IL6", "BMP4", "BMP5")
celltypes <- c("Fibro", "Endo", "Macro", "Mono")

# Add a granulocyte % barplot to the heatmap
granulo_tab <- read_csv("./Results_paper/Tables/MMR atlas granulo pcts.csv")

bulk_mmr_granulo <- bulk_mmr%>%
  select(-MMRStatus)%>%
  left_join(granulo_tab)

# Plot a heatmap of CXCL genes in niches and in T vs N. 
cxcl_hm <- plot_hm_mmr(genes = genes,bulk_anno = bulk_mmr_granulo,celltypes = celltypes,cpm = cpm)
cxcl_hm

# Read in the CPM and annotation from the niche level counts
niche_cpm <- qread("./Results_paper/analyses/Pseudobulk_niches_HQ/Pseudobulk_log2_CPMs.qs")
niche_cpm[1:5,1:5]

niche_anno <- read_csv("./Results_paper/analyses/Pseudobulk_niches_HQ/Pseudobulk_annotation.csv")

niches <- c("Granulocyte_rich", "TLS", "Epithelial_mass")

niche_anno <- niche_anno%>%
  filter(Sample == "118T" & niches == "Granulocyte_rich" | Sample == "73T" & niches == "TLS"| niches == "Epithelial_mass")%>%
  filter(Slide !="Run5850.2186.1")

cxcl_hm_niches <- plot_hm_niches(genes = genes, bulk_anno = niche_anno,
                                 cpm = niche_cpm,celltypes = celltypes,cellcols = cell_type_colors, 
                                 nichecols = nichecols_underscore, samplecols = Sample_cols,
                                 niches_keep = niches,
                                 donor_keep  = c("118", "73"))
cxcl_hm_niches

gb = grid.grabExpr(draw(cxcl_hm))
g2 = grid.grabExpr(draw(cxcl_hm_niches))

bottom <- plot_grid(c8,c12, labels = c("c", "d"), label_size = 8)

combined_cxcl <- plot_grid(gb, g2,bottom, nrow = 3, 
                                labels = c("a", "b"), 
                                label_size = 8)

combined_cxcl

ggsave(plot = combined_cxcl,"./Results_paper/Plots/Figure S granulo.pdf", 
       width = 170, height = 170, units = "mm")

# Plots some CD44 associated genes
plot_gene_celltype (cpm, "CD44", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
plot_gene_celltype (cpm, "HGF", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
plot_gene_celltype (cpm, "HLA-DRA", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)


# Plot some key tumour genes
t1 <- plot_gene_celltype (cpm, "EPCAM", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
t2 <- plot_gene_celltype (cpm, "TGFBI", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
t3 <- plot_gene_celltype (cpm, "FXYD5", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
t4 <- plot_gene_celltype (cpm, "LY6E", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
t5 <- plot_gene_celltype (cpm, "PIGR", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)
t6 <- plot_gene_celltype (cpm, "FCGBP", bulk_mmr, paste0(pseudobulk_dir,"Plots/"), axes = F)

# Make a combined plot of these key genes
combined_key_genes <- plot_grid(t1,t2,t3,t4,t5,t6,ncol = 2, 
                      labels = c("a", "b", "c", "d", "e", "f"), 
                      label_size = 8)

x.grob <- textGrob(expression('Mean log'[2]*' CPM'), 
                   gp=gpar(col="black", fontsize=7))

y.grob <- textGrob("Cell type", 
                   gp=gpar( col="black", fontsize=7),rot=90)

combined_key_genes <- grid.arrange(arrangeGrob(combined_key_genes, left = y.grob, bottom = x.grob))

ggsave(plot = combined_key_genes,"./Results_paper/Plots/Figure S17.pdf", 
       width = 170, height = 170, units = "mm")

# Compare the MMR gene expression to the CosMx pseudobulk
cpm[1:5,1:5]
mmr_median_cpm <- rowMeans(cpm)%>%
  data.frame()%>%
  dplyr::rename(MMR_CPM = 1)

mmr_median_cpm$Gene <-  rownames(cpm)

# Read in the V1 cpm
cpm_v1 <- read_csv("./Results_paper/analyses/Pseudobulk_T_vs_N_HQ/Pseudobulk_log2_CPMs.csv")

v1_median_cpm <- cpm_v1%>%
  select(-Gene)%>%
  as.matrix()%>%
  rowMeans()%>%
  data.frame()%>%
  dplyr::rename(V1_CPM = 1)

v1_median_cpm$Gene <- cpm_v1$Gene

# Read in the V1 cpm
cpm_v2 <- read_csv("./Results_paper/analyses/Pseudobulk_T_vs_N_G_HQ/Pseudobulk_log2_CPMs.csv")

v2_median_cpm <- cpm_v2%>%
  select(-Gene)%>%
  as.matrix()%>%
  rowMeans()%>%
  data.frame()%>%
  dplyr::rename(V2_CPM = 1)

v2_median_cpm$Gene <- cpm_v2$Gene

v1_median_cpm_joined <- left_join(v1_median_cpm, mmr_median_cpm)%>%
  left_join(v2_median_cpm)

cor <-cor.test(v1_median_cpm_joined$MMR_CPM, v1_median_cpm_joined$V1_CPM)
text <- paste0("Pearson's r = ", round(cor$estimate,2))

v1_mmr <- ggplot(data = v1_median_cpm_joined, aes(x = MMR_CPM, y = V1_CPM))+
  geom_point(size = 0.5, alpha = 0.5)+
  blank_theme+
  geom_smooth()+
  labs(x = expression('MMR atlas log'[2]*' CPM'), y = expression('CosMx V1 log'[2]*' CPM'))+
  annotate("text", x = 5, y =14, label = text, size = 2.5)

cor <-cor.test(v1_median_cpm_joined$MMR_CPM, v1_median_cpm_joined$V2_CPM)
text <- paste0("Pearson's r = ", round(cor$estimate,2))

v2_mmr <- ggplot(data = v1_median_cpm_joined, aes(x = MMR_CPM, y = V2_CPM))+
  geom_point(size = 0.5, alpha = 0.5)+
  blank_theme+
  geom_smooth()+
  labs(x = expression('MMR atlas log'[2]*' CPM'), y = expression('CosMx V2 log'[2]*' CPM'))+
  annotate("text", x = 5, y =14, label = text, size = 2.5)

cor <-cor.test(v1_median_cpm_joined$V1_CPM, v1_median_cpm_joined$V2_CPM)
text <- paste0("Pearson's r = ", round(cor$estimate,2))

v1_v2 <- ggplot(data = v1_median_cpm_joined, aes(x = V1_CPM, y = V2_CPM))+
  geom_point(size = 0.5, alpha = 0.5)+
  blank_theme+
  geom_smooth()+
  labs(x = expression('CosMx V1 log'[2]*' CPM'), y = expression('CosMx V2 log'[2]*' CPM'))+
  annotate("text", x = 7.5, y =14, label = text, size = 2.5)


# Compare the sample compositions
mmr_sample_comp <-read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Sample_compositions.csv")%>%
  group_by(Sample_type, Manual_toplevel_pred)%>%
  summarise(mmr_mean_prop = mean(Percent_total_cells))

v1_sample_comp <-read_csv("./Results_paper/analyses/Pseudobulk_T_vs_N_HQ/Sample_cell_type_compositions.csv")%>%
  group_by(Sample_type, Manual_toplevel_pred)%>%
  summarise(v1_mean_prop = mean(Percent_total_cells))

v2_sample_comp <-read_csv("./Results_paper/analyses/Pseudobulk_T_vs_N_G_HQ/Sample_cell_type_compositions.csv")%>%
  group_by(Sample_type, Manual_toplevel_pred)%>%
  summarise(v2_mean_prop = mean(Percent_total_cells))

# Join the sample comps
sample_comps_joined <- left_join(mmr_sample_comp, v1_sample_comp)%>%
  left_join(v2_sample_comp)%>%
  ungroup()

cor <-cor.test(sample_comps_joined$mmr_mean_prop, sample_comps_joined$v1_mean_prop)
text <- paste0("Pearson's r = ", round(cor$estimate,2))

prop_mmr_v1 <- ggplot(data = sample_comps_joined, aes(x = log2(mmr_mean_prop), y = log2(v1_mean_prop), label = Manual_toplevel_pred, colour =Sample_type))+
  geom_text_repel(size = 2)+ 
  blank_theme+
  guides(colour = "none")+
  scale_colour_manual(values = sample_type_cols)+
  geom_abline(intercept =0, linetype = 2, linewidth = 0.5)+
  labs(x = "MMR atlas mean cell type proportion", y = "CosMx V1 mean cell type proportion")+
  annotate("text", x = 0, y =5, label = text, size = 2.5)

cor <-cor.test(sample_comps_joined$mmr_mean_prop, sample_comps_joined$v2_mean_prop)
text <- paste0("Pearson's r = ", round(cor$estimate,2))

prop_mmr_v2 <- ggplot(data = sample_comps_joined, aes(x = log2(mmr_mean_prop), y = log2(v2_mean_prop), label = Manual_toplevel_pred, colour =Sample_type))+
  geom_text_repel(size = 2)+ 
  blank_theme+
  guides(colour = "none")+
  scale_colour_manual(values = sample_type_cols)+
  geom_abline(intercept =0, linetype = 2, linewidth = 0.5)+
  labs(x = "MMR atlas mean cell type proportion", y = "CosMx V2 mean cell type proportion")+
  annotate("text", x = 0, y =5, label = text, size = 2.5)

cor <-cor.test(sample_comps_joined$v1_mean_prop, sample_comps_joined$v2_mean_prop)
text <- paste0("Pearson's r = ", round(cor$estimate,2))

prop_v1_v2 <- ggplot(data = sample_comps_joined, aes(x = log2(v1_mean_prop), y = log2(v2_mean_prop), label = Manual_toplevel_pred, colour =Sample_type))+
  geom_text_repel(size = 2)+ 
  blank_theme+
  guides(colour = "none")+
  scale_colour_manual(values = sample_type_cols)+
  geom_abline(intercept =0, linetype = 2, linewidth = 0.5)+
  labs(x = "CosMx V1 mean cell type proportion", y = "CosMx V2 mean cell type proportion")+
  annotate("text", x = 0, y =5, label = text, size = 2.5)

# Combine the plots 
express_prop_qc <- plot_grid(prop_mmr_v1,prop_mmr_v2, prop_v1_v2,
                             v1_mmr, v2_mmr,v1_v2,
                             labels = c("a", "b", "c", "d", "e", "f"), label_size = 8, align = "hv")

ggsave(plot = express_prop_qc,"./Results_paper/Plots/Figure S CosMx vs MMR cor.pdf", 
       width = 170, height = 110, units = "mm")

# Look into the MMR epi cells a bit more
MMR_epi <- mmr_atlas[,mmr_atlas$clMidwayPr == "Epi"]

dim(MMR_epi)

MMR_epi <- FindVariableFeatures(MMR_epi, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MMR_epi), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MMR_epi)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2

# Scale the data
all.genes <- rownames(MMR_epi)
MMR_epi <- ScaleData(MMR_epi, features = all.genes)

MMR_epi <- RunPCA(MMR_epi, features = VariableFeatures(object = MMR_epi))

print(MMR_epi[["MMR_epi"]], dims = 1:5, nfeatures = 5)

# Plot the loadings (the genes most driving PC 1 and 2)
VizDimLoadings(MMR_epi, dims = 1:2, reduction = "pca")

# Load up the gut cell atlas as a reference
fetal_gut_atlas <- qread("/pvol/andrew/reference/Gut_Atlas_Seurat.qs")

# Find the epithelial cell types
pred_table <- table(fetal_gut_atlas$Integrated_05, fetal_gut_atlas$category)%>%
  data.frame()%>%
  filter(Var2 == "Epithelial")%>%
  filter(Freq > 0)

fetal_gut_atlas

# Focus on the epithelial cells
epi <- fetal_gut_atlas[,fetal_gut_atlas$category == "Epithelial"]

epi$type_cell <- paste0(epi$Age_group, "|", epi$Integrated_05)

# Drop the very low freq cell types
type_table <- table(epi$type_cell)%>%
  data.frame()%>%
  filter(Freq > 5)

epi <- epi[,epi$type_cell %in% type_table$Var1]

# Annotate cells with singleR
# Run singleR and compare against the human primary cell atlas
predictions <- SingleR(test=MMR_epi@assays$RNA$counts,
                       ref=epi@assays$RNA@counts, labels=epi$type_cell, 
                       num.threads = 30,aggr.ref = T)

MMR_epi$Gut_Atlas_pred_main <- predictions$labels

# Add on just a cell type prediction
MMR_epi$Gut_Atlas_pred_ct <- gsub(".*\\|", "", MMR_epi$Gut_Atlas_pred_main )

# Make the PCA plot
DimPlot(MMR_epi, reduction = "pca",group.by = "Gut_Atlas_pred_main", label = T) + NoLegend()

ElbowPlot(MMR_epi,ndims = 50)

MMR_epi <- FindNeighbors(MMR_epi, dims = 1:50)
MMR_epi <- FindClusters(MMR_epi)

MMR_epi <- RunUMAP(MMR_epi, dims = 1:50)
MMR_epi$donor <- gsub("_.*", "", MMR_epi$Sample_name)
DimPlot(MMR_epi, reduction = "umap", group.by = "Gut_Atlas_pred_main", label = T) + NoLegend()
p1 <- DimPlot(MMR_epi, reduction = "umap", group.by = "SPECIMEN_TYPE", label = T) + NoLegend()
p2 <- DimPlot(MMR_epi, reduction = "umap", group.by = "cl295v11SubFull", label = T) + NoLegend()
DimPlot(MMR_epi, reduction = "umap", group.by = "donor", label = T) + NoLegend()

preds_tab <- table(MMR_epi$Gut_Atlas_pred_main, MMR_epi$Sample_name)%>%
  data.frame()

p3 <- FeaturePlot(MMR_epi, features = c("FXYD5"))
p4 <- FeaturePlot(MMR_epi, features = c("PIGR"))
p5 <- FeaturePlot(MMR_epi, features = c("OTOP2"))

epi_split <- p1+p2+p3+p4

ggsave(plot = epi_split, filename = "./Results_paper/analyses/MMR_atlas_T_vs_N/Plots/FXYD5_PIGR_epi_split.pdf", width = 170, height = 150, units = "mm")

# Save the annotated object
qsave(MMR_epi, "/pvol/andrew/reference/GSE178341_Epi_clustered_annotated.qs")

# Look at some individual donors
donor <- MMR_epi[,grepl("C124",MMR_epi$Sample_name)]

dim(donor)

donor <- FindVariableFeatures(donor, selection.method = "vst", nfeatures = 2000)
# Scale the data
all.genes <- rownames(donor)
donor <- ScaleData(donor, features = all.genes)
donor <- RunPCA(donor, features = VariableFeatures(object = donor))
donor <- FindNeighbors(donor, dims = 1:50)
donor <- FindClusters(donor)
donor <- RunUMAP(donor, dims = 1:50)

p1 <- DimPlot(donor, reduction = "umap", group.by = "Gut_Atlas_pred_main", label = T) + NoLegend()
p2 <- DimPlot(donor, reduction = "umap", group.by = "SPECIMEN_TYPE", label = T) + NoLegend()
p3 <- DimPlot(donor, reduction = "umap", group.by = "Gut_Atlas_pred_ct", label = T) 
p4 <- DimPlot(donor, reduction = "umap", label = T)  + NoLegend()
p5 <- DimPlot(donor, reduction = "umap", group.by = "cl295v11SubFull", label = T) + NoLegend()
p6 <- DimPlot(donor, reduction = "umap", label = T)  + NoLegend()

p5+p2

FeaturePlot(donor, features = "PIGR")
FeaturePlot(donor, features = "FXYD5")
FeaturePlot(donor, features = "BEST4")
FeaturePlot(donor, features = "OLFM4")
FeaturePlot(donor, features = "CXCL1")

marks <- FindAllMarkers(donor)

donor_tab <- table(donor$Gut_Atlas_pred_main, donor$SPECIMEN_TYPE)%>%
  data.frame()

p1 +p2

FeaturePlot(donor, features = c("BEST4", "FXYD5", "PIGR", "EPCAM"))

# Try with just the tumour
C106_T <- MMR_epi[,grepl("C106_T",MMR_epi$Sample_name)]

dim(C106_T)

C106_T <- FindVariableFeatures(C106_T, selection.method = "vst", nfeatures = 2000)
# Scale the data
all.genes <- rownames(C106_T)
C106_T <- ScaleData(C106_T, features = all.genes)
C106_T <- RunPCA(C106_T, features = VariableFeatures(object = C106_T))
C106_T <- FindNeighbors(C106_T, dims = 1:50)
C106_T <- FindClusters(C106_T, resolution = 1)
C106_T <- RunUMAP(C106_T, dims = 1:50, min.dist	=0.001)

DimPlot(C106_T)

marks <- FindAllMarkers(C106_T)

p1 <- DimPlot(C106_T, reduction = "umap", group.by = "Gut_Atlas_pred_main", label = T) + NoLegend()
p2 <- DimPlot(C106_T, reduction = "umap", group.by = "cl295v11SubFull", label = T) + NoLegend()
p3 <- DimPlot(C106_T, reduction = "umap", group.by = "Gut_Atlas_pred_ct", label = T) + NoLegend()
p4 <- DimPlot(C106_T, reduction = "umap", label = T)  + NoLegend()

FeaturePlot(C106_T, "CLCA4")

p1+p2+p3+p4

# Do some ligand receptor filtering 
# Read all the genes in the CPDB database
genes <- read_csv("/pvol/andrew/reference/cpdb/v4.1.0/gene_input.csv")%>%
  dplyr::rename(SYMBOL = hgnc_symbol)%>%
  dplyr::select(SYMBOL, uniprot)

interaction <- read_csv("/pvol/andrew/reference/cpdb/v4.1.0/interaction_input.csv")

# Fibroblast T/N genes
t_n <- read_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/toptables/T_Fibro-N_Fibro.csv")%>%
  filter(adj.P.Val < 0.05)%>%
  filter(SYMBOL %in% genes$SYMBOL)%>%
  left_join(genes)%>%
  mutate(ligand = uniprot %in% interaction$partner_a)%>%
  mutate(receptor = uniprot %in% interaction$partner_b)%>%
  filter(ligand == T | receptor == T)%>%
  write_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Fibroblast_tumour_LR_genes.csv")

granulo_fibro <-  read_csv("./Results_paper/analyses/MMR_atlas_Neu/toptables/Granulo_high_Fibro-Granulo_low_Fibro.csv")%>%
  filter(adj.P.Val < 0.05)%>%
  filter(SYMBOL %in% genes$SYMBOL)%>%
  left_join(genes)%>%
  mutate(ligand = uniprot %in% interaction$partner_a)%>%
  mutate(receptor = uniprot %in% interaction$partner_b)%>%
  filter(ligand == T | receptor == T)%>%
  write_csv("./Results_paper/analyses/MMR_atlas_T_vs_N/Fibroblast_granulo_high_LR_genes.csv")

plt <- VlnPlot(mmr_atlas, features = "ANXA1", group.by = "clMidwayPr")
DotPlot(mmr_atlas, features = c("ANXA1", "CEACAM8", "CXCL8", "CXCR2", "CXCR1"), group.by = "clMidwayPr")

ggsave(plot = plt,filename =  paste0(pseudobulk_dir, "/Plots/ANXA1_violin.pdf"), width = 10, height = 6)

print(plt)

# Save the R session info for methods section
writeLines(capture.output(sessionInfo()), paste0("./", "Session_info.txt"))
