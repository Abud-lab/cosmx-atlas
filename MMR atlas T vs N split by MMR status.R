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

# Read in the MMR atlas with logtransformed counts
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

# Take a loot at the MMR atlas and see if you can do any better
# with processing
Idents(mmr_atlas) <- mmr_atlas$PatientTypeID
VlnPlot(mmr_atlas, features = c("nCount_RNA"),pt.size = -1)+NoLegend()

mmr_atlas$Barcode <- colnames(mmr_atlas)
mmr_atlas$Sample_name <- mmr_atlas$PatientTypeID

# Add the metdata to the seurat object
md <- mmr_atlas@meta.data%>%
  # Calculate the MAD values for counts features and mito %
  # Do this for each individual sample
  group_by(Sample_name)%>%
  mutate(m = median(nFeature_RNA))%>%
  mutate(s = mad(nFeature_RNA))%>%
  mutate(robzscore_nFeature_RNA = abs((nFeature_RNA - m) / (s)))%>%
  mutate(m = median(nCount_RNA))%>%
  mutate(s = mad(nCount_RNA))%>%
  mutate(robzscore_nCount_RNA = abs((nCount_RNA - m) / (s)))%>%
  mutate(m = median(percent.mt))%>%
  mutate(s = mad(percent.mt))%>%
  mutate(robzscore_percent.mt = abs((percent.mt - m) / (s)))%>%
  ungroup()%>%
  data.frame()%>%
  # Drop orig.ident since it sometimes is present and seems to cause an error
  select(-any_of(c("orig.ident")))

# Reset the rownames
rownames(md) <- md$Barcode
mmr_atlas@meta.data <- md

mmr_atlas@assays$RNA@data[1:20,1:20]

min_QC_robz <- 2

# Subset down based on QC cutoffs for each sample
mmr_atlas <- subset(mmr_atlas, subset = robzscore_nFeature_RNA < min_QC_robz & robzscore_percent.mt < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

plot_gene_celltype <- function(cpm, gene, bulk_anno, outdir){
  
  # Make a plot of a single gene
  single_gene <- cpm[gene,]%>%
    data.frame(check.rows = F)%>%
    rownames_to_column("Sample_cell")%>%
    dplyr::rename(CPM = 2)%>%
    left_join(bulk_anno)%>%
    group_by(Manual_toplevel_pred, Sample_type)%>%
    mutate(median_cpm = median(CPM))%>%
    ungroup()%>%
    arrange(-median_cpm)%>%
    mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))
  
  plt <- ggplot(data = single_gene, aes(x = Manual_toplevel_pred, y = CPM))+ 
    geom_jitter(width = 0.1, height = 0, aes(colour = Sample_name))+
    facet_wrap(~Sample_type)+
    guides(colour ="none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Cell type", y = "Mean log2 CPM", fill = "Sample", title = gene)+
    stat_summary(fun = "median", fun.min = "median", fun.max= "median", linewidth= 0.3, geom = "crossbar")
  
  ggsave(paste0(outdir, gene, " cell type expression.pdf"), width = 9, height = 5)
  
  print(plt)
  
}

# Figure out differential niche usage between conditions
# Make a pseudobulk object ----
# Look at left vs right DE analysis
pseudobulk_dir <- "/oldvol/apattison/Data/Spatial/Results_paper/analyses/MMR_atlas_T_split_vs_N/"

blank_theme <- theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

table(mmr_atlas$clMidwayPr)
table(mmr_atlas$Sex)
table(mmr_atlas$MMRStatus)
table(mmr_atlas$TumorStage)
table(mmr_atlas$Age)
table(mmr_atlas$SPECIMEN_TYPE)
mmr_atlas$Manual_toplevel_pred <- mmr_atlas$clMidwayPr

# Study annotations
unique(mmr_atlas$clMidwayPr)

to_bulk_anno <- mmr_atlas@meta.data%>%
  dplyr::select(Barcode = sampleID, Manual_toplevel_pred, Sample_name = PatientTypeID, 
                Sample_type = SPECIMEN_TYPE, MMRStatus, Sex, TumorStage)%>%
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
  labs(x = "Sample", y = "Percent of total", fill = "Cell type")+
  blank_theme+
  scale_fill_manual(values = col_vector)+
  guides(fill = guide_legend(ncol = 2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plt

save <- paste0(pseudobulk_dir, "/Plots/Per sample cell type composition.pdf")
ggsave(filename = save, plot = plt, width = 8, height = 5)

# Keep only the T samples for now
#mmr_atlas <- mmr_atlas[,mmr_atlas$PatientTypeID %in% to_bulk_anno$Sample_name]
mmr_atlas$MMRStatus <- replace(mmr_atlas$MMRStatus, mmr_atlas$SPECIMEN_TYPE == "N", "N")

mmr_atlas$type <-mmr_atlas$MMRStatus
mmr_atlas$type <- factor(mmr_atlas$type, levels = c("N", "MMRp","MMRd"))
table(mmr_atlas$type)

# Run sccomp with contrasts for tumour vs normal
sc_result <- mmr_atlas |>
  sccomp_glm( 
    formula_composition = ~type, 
    .sample = PatientTypeID,
    .cell_group = Manual_toplevel_pred, 
    bimodal_mean_variability_association = F,
    cores = 30
  )

# Save the results
#saveRDS(sc_result, paste0(pseudobulk_dir, "sccomp results.rds"))

plots <- plot_summary(sc_result) 

# Plot the DE results
CT_DA <- plots$boxplot[[1]]
CT_DA <- CT_DA+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(pseudobulk_dir, "Plots/sccomp DA boxplot.pdf"), 
    width = 12, height = 7)
CT_DA
dev.off()

# Get the raw counts
pb_counts <- mmr_atlas@assays$RNA@counts

# Remove the single cell object save RAM
rm(mmr_atlas)

# Function to get pseudobulk counts per sample_cell combo
get_pb <- function(i, to_bulk_anno){
  
  d_t <- unique(to_bulk_anno$Sample_cell)[i]
  
  # Filter the annotation to keep only cells from a certain sample/cluster
  anno_filt <- to_bulk_anno%>%
    filter(Sample_cell == d_t)
  
  # Keep only the counts that match the sample/cluster
  counts_dt <- pb_counts[,anno_filt$Barcode]
  
  # Skip over single cell groupings
  if(is.null(dim(counts_dt))){
    summed <- data.frame(counts_dt)%>%
      dplyr::rename(!!d_t := 1)
  }
  
  else{
    
    summed <- rowSums(counts_dt)%>%
      data.frame()%>%
      dplyr::rename(!!d_t := 1)
    
  }
  
  if(i == 1){
    summed <- summed%>%
      rownames_to_column("Gene")
    
  }
  
  return(summed)
  
}

# Get the number of samples to test
nums <- 1:length(unique(to_bulk_anno$Sample_cell))

# Parallelise making pseudobulk counts
pb_list <- mclapply(X = nums, FUN = get_pb, to_bulk_anno, mc.cores = 15)

# Bind the PB counts
bound <- bind_cols(pb_list)

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
  select(Sample_type, Sample_name,Donor)%>%
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

# Keep the MMR samples
bulk_anno <- bulk_anno%>%
  mutate(MMRStatus = replace(MMRStatus, Sample_type == "N", "N"))%>%
  mutate(Sample_type = MMRStatus)%>%
  mutate(Sample_type_cell = paste0(Sample_type, "_", Manual_toplevel_pred))

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
design <- model.matrix(~0 + Sample_type_cell + Sex, data = bulk_anno)

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
fit <- voomLmFit(rnaseq, design = design, plot = T)

# Voom for gsea
v <- voom(rnaseq, design, plot=TRUE)
#fit <- lmFit(v, design)

# Automate a contrast matrix
# Make sure baseline is the first group
unique(bulk_anno$Sample_type)

contrasts_manual <- c("MMRd-N","MMRp-N","MMRd-MMRp")

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
  cont_first <- gsub("MMRd_|MMRp", "", cont_first)
  
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
# system(paste0("rm -r /oldvol/apattison/Data/Single_cell/MMR_atlas/Results/T_vs_N/glimma/MDS_files"))

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

# Define a function to shut up some other functions
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

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

# Function to run GSEA for a contrast
run_GSEA <- function(contrast , collection, rnaseq, v, design, cont.matrix){
  
  collection_name <- gsub(".Hs.symbols.gmt|.gmt","", basename(collection))
  
  gene_set <- quiet(GSA.read.gmt(collection))
  gene_set_formatted <- gene_set$genesets
  names(gene_set_formatted) <- gene_set$geneset.names
  indexed <- ids2indices(gene.sets = gene_set_formatted, identifiers = rownames(rnaseq$counts), remove.empty=TRUE)
  
  camera_result <- camera(y = v ,index = indexed, design = design, contrast = cont.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    filter(FDR <= 0.05)%>%
    mutate(Contrast= contrast)
  
  write_csv(camera_result, paste0(pseudobulk_dir,"gsea/camera/",collection_name, "_", contrast,".csv"))
  
  
}


# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)
  
  lapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix)
  
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
  mutate(cell_type = gsub(".*-N_|.*-MMRp_", "", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

# Look at the MMR gene differences
d_vs_p <- camera_compiled%>%
  filter(grepl("MMRp", contrast) & grepl("MMRd", contrast))

no_granulo <- d_vs_p%>%
  filter(cell_type!="Granulo")

# Plot the expression of a given gene in a given celltype
cpm <- edgeR::cpm(rnaseq, log =T)

# Pull the HIF1A TFT gene sets and see what's changed
plot_gene_celltype (cpm, "IL23A", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "TNF", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCR2", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL8", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL2", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL3", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL5", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "CXCL13", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "VCAM1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "SPP1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "VEGFA", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "SPP1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "NRG1", bulk_anno, paste0(pseudobulk_dir,"Plots/"))
plot_gene_celltype (cpm, "HIF1A", bulk_anno, paste0(pseudobulk_dir,"Plots/"))



