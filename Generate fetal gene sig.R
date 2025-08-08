library(Seurat)
library(tidyverse)
library(qs)
# library(SeuratData)
library(SeuratDisk)
library(SingleR)
library(celldex)
#library(spatstat)
#library(viridis)
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

# Try and derive a fetal gene signature
outdir <- "/homevol/apattison/Data/Spatial/Results_paper/analyses/fetal_sig_atlas/Results"

# Make the direcotry
system(paste0("mkdir -p ", outdir))

fetal_gut_atlas <- qread("/pvol/andrew/reference/Gut_Atlas_Seurat.qs")

md <- fetal_gut_atlas@meta.data

# Looks like I do indeed have the raw counts here
fetal_gut_atlas@assays$RNA@counts[1:20,1:50]

# Grab just the 
unique(md$Fraction)
unique(md$category)
unique(md$Age)
# Focus on the epithelial cells
epi <- fetal_gut_atlas[,fetal_gut_atlas$category == "Epithelial"]

dim(epi)
rownames(epi)[1:5]

# Get robust Z scores for the metadata
epi[["percent.mt"]] <- PercentageFeatureSet(epi, pattern = "^MT-")
# Add in percent ribo
epi[["percent.ribo"]] <- PercentageFeatureSet(epi, pattern = "^RPL|^RPS")

# Figure out QC per-sample
md <- epi@meta.data%>%
  # Calculate the MAD values for counts features and mito %
  # Do this for each individual sample
  group_by(sample.name)%>%
  mutate(m = median(nFeature_RNA))%>%
  mutate(s = mad(nFeature_RNA))%>%
  mutate(robzscore_nFeature_RNA = abs((nFeature_RNA - m) / (s)))%>%
  mutate(m = median(nCount_RNA))%>%
  mutate(s = mad(nCount_RNA))%>%
  mutate(robzscore_nCount_RNA = abs((nCount_RNA - m) / (s)))%>%
  mutate(m = median(percent.mt))%>%
  mutate(s = mad(percent.mt))%>%
  mutate(robzscore_percent.mt = abs((percent.mt - m) / (s)))%>%
  mutate(m = median(percent.ribo))%>%
  mutate(s = mad(percent.ribo))%>%
  mutate(robzscore_percent.ribo = abs((percent.ribo - m) / (s)))%>%
  ungroup()%>%
  data.frame()

rownames(md) <- md$Barcode
epi@meta.data <- md

# Filter down the data
total_cells <- nrow(md)
min_features <- 50
sum(md$nFeature_RNA < min_features)
min_QC_robz <- 2.5

VlnPlot(epi, features = c("percent.mt"),
        pt.size = -1, group.by = "sample.name")

VlnPlot(epi, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "sample.name")+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

epi <- subset(epi, subset = !is.na(robzscore_nFeature_RNA) & nFeature_RNA > min_features & robzscore_nFeature_RNA < min_QC_robz & robzscore_percent.mt < min_QC_robz & robzscore_nCount_RNA < min_QC_robz)

VlnPlot(epi, features = c("robzscore_percent.mt", "robzscore_nFeature_RNA", "robzscore_nCount_RNA"),
        pt.size = -1, group.by = "sample.name")+
  geom_hline(yintercept = min_QC_robz, linetype = 2)

md <- md %>%
  filter(!is.na(robzscore_percent.mt))

# Get a summary of the filtering
total_cells_filtered <- ncol(epi)
filter_summary <- data.frame("Minimum features" = min_features,
                             "Robust Z score cutoff" = min_QC_robz,
                             "Cells less than min features" = sum(md$nFeature_RNA < min_features),
                             "Cells less than min features robust Z score" = sum(md$robzscore_nFeature_RNA > min_QC_robz),
                             "Cells less than min percent mt robust Z score" = sum(md$robzscore_percent.mt > min_QC_robz),
                             "Cells less than min count robust Z score" = sum(md$robzscore_nCount_RNA > min_QC_robz),
                             "Total cells" = total_cells,
                             "Total cells after filtering" = total_cells_filtered
)%>%
  gather("Description", "Count")%>%
  mutate(Description = gsub("\\.", " ", Description))%>%
  write_csv(paste0(outdir, "/cell filtering summary.csv"))

filter_summary

pseudobulk_dir <- "/homevol/apattison/Data/Single_cell/fetal_sig_atlas/Results/Fetal_vs_adult_pseudobulk/"

blank_theme <- theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

# Set the celltype level
epi$Manual_toplevel_pred <- epi$Integrated_05

# Study annotations
unique(epi$Manual_toplevel_pred)

table(epi$Region, epi$Age_group)

to_bulk_anno <- epi@meta.data%>%
  dplyr::select(Barcode, Manual_toplevel_pred, Sample_name = sample.name, 
                Sample_type = Age_group, Region, Sex = Gender)%>%
  mutate(Sample_cell = paste0(Sample_name, "_", Manual_toplevel_pred))

# Keep only conditions where we have a decent number of cells
to_drop <- to_bulk_anno %>%
  group_by(Sample_cell,Sample_type)%>%
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
plt <- ggplot(data = cell_type_counts, aes(x = Sample_name, y = Percent_total_cells, fill = Manual_toplevel_pred))+
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

# Conver to character to drop empty factor levels
epi$type <-as.character(epi$Age_group)
epi$Manual_toplevel_pred <- as.character(epi$Manual_toplevel_pred)
epi$sample.name <- as.character(epi$sample.name)
table(epi$sample.name)

# Run sccomp with contrasts for age groups
sc_result <- epi |>
  sccomp_glm( 
    formula_composition = ~type, 
    .sample = sample.name,
    .cell_group = Manual_toplevel_pred, 
    bimodal_mean_variability_association = F,
    cores = 5
  )

plots <- plot_summary(sc_result) 

# Plot the DE results
plots$boxplot

pdf(paste0(pseudobulk_dir, "Plots/sccomp DA boxplot.pdf"), 
    width = 7, height = 4.5)
plots$boxplot
dev.off()

# Get the raw counts
pb_counts <- epi@assays$RNA@counts

# Remove the single cell object save RAM
rm(epi)
gc()

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

bulk_anno <- data.frame(Sample_cell = colnames(bound_counts), check.names = F)%>%
  left_join(condensed_SC_anno)%>%
  mutate(Sample_type = gsub(" ", "_", Sample_type))%>%
  mutate(Manual_toplevel_pred = gsub(" |/", "_", Manual_toplevel_pred))%>%
  mutate(Sample_name = as.character(Sample_name))%>%
  filter(Sample_name != "nan")%>%
  mutate(Sample_type_cell = paste0(Sample_type, "_", Manual_toplevel_pred))%>%
  # Save the Pseudobulk annotation
  write_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# bulk_anno <- read_csv(paste0(pseudobulk_dir, "Pseudobulk_annotation.csv"))

# Make a DGElist
rnaseq <- DGEList(bound_counts[,bulk_anno$Sample_cell])

# Remove weird characters
rownames(rnaseq) <- gsub("\\+", "_pos_", rownames(rnaseq))
colnames(rnaseq) <- gsub("\\+", "_pos_", colnames(rnaseq))
colnames(rnaseq) <- gsub("\\'|-|\\(|\\)| ", "_", colnames(rnaseq))
bulk_anno$Sample_cell <- gsub("\\+", "_pos_", bulk_anno$Sample_cell)
bulk_anno$Sample_cell <- gsub("\\'|-|\\(|\\)| ", "_", bulk_anno$Sample_cell)

# Sanity check
sum(colnames(rnaseq) == bulk_anno$Sample_cell) == length(bulk_anno$Sample_cell)

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
system(paste0("mkdir -p ", pseudobulk_dir,"glimma/"))
mds_save <- paste0(paste0(pseudobulk_dir,"glimma/", "MDS.html"))
htmlwidgets::saveWidget(glimmaMDS(rnaseq, groups = bulk_anno, labels = bulk_anno$Sample_cell), mds_save)

# Normalise and fit linear model
v <- voom(rnaseq, design, plot=TRUE)
fit <- lmFit(v, design)

# Automate a contrast matrix
# Make sure baseline is the first group
unique(bulk_anno$Sample_type)

# Try and derive a fet sig here. 
# Might be a good idea to also compare within fetal cell types
contrasts_manual <- "Second_trim-Adult"

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

unique(bulk_anno$Sample_type_cell)

# Make some new contrast lines that compares everything to everything else
contrast_lines2 <- character()
all_cells_all_groups <- unique(bulk_anno$Sample_type_cell)
for(i in 1:length(all_cells_all_groups)){
  
  cell_type_group <- all_cells_all_groups[i]
  
  rest <- all_cells_all_groups[all_cells_all_groups != cell_type]
  
  restlen <- length(rest)
  
  rest <- paste(rest, collapse = "+")
  
  contrast_line <- paste0(cell_type_group, " = ",cell_type_group, "-(", rest, ")/", restlen)
  
  contrast_lines2[i] <- c(contrast_line)
}

# Add in some within age group contrasts
contrast_lines_use <- c(contrast_lines, contrast_lines2,
                    "Adult_Stem_cells-Adult_Colonocyte",
                    "First_trim_Stem_cells-First_trim_Colonocyte",
                    "Second_trim_Stem_cells-Second_trim_Colonocyte",
                    "Adult_Stem_cells-Adult_Enterocyte",
                    "First_trim_Stem_cells-First_trim_Enterocyte",
                    "Second_trim_Stem_cells-Second_trim_Enterocyte")  

cont.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                              (contrast_lines_use),levels=list(design))))

colnames(cont.matrix) <- gsub(" = .*", "", colnames(cont.matrix))
colnames(cont.matrix) <- gsub("__$", "", colnames(cont.matrix))
colnames(cont.matrix) <- gsub("__|___", "_", colnames(cont.matrix))

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

treat_fit <- treat(fit.cont, lfc = 1)

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

system(paste0("mkdir -p ", pseudobulk_dir, "/toptables/"))
system(paste0("mkdir -p ", pseudobulk_dir, "/treat_lfc_1/"))

system(paste0("rm ", pseudobulk_dir, "/toptables/*.csv"))

VOL <- paste0(pseudobulk_dir, "glimma/volcano/")
system(paste0("mkdir -p ", VOL))

# Make a vector of all contrast types to remove
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
  
  output_treat <- paste0(pseudobulk_dir, "treat_lfc_1/", contrast, ".csv")
  
  toptreat <- topTreat(treat_fit,coef=contrast,sort.by="p",number = Inf)%>%
    rownames_to_column("SYMBOL")%>%
    write_csv(output_treat)
  
  conts <- gsub(all_cty, "", contrast)
  
  # Get back to just the cell type
  cont_first <- gsub("-.*", "", conts)
  cont_second <- gsub(".*-", "", conts)
  
  conts_both <- paste0(cont_first, "|", cont_second)
  
  rnaseq_filt <- rnaseq[,grepl(conts_both, colnames(rnaseq))]
  
  bulk_filt <- bulk_anno%>%
    filter(Sample_cell %in% colnames(rnaseq_filt))
  
  MA_fac <- factor(bulk_filt$Sample_type_cell, levels = unique(bulk_filt$Sample_type_cell))
  
  vol_save <- paste0(pseudobulk_dir, "glimma/volcano/",contrast,"_", "glimma_volcano.html")
  htmlwidgets::saveWidget(glimmaVolcano(fit.cont, coef = contrast,main = gsub("_"," ",contrast),
                                        counts = round(rnaseq_filt$counts),
                                        dge = rnaseq_filt, groups = MA_fac), vol_save)
  
}

# Drop the glimma files
system(paste0("rm -r /homevol/apattison/Data/Single_cell/fetal_sig_atlas/Results/Fetal_vs_adult_pseudobulk/glimma/*/*_files"))
system(paste0("rm -r /homevol/apattison/Data/Single_cell/fetal_sig_atlas/Results/Fetal_vs_adult_pseudobulk/glimma/MDS_files"))

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
  mutate(cell_type = gsub(all_cty, "", contrast))

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

# Compile the top treat tables
all_toptreat <- list.files(paste0(pseudobulk_dir, "treat_lfc_1/"), full.names = T)

tt_list <- list()
for(i in 1:length(all_toptreat)){
  
  contrast <- gsub(".csv", "", basename(all_toptreat[i]))
  
  tt <- read_csv(all_toptreat[i])%>%
    mutate(contrast = contrast)
  
  tt_list[[i]] <- tt
  
  
}

# Compile toptables and save the significant results
toptreat_compiled <- bind_rows(tt_list)%>%
  mutate(cell_type = gsub(all_cty, "", contrast))%>%
  mutate(cell_type = gsub(".*-", "", cell_type))

toptreat_signif <- toptreat_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptreat_lfc_1_significant_genes.csv"))

CLU <- toptables_compiled%>%
  filter(SYMBOL == "CLU")

# Define a function to shut up some other functions
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Gene set collections
gsea_gmt_dir <- "/olvol/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
keep <- c(5,1,20, 28,31,9,16)
#keep <- c(31)
collections <- collections[keep]

# Make a directory
system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))

# Function to run GSEA for a contrast
run_GSEA <- function(contrast , collection, rnaseq, v, design, cont.matrix){
  
  collection_name <- gsub(".Hs.symbols.gmt","", basename(collection))
  
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

# Remove the full sc object to save space
rm(fetal_gut_atlas)
gc()

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

remove_names <- all_cty

# Compile toptables and save the significant results
camera_compiled <- bind_rows(clist)%>%
  mutate(FDR = as.numeric(FDR))%>%
  arrange(FDR)%>%
  # Fix the cell type naming
  mutate(cell_type = gsub(".*-", "", contrast))%>%
  #mutate(cell_type = gsub(to_remove, "", cell_type))%>%
  mutate(cell_type = gsub(remove_names, "", cell_type))%>%
  mutate(cell_type = gsub("^_", "", cell_type))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

unique(toptables_signif$contrast)

# Make some gene sets to use with sigle cell/spatial data
Adult_stem_vs_Adult_Enterocyte_up <- toptreat_signif %>%
  filter(logFC > 0)%>%
  filter(contrast == "Adult_Stem_cells-Adult_Enterocyte")

Adult_stem_vs_Adult_Enterocyte_up <- list(id = "Adult_stem_vs_Adult_Enterocyte_up", name = "Adult_stem_vs_Adult_Enterocyte_up", genes = Adult_stem_vs_Adult_Enterocyte_up$SYMBOL)

Adult_stem_vs_Adult_Enterocyte_down <- toptreat_signif %>%
  filter(logFC < 0)%>%
  filter(contrast == "Adult_Stem_cells-Adult_Enterocyte")

Adult_stem_vs_Adult_Enterocyte_down <- list(id = "Adult_stem_vs_Adult_Enterocyte_down", name = "Adult_stem_vs_Adult_Enterocyte_down", genes = Adult_stem_vs_Adult_Enterocyte_down$SYMBOL)

# Change the P value to get a smaller set of genes
Second_trim_Enterocyte_vs_Adult_Enterocyte_up <- toptreat_signif %>%
  filter(logFC > 0)%>%
  filter(contrast == "Second_trim_Enterocyte-Adult_Enterocyte")

Second_trim_Enterocyte_vs_Adult_Enterocyte_up <- list(id = "Second_trim_Enterocyte_vs_Adult_Enterocyte_up", name = "Second_trim_Enterocyte_vs_Adult_Enterocyte_up", genes = Second_trim_Enterocyte_vs_Adult_Enterocyte_up$SYMBOL)

Second_trim_Enterocyte_vs_Adult_Enterocyte_down <- toptreat_signif %>%
  filter(logFC < 0)%>%
  filter(contrast == "Second_trim_Enterocyte-Adult_Enterocyte")

Second_trim_Enterocyte_vs_Adult_Enterocyte_down <- list(id = "Second_trim_Enterocyte_vs_Adult_Enterocyte_down", name = "Second_trim_Enterocyte_vs_Adult_Enterocyte_down", genes = Second_trim_Enterocyte_vs_Adult_Enterocyte_down$SYMBOL)

# Make the gene sets into a GMT file
gmt <- list("Adult_stem_vs_Adult_Enterocyte_up" = Adult_stem_vs_Adult_Enterocyte_up, 
            "Adult_stem_vs_Adult_Enterocyte_down" = Adult_stem_vs_Adult_Enterocyte_down,
            "Second_trim_Enterocyte_vs_Adult_Enterocyte_up" = Second_trim_Enterocyte_vs_Adult_Enterocyte_up,
            "Second_trim_Enterocyte_vs_Adult_Enterocyte_down" = Second_trim_Enterocyte_vs_Adult_Enterocyte_down)

class(gmt) <- "GMT"

is.GMT(gmt)

write.GMT(gmt, "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt")

# Make a gene set list for all the different cell types in different conditions
toptreat_signif_celltype <- toptreat_signif%>%
  filter(!grepl("-", contrast))%>%
  filter(adj.P.Val < 1E-5)
  
contrasts <- unique(toptreat_signif_celltype$contrast)
gmt_list <- list()
for(i in 1:length(contrasts)){
  
  contrast_use <- contrasts[i]
  
  cont_genes <- toptreat_signif_celltype%>%
    filter(contrast == contrast_use)
  
  cont_genes_up <- cont_genes %>%
    filter(logFC > 0)
  
  cont_genes_down <- cont_genes %>%
    filter(logFC < 0)
  
  cont_genes_up_list <- list(id =paste0(contrast_use, "_up"), name = paste0(contrast_use, "_up"), genes = cont_genes_up$SYMBOL)
  cont_genes_down_list <- list(id =paste0(contrast_use, "_down"), name = paste0(contrast_use, "_down"), genes = cont_genes_down$SYMBOL)
  
  if(nrow(cont_genes_up) > 0 ){
    gmt_list[[i]] <- cont_genes_up_list
    
  }
  
  i2 <- i+length(contrasts)
  
  if(nrow(cont_genes_down) > 0 ){
    gmt_list[[i2]] <- cont_genes_down_list
  }
}
length(gmt_list)

gmt_dropped <- vctrs::list_drop_empty(gmt_list) 

gmt_names <- character()

for(i in 1:length(gmt_dropped)){
  
  gmt_names[i] <- gmt_dropped[[i]]$name
  
}

names(gmt_dropped) <- gmt_names

# Make into a GMT
class(gmt_dropped) <- "GMT"

is.GMT(gmt_dropped)

write.GMT(gmt_dropped, "/pvol/andrew/reference/msigdb/msigdb_v2023.1.Hs_GMTs/Gut_cell_atlas_epithelial_cell_types.gmt")

# Compare the fetal and stem gene sets
fetal <- read_csv("/olvol/Data/Single_cell/fetal_sig_atlas/Results/Fetal_vs_adult_pseudobulk/toptables/Second_trim_Colonocyte-Adult_Colonocyte.csv")%>%
  select(SYMBOL, fetal_fc = logFC, fetal_FDR = adj.P.Val)

stem <-  read_csv("/olvol/Data/Single_cell/fetal_sig_atlas/Results/Fetal_vs_adult_pseudobulk/toptables/Adult_Stem_cells-Adult_Colonocyte.csv")%>%
  select(SYMBOL, stem_fc = logFC, stem_FDR = adj.P.Val)

sum(stem$SYMBOL %in% fetal$SYMBOL)

both <- full_join(fetal, stem)%>%
  filter(stem_FDR < 0.05 | fetal_FDR < 0.05)%>%
  mutate(Direction_stem = ifelse(stem_fc >0, "Up", "Down"))%>%
  mutate(Direction_fetal = ifelse(fetal_fc  >0, "Up", "Down"))%>%
  mutate(Direction_both = paste0(Direction_stem, Direction_fetal))%>%
  write_csv("/olvol/Data/Single_cell/fetal_sig_atlas/Results/Fetal_vs_adult_pseudobulk/Fetal_vs_stem_DE_comparison.csv")

table(both$Direction_both)

cor.test(both$fetal_fc, both$stem_fc)

ggplot(data = both, aes(x = fetal_fc, y = stem_fc, colour = Direction_both))+
  geom_point()

