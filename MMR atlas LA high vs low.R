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

# Source the R functions I have created
source("~/Code/Spatial/Paper_code_V3/Functions.R")

# Set a homedir for the analysis so as to no use hardcoded file paths
setwd("/oldvol/apattison/Data/Spatial")

# Read in the MMR atlas with logtransformed counts
mmr_atlas <- qread("/pvol/andrew/reference/GSE178341_lognorm_annotated.qs")

pseudobulk_dir <- "./Results_paper/analyses/MMR_atlas_TLS/"

system(paste0("mkdir -p ", pseudobulk_dir, "/Plots"))

t_tab <- mmr_atlas@meta.data%>%
  select(Sample_name = PatientTypeID, cl295v11SubFull, MMRStatus, Age)%>%
  mutate(T_N = gsub(".*_", "", Sample_name))%>%
  mutate(MMRStatus = replace(MMRStatus, is.na(MMRStatus), "N"))%>%
  mutate(iscTNI02 = ifelse(cl295v11SubFull == "cTNI02 (CD4+ IL7R+SELL+)", 1, 0))%>%
  group_by(Sample_name, MMRStatus)%>%
  mutate(Total = n())%>%
  group_by(Sample_name, Total, MMRStatus, Age)%>%
  summarise(Total_cTNI02 = sum(iscTNI02))%>%
  ungroup()%>%
  mutate(Pct_cTNI02 = Total_cTNI02/Total*100)%>%
  mutate(MMRStatus = factor(MMRStatus, levels = c("N", "MMRp", "MMRd")))%>%
  group_by(MMRStatus)%>%
  mutate(Mean_cTNI02_pct = mean(Pct_cTNI02))%>%
  ungroup()

# The % of DCs that are mRegs is not that different across samples
ggplot(data = t_tab, aes(x = MMRStatus, y = Pct_cTNI02))+
  geom_boxplot()+
  geom_jitter(height = 0)

ggplot(data = t_tab, aes(x = Age, y = log2(Pct_cTNI02)))+
  geom_point()

# Find the cell type most correlated with cTNI02+ T cells
hist(t_tab$Pct_cTNI02, breaks = 100)

median_cTNI02 <- quantile(t_tab$Pct_cTNI02, .5)

# Get the TLS high tumours
TLS_high <- t_tab$Sample_name[t_tab$Pct_cTNI02 > 1]
  
to_bulk_anno <- mmr_atlas@meta.data%>%
  dplyr::select(Barcode = sampleID, Manual_toplevel_pred = clMidwayPr, Sample_name = PatientTypeID,
                Sample_type = SPECIMEN_TYPE, MMRStatus, Sex)%>%
  mutate(T_N = gsub(".*_", "", Sample_name))%>%
  mutate(T_N = gsub("A$|B$", "", T_N))%>%
  mutate(Sample_cell = paste0(Sample_name, "_", Manual_toplevel_pred))%>%
  mutate(Sample_type = ifelse(Sample_name %in% TLS_high, "TLS_high", "TLS_low"))%>%
  mutate(Sample_type = paste0(Sample_type, "_", T_N))

# Keep only conditions where we have a decent number of cells
to_drop <- to_bulk_anno %>%
  group_by(Sample_cell,Sample_type, T_N)%>%
  summarise(count = n())%>%
  ungroup()

keep <- to_drop%>%
  filter(count >= 10)

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

cell_type_counts_unique <- cell_type_counts %>%
  filter(!duplicated(cell_type_counts$Sample_name))

table(cell_type_counts_unique$Sample_type)

# Get a big vector of different colours
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Cell type composition plot
plt <- ggplot(data = cell_type_counts, aes(x = Sample_name, y = Percent_total_cells, fill =Manual_toplevel_pred))+
  geom_bar(stat = "identity")+
  facet_wrap(~Sample_type, scales = "free_x")+
  labs(x = "Sample", y = "Percent of total", fill = "Niche")+
  blank_theme+
  scale_fill_manual(values = cell_type_colors)+
  guides(fill = guide_legend(ncol = 2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plt

save <- paste0(pseudobulk_dir, "/Plots/Per sample cell type composition.pdf")
ggsave(filename = save, plot = plt, width = 8, height = 5)

sample_ct <- cell_type_counts%>%
  dplyr::select(Sample_name, Sample_type)%>%
  filter(!duplicated(Sample_name))

# Look at TLS in T vs N with sccomp
cell_type_table <- mmr_atlas@meta.data%>%
  left_join(sample_ct)

mmr_atlas$TLS_type <- cell_type_table$Sample_type

unique(mmr_atlas$TLS_type)

mmr_atlas$TLS_type <- factor(mmr_atlas$TLS_type, 
                             levels = c("TLS_low_N", "TLS_low_T", "TLS_high_N", "TLS_high_T"))

# Shorten full cell names
mmr_atlas$cl295v11SubFull_short <- gsub(".*\\(|)", "",mmr_atlas$cl295v11SubFull)
unique(mmr_atlas$cl295v11SubFull_short )

mmr_atlas$MMRStatus <- replace(mmr_atlas$MMRStatus, is.na(mmr_atlas$MMRStatus) ,"N")

# Run sccomp with contrasts for tumour vs normal
sc_result <- mmr_atlas |>
  sccomp_glm(
    formula_composition = ~TLS_type,
    .sample = PatientTypeID,
    .cell_group = cl295v11SubFull_short,
    bimodal_mean_variability_association = T,
    cores = 30
  )

plots <- sccomp_test(sc_result) 

CT_DA <- plots |> 
  sccomp_boxplot(factor = "TLS_type")

CT_DA <- CT_DA+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(paste0(pseudobulk_dir, "Plots/sccomp DA boxplot.pdf"),
    width = 30, height = 10)
CT_DA
dev.off()

sc_result_signif <- plots%>%
  filter(c_FDR <0.025)%>%
  filter(parameter == "TLS_typeTLS_high_T")

plots2 <- sccomp_boxplot(sc_result[sc_result$cl295v11SubFull_short %in% sc_result_signif$cl295v11SubFull_short,],factor = "TLS_type")

plots2 <- plots2+
  blank_theme+
  labs(title = NULL, shape = "Outlier", fill ='Signif')+
  theme(legend.key.size = unit(2, 'mm'),
        axis.text.y =element_text(size=4),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "none")+
  facet_wrap(vars(cl295v11SubFull_short), nrow = 4, scales = "free_y")

pdf(paste0(pseudobulk_dir, "Plots/sccomp DA boxplot signif.pdf"),
    width = 5, height = 6)
plots2
dev.off()


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
design <- model.matrix(~0 + Sample_type_cell + Sex + MMRStatus, data = bulk_anno)

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

contrasts_manual <- "TLS_high-TLS_low"

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
  scale_fill_manual(values = c("blue", "red"))

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

  #Get back to just the cell type
  cont_first <- gsub("-.*", "", contrast)
  cont_first <- gsub("TLS_high_", "", cont_first)

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
  mutate(cell_type = gsub(".*TLS_low_", "", contrast))

toptables_signif <- toptables_compiled %>%
  filter(adj.P.Val < 0.05)%>%
  arrange(adj.P.Val)%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_toptables_significant_genes.csv"))

# Gene set collections
gsea_gmt_dir <- "~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/"
collections <- list.files(gsea_gmt_dir, full.names = T,pattern = "*.symbols.gmt")
# Keep only some collections
keep <- c(5,1,20, 28,31,9,16)
#keep <- c(31)
collections <- collections[keep]
# Add on my fetal signature
collections <- c(collections,"~/Data/Reference/msigdb/msigdb_v2023.1.Hs_GMTs/Fetal_gene_sigs.gmt")

# Loop over gene sets and run GSEA for each contrast
for(collection in collections){
  print(collection)

  mclapply(X = colnames(cont.matrix),run_GSEA, collection, rnaseq, v, design, cont.matrix,pseudobulk_dir, mc.cores = 10)

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
  mutate(cell_type = gsub(".*TLS_low_", "", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_camera.csv"))

GO <- camera_compiled%>%
  filter(grepl("^GO", `Gene set`))

HALLMARK <- camera_compiled%>%
  filter(grepl("^HALLMARK", `Gene set`))

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
  mutate(cell_type = gsub(".*TLS_low_", "", contrast))%>%
  write_csv(paste0(pseudobulk_dir, "Compiled_significant_gene_sets_fry.csv"))

GO_fry <- fry_compiled%>%
  filter(grepl("^GO", `Gene set`))

HALLMARK_fry <- fry_compiled%>%
  filter(grepl("^HALLMARK", `Gene set`))

HALLMARK_b <- HALLMARK%>%
  filter(cell_type == "B")

hallmark_b <- gsea_plot(HALLMARK_b, title = "CCR7/SELL CD4+T high vs low (B cells)", top_n = 10)

HALLMARK_mono <- HALLMARK%>%
  filter(cell_type == "Mono")

hallmark_mo <- gsea_plot(HALLMARK_mono, title = "CCR7/SELL CD4+T high vs low (Mono)", top_n = 10)

HALLMARK_Macro <- HALLMARK%>%
  filter(cell_type == "Macro")

hallmark_ma <- gsea_plot(HALLMARK_Macro, title = "CCR7/SELL CD4+T high vs low (Macro)", top_n = 10)

# Make a volcano plot of B cell DE genes
b_cell_de <- read_csv("./Results_paper/analyses/MMR_atlas_TLS/toptables/TLS_high_B-TLS_low_B.csv")
b_cell_de_vol <- volcano_plot(b_cell_de, title = "CCR7/SELL CD4+T high vs low (B cells)", top_n = 10)

b_cell_de_vol <- rasterize(b_cell_de_vol, layers='Point', dpi=100)

# Make a volcano plot of B cell DE genes
EC_de <- read_csv("./Results_paper/analyses/MMR_atlas_TLS/toptables/TLS_high_Endo-TLS_low_Endo.csv")
EC_de_cell_de_vol <- volcano_plot(EC_de, title = "CCR7/SELL CD4+T high vs low (Endo cells)", top_n = 10, )

EC_de_cell_de_vol <- rasterize(EC_de_cell_de_vol, layers='Point', dpi=100)

# Make a plot using TLS high/low subsets
top_left <- plot_grid(de_genes, hallmark_b, labels = c("b", "c"), 
                      label_size = 8, ncol = 1)
tls_toprow <- plot_grid(CT_DA, top_left,rel_widths = c(1.2,1),
                        labels = c("a"), 
                        label_size = 8)

bottom <- plot_grid(b_cell_de_vol,EC_de_cell_de_vol,b_cell_de_vol,EC_de_cell_de_vol,
                                  labels = c("d", "e", "f","g"),nrow = 1, 
                                  label_size = 8)
TLS_2 <- plot_grid(tls_toprow, bottom, nrow = 2, rel_heights = c(1.2,0.8))

ggsave(plot = TLS_2,"./Results_paper/Plots/Figure TLS 2.pdf", width = 170, height = 120, units = "mm")

