
# Set ggplot2 themes for the paper
blank_theme <- theme_bw(base_size = 7)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(size=7))

cell_types <- c("Epi", "Plasma", "Fibro", "Peri", "Macro", "Granulo", "TCD4", "SmoothMuscle", "Endo",
                "B", "DC", "ILC", "Schwann", "Mast", "Mono", "TCD8", "NK", "TZBTB16", "Tgd","QC_fail")

colors <- c("#DEA0FD", "green", "#f99379", "yellowgreen","#654522" ,"#dcd300", "#fdbf6f", "#b2df8a" ,"#CC79A7","#cab2d6",
             "#6a3d9a", "#1f78b4", "#85660D", "#2b3d26","#D55E00", "#a6cee3","darkblue","lightpink", "#1C7F93",
            "grey")

gsea_cols <- c("Up" = "#ff7f00", "Down" = "#a6cee3")

slide_names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
slide_alias <- data.frame(Slide = c("Tumour_A", "Tumour_B", "Normal_A", 
                                    "Run5850.2186.1", "Run5850.2186.2", "10501A", "120131", "148191A", "3461G"),
                          slide_names = slide_names,
                          slide_names_long = c("Slide 1 (Tumour A)", "Slide 2 (Tumour B)", "Slide 3 (Normal A)", "Slide 4 (5850.2186.1)", 
                            "Slide 5 (5850.2186.2)", "Slide 6 (10501A)", "Slide 7 (120131)", 
                            "Slide 8 (148191A)", "Slide 9 (3461G)"))

cell_type_colors <- setNames(colors, cell_types)
print(cell_type_colors)

cell_type_colors[duplicated(cell_type_colors)]

# Set colours for donors and samples
Donor_cols <- c("118" = "#875692", "150" = "#f99379", "73" = "magenta",
                "64" = "lightblue",  "170" = "red", "7" = "blue",
                "142" ="black",  "67" = "darkgreen",
                "278" =  "#40E0D0",
                "NC" = "#708090",
                "130" = "#800000", 
                "280" = "#008080",
                "46" = "#FFA500",
                "7" = "#4B0082")

Sample_cols <- c("118T" = "#875692", "150T" = "#f99379", 
                 "73T" = "magenta","73N" = "purple",
                 "64T" = "lightblue", "64N" = "cyan1", 
                 "170T" = "red", "7N" = "blue","170N" = "pink",
                 "142N" ="black",  "67T" = "darkgreen",
                 "278N" = "#40E0D0",
                 "NC" = "#708090",
                 "130T" = "#800000",
                 "278T" ="#00FF00", 
                 "280T" = "#008080", 
                 "46T" = "#FFA500", 
                 "7T" = "#4B0082")

sample_type_cols <- c("N" = "lightblue", "T" = "orange")

nichecols <- c("Vasculature" = "darkred", 
               "Epithelial mass" = "violet",
               "SMC epithelial rich" = "#875692",
               "Normal like 2" = "darkblue",
               "Macrophage rich" = "#654522",
               "Granulocyte rich" = "#dcd300",
               "Stroma" = "#f99379",
               "Normal like" = "blue",
               "TLS" = "darkgreen",
               "QC_fail" = "grey")
               
nichecols_underscore <- nichecols
names(nichecols_underscore) <- gsub(" ", "_", names(nichecols_underscore))

# Define a function to shut up some other functions
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# Make a pseudobulk heatmap
plot_hm_mmr <- function(genes, bulk_anno, celltypes, cpm){
  
  # Make a heatmap
  bulk_ordered <- bulk_anno%>%
    arrange(Sample_type, Manual_toplevel_pred,Pct_granulo)%>%
    filter(Manual_toplevel_pred %in% celltypes)
  
  genes <- genes[genes%in%rownames(cpm)]
  
  # Scale the data
  scaled <- cpm[genes,bulk_ordered$Sample_cell]%>%
    t()%>%
    scale()%>%
    t()
  
  # mac_mono_cols <- c(Macro = "#654522", Mono = "brown")
  sample_type_cols <- c("N " = "lightgreen", "T MMRp" = "brown", "T MMRd" = "orange")
  
  ha = HeatmapAnnotation(Sample_type = bulk_ordered$Sample_type,
                         Celltype = bulk_ordered$Manual_toplevel_pred,
                         Granulo = anno_barplot(bulk_ordered$Pct_granulo),
                         col = list(Sample_type=sample_type_cols,
                                    Celltype = cell_type_colors))
  
  hm <- Heatmap(scaled, top_annotation = ha, name = "Z score",
                show_column_names = F, 
                cluster_columns = F,
                cluster_rows = F)
  
}

# Function to read the griffith nanostring data
get_nano_counts_FIX <- function(data.dir, sample_name){
  
  slide <- LoadNanostring.FIX(data.dir = data.dir, fov = "fov_all")
  
  slide <- slide[["Nanostring"]]$counts
  
  colnames(slide) <- paste0(sample_name, "_", colnames(slide))
  
  slide <- CreateSeuratObject(slide)
  
}

LoadNanostring.FIX <- function(data.dir, fov, assay = 'Nanostring') {
  data <- ReadNanostring.FIX(
    data.dir = data.dir,
    type = c("centroids", "segmentations")
  )
  segs <- CreateSegmentation(data$segmentations)
  cents <- CreateCentroids(data$centroids)
  segmentations.data <- list(
    "centroids" = cents,
    "segmentation" = segs
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$pixels,
    assay = assay
  )
  obj <- CreateSeuratObject(counts = data$matrix, assay = assay)
  
  # subset both object and coords based on the cells shared by both
  cells <- intersect(
    Cells(x = coords, boundary = "segmentation"),
    Cells(x = coords, boundary = "centroids")
  )
  cells <- intersect(Cells(obj), cells)
  coords <- subset(x = coords, cells = cells)
  obj[[fov]] <- coords
  return(obj)
}

# Get the expression data out of the seurat nanostring object
get_nano_counts <- function(data.dir, sample_name){
  
  slide <- LoadNanostring(data.dir = data.dir, fov = "fov_all")
  
  slide <- slide[["Nanostring"]]$counts
  
  colnames(slide) <- paste0(sample_name, "_", colnames(slide))
  
  slide <- CreateSeuratObject(slide)
  
}

# Read nanostring functions but updated by Lochlan Fennell to work with the new formats
ReadNanostring.FIX <- function(
    data.dir,
    mtx.file = NULL,
    metadata.file = NULL,
    molecules.file = NULL,
    segmentations.file = NULL,
    type = 'centroids',
    mol.type = 'pixels',
    metadata = NULL,
    mols.filter = NA_character_,
    genes.filter = NA_character_,
    fov.filter = NULL,
    subset.counts.matrix = NULL,
    cell.mols.only = TRUE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('centroids', 'segmentations'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels'),
    several.ok = TRUE
  )
  if (!is.null(metadata)) {
    metadata <- match.arg(
      arg = metadata,
      choices = c(
        "Area", "fov", "Mean.MembraneStain", "Mean.DAPI", "Mean.G",
        "Mean.Y", "Mean.R", "Max.MembraneStain", "Max.DAPI", "Max.G",
        "Max.Y", "Max.R"
      ),
      several.ok = TRUE
    )
  }
  
  use.dir <- all(vapply(
    X = c(mtx.file, metadata.file, molecules.file),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  
  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Nanostring directory ", data.dir)
  }
  # Identify input files
  files <- c(
    matrix = mtx.file %||% '[_a-zA-Z0-9]*_exprMat_file.csv',
    metadata.file = metadata.file %||% '[_a-zA-Z0-9]*_metadata_file.csv',
    molecules.file = molecules.file %||% '[_a-zA-Z0-9]*_tx_file.csv',
    segmentations.file = segmentations.file %||% '[_a-zA-Z0-9]*-polygons.csv'
  )
  
  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  
  if (all(is.na(x = files))) {
    stop("Cannot find Nanostring input files in ", data.dir)
  }
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['metadata.file']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    md <- data.table::fread(
      file = files[['metadata.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
    )
    
    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      md <- md[md$fov %in% fov.filter,]
    }
    pprecoord(type = 'finish')
  }
  if (!is.na(x = files[['segmentations.file']])) {
    ppresegs <- progressor()
    ppresegs(
      message = "Preloading cell segmentation vertices",
      class = 'sticky',
      amount = 0
    )
    segs <- data.table::fread(
      file = files[['segmentations.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
    )
    
    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      segs <- segs[segs$fov %in% fov.filter,]
    }
    ppresegs(type = 'finish')
  }
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules.file']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    mx <- data.table::fread(
      file = files[['molecules.file']],
      sep = ',',
      verbose = FALSE
    )
    
    # filter molecules file by FOVs
    if (!is.null(x = fov.filter)) {
      mx <- mx[mx$fov %in% fov.filter,]
    }
    
    # Molecules outside of a cell have a cell_ID of 0
    if (cell.mols.only) {
      mx <- mx[mx$cell_ID != 0,]
    }
    
    if (!is.na(x = mols.filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", mols.filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = mols.filter, x = mx$target), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules.file']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules.file')]
  }
  files <- files[!is.na(x = files)]
  
  outs <- list("matrix"=NULL, "pixels"=NULL, "centroids"=NULL)
  if (!is.null(metadata)) {
    outs <- append(outs, list("metadata" = NULL))
  }
  if ("segmentations" %in% type) {
    outs <- append(outs, list("segmentations" = NULL))
  }
  
  for (otype in names(x = outs)) {
    outs[[otype]] <- switch(
      EXPR = otype,
      'matrix' = {
        ptx <- progressor()
        ptx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        if (!is.null(subset.counts.matrix)) {
          tx <- build.cellcomp.matrix(mols.df=mx, class=subset.counts.matrix)
        } else {
          tx <- data.table::fread(
            file = files[[otype]],
            sep = ',',
            data.table = FALSE,
            verbose = FALSE
          )
          # Combination of Cell ID (for non-zero cell_IDs) and FOV are assumed to be unique. Used to create barcodes / rownames.
          bcs <- paste0(as.character(tx$cell_ID), "_", tx$fov)
          rownames(x = tx) <- bcs
          # remove all rows which represent counts of mols not assigned to a cell for each FOV
          tx <- tx[!tx$cell_ID == 0,]
          # filter fovs from counts matrix
          if (!is.null(x = fov.filter)) {
            tx <- tx[tx$fov %in% fov.filter,]
          }
          # atomx export has an additional column called 'cell' Simply remove it!
          # This casuses Seurat function to fail
          # tx <- subset(tx, select = -c(fov, cell_ID))
          tx <- subset(tx, select = -c(fov, cell, cell_ID))
        }
        
        tx <- as.data.frame(t(x = as.matrix(x = tx)))
        if (!is.na(x = genes.filter)) {
          ptx(
            message = paste("Filtering genes with pattern", genes.filter),
            class = 'sticky',
            amount = 0
          )
          tx <- tx[!grepl(pattern = genes.filter, x = rownames(x = tx)), , drop = FALSE]
        }
        # only keep cells with counts greater than 0
        tx <- tx[, which(colSums(tx) != 0)]
        ratio <- getOption(x = 'Seurat.input.sparse_ratio', default = 0.4)
        
        if ((sum(tx == 0) / length(x = tx)) > ratio) {
          ptx(
            message = 'Converting counts to sparse matrix',
            class = 'sticky',
            amount = 0
          )
          tx <- as.sparse(x = tx)
        }
        
        ptx(type = 'finish')
        
        tx
      },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Creating centroid coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = md$CenterX_global_px,
          y = md$CenterY_global_px,
          cell = paste0(as.character(md$cell_ID), "_", md$fov),
          stringsAsFactors = FALSE
        )
      },
      'segmentations' = {
        pcents <- progressor()
        pcents(
          message = 'Creating segmentation coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = segs$x_global_px,
          y = segs$y_global_px,
          cell = paste0(as.character(segs$cellID), "_", segs$fov),  # cell_ID column in this file doesn't have an underscore
          stringsAsFactors = FALSE
        )
      },
      'metadata' = {
        pmeta <- progressor()
        pmeta(
          message = 'Loading metadata',
          class = 'sticky',
          amount = 0
        )
        pmeta(type = 'finish')
        df <- md[,metadata]
        df$cell <- paste0(as.character(md$cell_ID), "_", md$fov)
        df
      },
      'pixels' = {
        ppixels <- progressor()
        ppixels(
          message = 'Creating pixel-level molecule coordinates',
          class = 'sticky',
          amount = 0
        )
        df <- data.frame(
          x = mx$x_global_px,
          y = mx$y_global_px,
          gene = mx$target,
          stringsAsFactors = FALSE
        )
        ppixels(type = 'finish')
        df
      },
      # 'microns' = {
      #   pmicrons <- progressor()
      #   pmicrons(
      #     message = "Creating micron-level molecule coordinates",
      #     class = 'sticky',
      #     amount = 0
      #   )
      #   df <- data.frame(
      #     x = mx$global_x,
      #     y = mx$global_y,
      #     gene = mx$gene,
      #     stringsAsFactors = FALSE
      #   )
      #   pmicrons(type = 'finish')
      #   df
      # },
      stop("Unknown Nanostring input type: ", outs[[otype]])
    )
  }
  return(outs)
}

# Make a heatmap for the niches
plot_hm_niches <- function(genes, bulk_anno, celltypes, cpm, nichecols, cellcols, samplecols, niches_keep,donor_keep){

  bulk_ordered <- bulk_anno%>%
    arrange(niches,Sample, Manual_toplevel_pred)%>%
    filter(Manual_toplevel_pred %in% celltypes)%>%
    filter(niches %in% niches_keep)%>%
    filter(Donor %in% donor_keep)
  
  # Scale the data
  scaled <- cpm[genes,bulk_ordered$Niche_sample_cell]%>%
    t()%>%
    scale()%>%
    t()
  
  ha = HeatmapAnnotation(Niche = bulk_ordered$niches,
                         Celltype = bulk_ordered$Manual_toplevel_pred,
                         col = list(Niche = nichecols,
                                    Celltype = cellcols))
  
  hm <- Heatmap(scaled, top_annotation = ha, name = "Z score",
                show_column_names = F, 
                cluster_columns = F,
                cluster_rows = F,
                column_split = bulk_ordered$Sample_name)
  
}

# A function to get the niche counts out of the Seurat objects
grab_nichecounts <- function(nanostring_object, slide_name, md_all, k_param){
  
  md_join <- nanostring_object@meta.data%>%
    rownames_to_column("Barcode_old")%>%
    mutate(Barcode = paste0(slide_name, Barcode_old))%>%
    select(-orig.ident)%>%
    left_join(md_all)
  
  nanostring_object$Manual_toplevel_pred <- md_join$Manual_toplevel_pred
  nanostring_object$Barcode <- md_join$Barcode
  nanostring_object$fov <- md_join$fov
  nanostring_object$Sample <- md_join$Sample
  
  # Drop missing cells that were filtered out during clustering
  nanostring_object <- nanostring_object[,!is.na(nanostring_object$Manual_toplevel_pred)]
  
  # Steal from the Seurat BuildNicheAssay function to get nearest neighbours
  
  coords <- GetTissueCoordinates(nanostring_object[["fov_all"]], which = "centroids")
  cells <- coords$cell
  rownames(coords) <- cells
  coords <- as.matrix(coords[, c("x", "y")])
  neighbors <- FindNeighbors(coords, k.param = k_param)
  ct.mtx <- matrix(data = 0, nrow = length(cells), ncol = length(unlist(unique(nanostring_object[["Manual_toplevel_pred"]]))))
  rownames(ct.mtx) <- cells
  colnames(ct.mtx) <- unique(unlist(nanostring_object[["Manual_toplevel_pred"]]))
  cts <- nanostring_object[["Manual_toplevel_pred"]]
  for (i in 1:length(cells)) {
    ct <- as.character(cts[cells[[i]], ])
    ct.mtx[cells[[i]], ct] <- 1
  }
  niche_counts <- as.matrix(neighbors$nn %*% ct.mtx)%>%t()
  niche_counts[1:5,1:5]
  
  # Add the slide name onto the niche counts
  colnames(niche_counts) <- paste0(slide_name, colnames(niche_counts))
  
  return(list(nanostring_object, niche_counts))
  
}

# Function to run GSEA for a contrast
run_GSEA <- function(contrast , collection, rnaseq, v, design, cont.matrix, pseudobulk_dir){
  
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
  
  fry_result <- fry(y = v ,index = indexed, design = design, contrast = cont.matrix[,contrast])%>%
    rownames_to_column("Gene set")%>%
    dplyr::select(`Gene set`,"NGenes" , "Direction", "PValue", "FDR")%>%
    filter(FDR <= 0.05)%>%
    mutate(Contrast= contrast)
  
  # Make a directory
  system(paste0("mkdir -p ", pseudobulk_dir,"gsea/camera/"))
  system(paste0("mkdir -p ", pseudobulk_dir,"gsea/fry/"))
  
  write_csv(camera_result, paste0(pseudobulk_dir,"gsea/camera/",collection_name, "_", contrast,".csv"))
  write_csv(fry_result, paste0(pseudobulk_dir,"gsea/fry/",collection_name, "_", contrast,".csv"))
  
}

# A function to get the niche counts out of the Seurat objects
annotate_obj <- function(nanostring_object, slide_name, md_all){
  
  # The full unfiltered nanostring metadata
  md_all_join<- md_all %>%
    select(Barcode,fov, Sample)
  
  md_join <- nanostring_object@meta.data%>%
    rownames_to_column("Barcode_old")%>%
    mutate(Barcode = paste0(slide_name, Barcode_old))%>%
    select(-orig.ident)%>%
    left_join(md_all)%>%
    select(-Sample, -fov)%>%
    left_join(md_all_join)
  
  rownames(md_join) <- md_join$Barcode_old
  nanostring_object@meta.data <- md_join
  
  # Drop missing cells that were filtered out during clustering
  nanostring_object$Manual_toplevel_pred <- replace(nanostring_object$Manual_toplevel_pred, is.na(nanostring_object$Manual_toplevel_pred), "QC_fail")
  nanostring_object$niches <- replace(nanostring_object$niches, is.na(nanostring_object$niches), "QC_fail")
  
  return(nanostring_object)
  
}

# Plot some molecules for a given celltype/FOV
plot_molecules <- function(niche_match, nanostring_obj, sample_name, FOV, mycols, mols, outdir, cell_types = "all", plot_name = "plot"){
  
  system(paste0("mkdir -p ", outdir, sample_name))
  
  # Match the sample barcode to the niche data frame
  matched <- match(nanostring_obj$Barcode, niche_match$Barcode)
  
  # Get the correct rows of the dataframe containing the niches
  niche_df_nanostring <- niche_match[matched,]
  
  # Add the correctly ordered niches to the nanostring object
  nanostring_obj$niches <- niche_df_nanostring$niches
  
  nanostring_obj$niches<- factor(nanostring_obj$niches, levels = unique(niche_match$niches))
  
  # Set the defualt boundry to be segmentation
  DefaultBoundary(nanostring_obj[["fov_all"]]) <- "segmentation"
  
  # Subset the full object down to just that FOV
  fov_1_cells <- colnames(nanostring_obj)[nanostring_obj$fov == FOV]
  
  # Get the sample name
  sample <- nanostring_obj$Sample[nanostring_obj$fov == FOV]%>%
    unique()
  
  # Pull out the coordinates for that object
  coords <- GetTissueCoordinates(nanostring_obj)
  coords_fov <- coords%>%
    filter(cell %in% fov_1_cells)
  
  # Add an FOV to as a zoom back onto the original data
  crop <- Crop(nanostring_obj[["fov_all"]], y = c(max(coords_fov$x), min(coords_fov$x)), 
               x = c(max(coords_fov$y), min(coords_fov$y))
  )
  nanostring_obj[["zoom1"]] <- crop
  DefaultBoundary(nanostring_obj[["zoom1"]]) <- "segmentation"
  
  p1 <- ImageDimPlot(nanostring_obj, fov = "zoom1", cols = mycols,
                     group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                     mols.size = 0.3, border.color = "black", coord.fixed = FALSE)+
    coord_flip()+
    theme(aspect.ratio = 1)+
    ggtitle(paste0(sample_name, " ", sample, " FOV ", FOV))
  
  # Figure out the cells to plot based on the cell types
  if(cell_types[1] != "all"){
    cells_plot <- colnames(nanostring_obj)[nanostring_obj$Manual_toplevel_pred %in% cell_types] 
    
    # Plot the FOV
    p2 <- ImageDimPlot(nanostring_obj, fov = "zoom1", 
                       cols = mycols,
                       group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                       cells = cells_plot, 
                       molecules = mols,
                       mols.size = 0.3, 
                       border.color = "black", coord.fixed = FALSE,
                       nmols = 60000)+
      coord_flip()+
      theme(aspect.ratio = 1)
  }else{
    # Plot the FOV
    p2 <- ImageDimPlot(nanostring_obj, fov = "zoom1", 
                       #cols = rep("grey", length(unique(nanostring_obj$Manual_toplevel_pred))),
                       cols = mycols,
                       group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                       molecules = mols,
                       mols.size = 0.3, 
                       border.color = "black", coord.fixed = FALSE,
                       nmols = 60000)+
      coord_flip()+
      theme(aspect.ratio = 1)+
      guides(fill = "none")
  }
  
  plot_save <- paste0(outdir, sample_name, "/moecules_FOV_", FOV,"_",plot_name, ".pdf")
  
  pdf(plot_save, width = 16, height = 9)
  print(p1 | p2)
  dev.off()
  
}

# Cache the FOVs for quick plotting
cache_FOV <- function(FOV, niche_match, nanostring_obj, slide_name, outdir){
  
  system(paste0("mkdir -p ", outdir, slide_name))
  
  # Match the sample barcode to the niche data frame
  matched <- match(nanostring_obj$Barcode, niche_match$Barcode)
  
  # Get the correct rows of the dataframe containing the niches
  niche_df_nanostring <- niche_match[matched,]
  
  # Add the correctly ordered niches to the nanostring object
  nanostring_obj$niches <- niche_df_nanostring$niches
  
  # Set the defualt boundry to be segmentation
  DefaultBoundary(nanostring_obj[["fov_all"]]) <- "segmentation"
  
  # Subset the full object down to just that FOV
  fov_1_cells <- colnames(nanostring_obj)[nanostring_obj$fov == FOV]
  
  # Get the sample name
  sample <- nanostring_obj$Sample[nanostring_obj$fov == FOV]%>%
    unique()
  
  nanostring_obj <- nanostring_obj[,colnames(nanostring_obj) %in% fov_1_cells]
  
  # Pull out the coordinates for that object
  coords <- GetTissueCoordinates(nanostring_obj)
  coords_fov <- coords%>%
    filter(cell %in% fov_1_cells)
  
  # Add an FOV as a zoom back onto the original data
  crop <- Crop(nanostring_obj[["fov_all"]], y = c(max(coords_fov$x), min(coords_fov$x)), 
               x = c(max(coords_fov$y), min(coords_fov$y))
  )
  nanostring_obj[["zoom1"]] <- crop
  DefaultBoundary(nanostring_obj[["zoom1"]]) <- "segmentation"
  
  nanostring_obj$niches <- as.character(nanostring_obj$niches)
  nanostring_obj$niches <- replace(nanostring_obj$niches, is.na(nanostring_obj$niches), "QC_fail")
  
  # Make the object a lot smaller
  nanostring_obj[['fov_all']] <- NULL
  
  # Qsave the object
  qsave(nanostring_obj, paste0(outdir, slide_name,"/FOV_", FOV, ".qs"))
  
  return("Done")
  
}

plot_gene_set <- function(niche_match, nanostring_obj, sample_name, FOV, mycols, gene_set, outdir, cell_types = "all", plot_name = "plot"){
  
  system(paste0("mkdir -p ", outdir, sample_name))
  
  # Match the sample barcode to the niche data frame
  matched <- match(nanostring_obj$Barcode, niche_match$Barcode)
  
  # Get the correct rows of the dataframe containing the niches
  niche_df_nanostring <- niche_match[matched,]
  
  # Add the correctly ordered niches to the nanostring object
  nanostring_obj$niches <- niche_df_nanostring$niches
  
  nanostring_obj$niches<- factor(nanostring_obj$niches, levels = unique(niche_match$niches))
  
  # Set the defualt boundry to be segmentation
  DefaultBoundary(nanostring_obj[["fov_all"]]) <- "segmentation"
  
  # Subset the full object down to just that FOV
  fov_1_cells <- colnames(nanostring_obj)[nanostring_obj$fov == FOV]
  
  # Get the sample name
  sample <- nanostring_obj$Sample[nanostring_obj$fov == FOV]%>%
    unique()
  
  # Pull out the coordinates for that object
  coords <- GetTissueCoordinates(nanostring_obj)
  coords_fov <- coords%>%
    filter(cell %in% fov_1_cells)
  
  # Add an FOV to as a zoom back onto the original data
  crop <- Crop(nanostring_obj[["fov_all"]], y = c(max(coords_fov$x), min(coords_fov$x)), 
               x = c(max(coords_fov$y), min(coords_fov$y))
  )
  nanostring_obj[["zoom1"]] <- crop
  DefaultBoundary(nanostring_obj[["zoom1"]]) <- "segmentation"
  
  p1 <- ImageDimPlot(nanostring_obj, fov = "zoom1", cols = mycols,
                     group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                     mols.size = 0.3, border.color = "black", coord.fixed = FALSE)+
    coord_flip()+
    theme(aspect.ratio = 1)+
    ggtitle(paste0(sample_name, " ", sample, " FOV ", FOV))
  
  p2 <- ImageFeaturePlot(nanostring_obj, fov = "zoom1", 
                         features = gene_set,
                         max.cutoff = "q95",
                         boundaries = "segmentation",
                         border.color = "black", 
                         coord.fixed = FALSE)+
    coord_flip()+
    theme(aspect.ratio = 1)
  
  p2 <- ImageDimPlot(nanostring_obj, fov = "zoom1", 
                     cols = mycols,
                     group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                     cells = cells_plot, 
                     molecules = mols,
                     mols.size = 0.3, 
                     border.color = "black", coord.fixed = FALSE,
                     nmols = 60000)+
    coord_flip()+
    theme(aspect.ratio = 1)
  
  # Figure out the cells to plot based on the cell types
  if(cell_types[1] != "all"){
    cells_plot <- colnames(nanostring_obj)[nanostring_obj$Manual_toplevel_pred %in% cell_types] 
    
    # Plot the FOV
    p2 <- ImageDimPlot(nanostring_obj, fov = "zoom1", 
                       cols = mycols,
                       group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                       cells = cells_plot, 
                       molecules = mols,
                       mols.size = 0.3, 
                       border.color = "black", coord.fixed = FALSE,
                       nmols = 60000)+
      coord_flip()+
      theme(aspect.ratio = 1)
  }else{
    # Plot the FOV
    p2 <- ImageDimPlot(nanostring_obj, fov = "zoom1", 
                       #cols = rep("grey", length(unique(nanostring_obj$Manual_toplevel_pred))),
                       cols = mycols,
                       group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                       molecules = mols,
                       mols.size = 0.3, 
                       border.color = "black", coord.fixed = FALSE,
                       nmols = 60000)+
      coord_flip()+
      theme(aspect.ratio = 1)+
      guides(fill = "none")
  }
  
  plot_save <- paste0(outdir, sample_name, "/moecules_FOV_", FOV,"_",plot_name, ".pdf")
  
  pdf(plot_save, width = 19, height = 9)
  print(p1 | p2)
  dev.off()
  
}

# Plot gene T vs N
t_vs_n_gene <- function(cpm, gene, bulk_anno, outdir){
  
  # Make a plot of a single gene
  single_gene <- cpm[gene,]%>%
    data.frame(check.rows = F)%>%
    rownames_to_column("Sample_cell")%>%
    dplyr::rename(CPM = 2)%>%
    left_join(bulk_anno)%>%
    group_by(Manual_toplevel_pred, Sample_type)
  
  plt <- ggplot(data = single_gene, aes(x = Sample_type, y = CPM, group = donor, colour = donor))+ 
    geom_point(alpha = 0.5)+
    geom_line(alpha = 0.5)+
    facet_wrap(~Manual_toplevel_pred)+
    guides(colour ="none")+
    blank_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Sample type", y = "Log2 CPM", fill = "Sample", title = gene)
  
  ggsave(paste0(outdir, gene, " cell type T vs N.pdf"), width = 5, height = 3)
  
  print(plt)
  
}

# Get the pseudobulk counts
get_pb <- function(i, to_bulk_anno, aggr_col){
  
  d_t <- unique(to_bulk_anno[,aggr_col])[i]
  
  # Filter the annotation to keep only cells from a certain sample/cluster
  anno_filt <- to_bulk_anno[to_bulk_anno[,aggr_col] == d_t,]
  
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

# Get only the nuclear reads for each cell
get_nuclear_reads <- function(location_file){
  
  slide <- basename(dirname(location_file))
  
  # Transcripts not assigned to cells have a value of 0
  tx_file <- read_csv(location_file)%>%
    filter(CellComp == "Nuclear")%>%
    filter(cell_ID != 0)%>%
    mutate(Barcode = paste0(slide, "_", cell_ID, "_", fov))%>%
    group_by(Barcode, target)%>%
    summarise(Gene_count = n())%>%
    ungroup()%>%
    pivot_wider(names_from = Barcode, values_from = Gene_count)
  
  # Reprocess the couts out of the TX_file
  head(tx_file)
  
  tx_mat <- as.matrix(tx_file[,2:ncol(tx_file)])
  
  rownames(tx_mat) <- tx_file$target
  
  tx_mat[is.na(tx_mat)] <- 0
  
  return(tx_mat)
  
}

# Plot a gene in a pseudobulk dataset
plot_gene <- function(cpm, gene, bulk_anno, outdir){
  
  # Make a plot of a single gene
  single_gene <- cpm[gene,]%>%
    data.frame(check.rows = F)%>%
    rownames_to_column("Niche_sample_cell")%>%
    dplyr::rename(CPM = 2)%>%
    left_join(bulk_anno)%>%
    group_by(Manual_toplevel_pred, niches)%>%
    summarise(mean_cpm = mean(CPM))%>%
    ungroup()%>%
    arrange(-mean_cpm)%>%
    mutate(Manual_toplevel_pred = factor(Manual_toplevel_pred, levels = unique(Manual_toplevel_pred)))
  
  plt <- ggplot(data = single_gene, aes(x = Manual_toplevel_pred, y = mean_cpm, fill = niches))+ 
    geom_bar(stat = "identity", position = "dodge")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = nichecols)+
    labs(x = "Cell type", y = "Mean log2 CPM", fill = "Niche", title = gene)
  
  ggsave(paste0(outdir, gene, " cell type and niche expression.pdf"), width = 9, height = 5)
  
  print(plt)
  
}

# Plot a gene across niches
niche_vs_niche_gene <- function(cpm, gene, bulk_anno, outdir, ct_filter = NULL){
  
  # Make a plot of a single gene
  single_gene <- cpm[gene,]%>%
    data.frame(check.rows = F)%>%
    rownames_to_column("Niche_sample_cell")%>%
    dplyr::rename(CPM = 2)%>%
    left_join(bulk_anno)%>%
    group_by(Manual_toplevel_pred, niches)%>%
    mutate(niches = gsub("_", " ", niches))
  
  if(!is.null(ct_filter)){
    
    single_gene <- single_gene %>%
      filter(Manual_toplevel_pred %in% ct_filter)
    
  }
  
  plt <- ggplot(data = single_gene, aes(x = niches, y = CPM, group = Donor, colour = niches))+ 
    geom_point(alpha = 0.8)+
    stat_summary(fun = "median", fun.min = "median", fun.max= "median", linewidth= 0.3, geom = "crossbar", 
                 aes(group = niches, colour = niches))+
    scale_colour_manual(values = nichecols)+
    geom_line(alpha = 0.1)+
    facet_wrap(~Manual_toplevel_pred)+
    guides(colour ="none")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Niche", y = "Log2 CPM", fill = "Sample", title = gene)+
    blank_theme
  
  ggsave(paste0(outdir, gene, " cell type T vs N.pdf"), width = 9, height = 5)
  
  print(plt)
  
}


# Plt the expression ofa gene in all cellypes
plot_gene_celltype <- function(cpm, gene, bulk_anno, outdir, axes = T){
  
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
    geom_jitter(width = 0.1, height = 0, aes(colour = Manual_toplevel_pred), size = 1)+
    facet_wrap(~Sample_type)+
    guides(colour ="none")+
    blank_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Cell type", y = expression('Mean log'[2]*' CPM'), fill = "Sample")+
    stat_summary(fun = "median", fun.min = "median", fun.max= "median", linewidth= 0.3, geom = "crossbar")+
    scale_colour_manual(values = cell_type_colors)+
    ggtitle(gene)+
    theme(plot.title = element_text(face = "italic"))+
    coord_flip()
  
  ggsave(paste0(outdir, gene, " cell type expression.pdf"), width = 9, height = 5)
  
  if(axes == F){
    plt <- plt+labs(x = NULL, y = NULL)
  }
  
  print(plt)
  
}

# This functions finds FOVs where a gene is highly expressed
gene_fov_finder <- function(gene, seurat_obj){
  
  gene_expr <- seurat_obj@assays$RNA$counts[gene,]
  seurat_obj$Gene <- gene_expr
  gene_plot <- seurat_obj@meta.data%>%
    group_by(Sample, Slide, fov, Manual_toplevel_pred)%>%
    summarise(Total = sum(Gene))%>%
    mutate(Gene = gene)%>%
    ungroup()%>%
    arrange(-Total)
  
  return(gene_plot)
  
}

# List the order of gene expression for a given slide.
slide_list <- function(object, slide, FOV){
  
  # Get the right cells
  slide_fov <- object[,object$Slide == slide & object$fov == FOV]
  
  # Return expression
  gene_summary <- slide_fov@assays$RNA$counts
  
  gene_mat <- data.frame(gene_summary)%>%
    rownames_to_column("Gene")%>%
    gather(Barcode, Count, -Gene)%>%
    left_join(slide_fov@meta.data)%>%
    group_by(Gene, niches, Manual_toplevel_pred)%>%
    summarise(Total_probecount = sum(Count))%>%
    arrange(-Total_probecount)
  
}

# Make a volcano plot of limma toptables
volcano_plot <- function(toptable, title, top_n =20, text_size = 2.5){
  
  tt <- toptable %>%
    mutate(FC = ifelse(logFC > 0, "Up (FDR < 0.05)", "Down (FDR < 0.05)"))%>%
    mutate(FC = replace(FC, adj.P.Val >= 0.05, "FDR > 0.05"))%>%
    mutate(FC = factor(FC, levels = c("FDR > 0.05", "Down (FDR < 0.05)","Up (FDR < 0.05)")))%>%
    mutate(Count = 1:n())%>%
    mutate(label = replace(SYMBOL, Count >top_n, NA))
  
  max <- abs(toptable$logFC)%>%
    max()
  min <- max *-1
  
  ggplot(data = tt, aes(x = logFC, y = -log10(adj.P.Val), label = label,  colour= FC))+
    geom_point(size = 0.5)+
    blank_theme+
    theme(legend.key.size = unit(4, 'mm'))+
    geom_text_repel(max.overlaps =40,size=text_size, colour = "black", aes(fontface=3))+
    geom_vline(xintercept = 0,linetype = "dashed")+
    scale_colour_manual(values = c("grey", "blue", "red"), drop = F)+
    guides(alpha = "none", colour = "none", size = "none")+
    labs(x = expression('Log'[2]*' fold change'), y = expression('-Log'[10]*' FDR'), colour = "Significance")+
    ggtitle(title)+
    scale_x_continuous(limit = c(min, max))
  
}

# Plot gene set enrichment analysis results
gsea_plot <- function(gsea_table, top_n, title, xlab = "Hallmark gene set"){
  
  gsea_tab <- gsea_table%>%
    mutate(count = 1:n())%>%
    filter(count <= top_n)%>%
    mutate(dir = ifelse(Direction == "Up", 1, -1))%>%
    mutate(negl10 = -log10(FDR))%>%
    mutate(negl10 = negl10 * dir)%>%
    mutate(`Gene set` = gsub("_", " ",`Gene set`))%>%
    mutate(`Gene set` = gsub("HALLMARK", "",`Gene set`))%>%
    mutate(`Gene set` = factor(`Gene set`, levels = unique(`Gene set`)))
  
  gsea_1 <- ggplot(data = gsea_tab, aes(x = `Gene set`, y = negl10, fill = Direction))+
    geom_bar(stat = "identity")+
    coord_flip()+
    blank_theme+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(4, 'mm'))+
    scale_fill_manual(values = gsea_cols)+
    labs(title = title, x = xlab, y = expression('-log'[10]*' FDR'),
         fill = "Direction\n")
  
}
