#shinyOptions(cache = cachem::cache_disk("/homevol/shiny/app_cache/cache/"))
# Secure shiny login
library(shinymanager)
library(qs)
library(tidyverse)
library(Seurat)
# Shiny themes
library(bslib)
# Extra shiny UI widgets
library(shinyWidgets)
library(shiny)
library(DT)

toplevel_dir <- "/pvol/shiny_data/cached_fovs_smaller/"

donor_stats <- read_csv("/pvol/shiny_data/Slide-FOV sample metadata.csv")

# List the samples
#slides <- list.files(toplevel_dir)

slides <- c("Tumor_A", "Tumor_B", "Normal_A", 
            "Run5850.2186.1", "Run5850.2186.2", "10501A_", "120131_", "148191A_", "3461G_")
names(slides) <- c("Slide 1 (Tumour A)", "Slide 2 (Tumour B)", "Slide 3 (Normal A)", "Slide 4 (5850.2186.1)", 
                   "Slide 5 (5850.2186.2)", "Slide 6 (10501A)", "Slide 7 (120131)", 
                   "Slide 8 (148191A)", "Slide 9 (3461G)")


cell_types <- c("Epi", "Plasma", "Fibro", "Peri", "Macro", "Granulo", "TCD4", "SmoothMuscle", "Endo",
                "B", "DC", "ILC", "Schwann", "Mast", "Mono", "TCD8", "NK", "TZBTB16", "Tgd","QC_fail")

colors <- c("#DEA0FD", "green", "#f99379", "yellowgreen","#654522" ,"#dcd300", "#fdbf6f", "#b2df8a" ,"#CC79A7","#cab2d6",
            "#6a3d9a", "#1f78b4", "#85660D", "#2b3d26","#D55E00", "#a6cee3","darkblue","lightpink", "#1C7F93",
            "grey")

mycols <- setNames(colors, cell_types)

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

overlay_columns <- c("niches", 
                     "nCount_Nanostring",
                     "nFeature_Nanostring","Area",
                     "Mean.MembraneStain",
                     "Mean.MembraneStain", "Max.MembraneStain", "Mean.PanCK",                 
                     "Max.PanCK", "Mean.CD45", "Max.CD45", "Mean.CD3",
                     "Mean.CD68", "Max.CD68",
                     "Mean.MembraneStain_B2M",     
                     "Max.MembraneStain_B2M",
                     "Max.CD3", "Mean.DAPI", "Max.DAPI",
                     "IFNG_genes", "MHC_2", 
                     "PrimaryCellAtlas_pred",
                     "Main_PrimaryCellAtlas_pred",
                     "MMR_pred",
                     "Main_MMR_Atlas_pred", "predicted.celltype.l1",
                     "predicted.celltype.l2",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                     "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                     "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
                     "Adult_stem_vs_Adult_Enterocyte_up",
                     "Adult_stem_vs_Adult_Enterocyte_down",
                     "Second_trim_Enterocyte_vs_Adult_Enterocyte_up",
                     "Second_trim_Enterocyte_vs_Adult_Enterocyte_down",
                     "GOBP_NEUTROPHIL_CHEMOTAXIS",
                     "HALLMARK_G2M_CHECKPOINT"
                     )

blank_theme <- theme_bw(base_size = 14)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

# The Seurat plotting code sometimes needs this for some reason
options(future.globals.maxSize = 8000 * 1024^2)

# Shiny UI (defines how the app looks) ----
ui <- fluidPage(
  
  # Set the theme for the app
  theme = bs_theme(version = 5, bootswatch = "vapor"),
  
  titlePanel("cosmxos (the x is silent)"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput("slide", label ="Select a slide", 
                  choices = slides, 
                  width = "100%",
                  multiple = F,
                  selected = 1),
      
      selectInput("FOV", "Select a FOV to plot", 
                  width = "100%",
                  choices = NULL,
                  multiple = F),
      
      selectInput("niche_overlay", "Select data type to plot", 
                  width = "100%",
                  choices = overlay_columns,
                  multiple = F),
      
      pickerInput("cell_types","Select cell types to plot (overlay and molecule plots)", 
                  choices=NULL, 
                  options = list(`actions-box` = TRUE),
                  multiple = T,
                  width = "100%"
      ),
      selectizeInput("mols", label ="Select genes to plot on the molecule plot", 
                     choices = NULL, 
                     width = "100%",
                     multiple = T),
      checkboxInput("legend", label = "Show legend on plot",value = T),
      checkboxInput("title", label = "Show title on plot",value = T),
      h5("Slide guide"),
      tags$iframe(style="height:600px; width:100%", src="pics/FOV_niches_long_2.pdf"),
      h5("Sample quailty guide"),
      tags$iframe(style="height:600px; width:100%", src="pics/Silhouette plot.pdf")
    ),
    mainPanel(
      # The plot of each cell type
      plotOutput("celltype_plot_print", width = "750px", height = "650px"),
      fluidRow(
        column(2, numericInput("celltype_plot_width", label = "Plot download width (mm)", value = 250)),
        column(2, numericInput("celltype_plot_height", label = "Plot download height (mm)", value = 200)),
        column(2, downloadButton('celltype_plot_download', label = "Download celltype plot (pdf)")),
        column(2, downloadButton('celltype_plot_download_png', label = "Download celltype plot (png)"))
      ),
      plotOutput("niche_plot_print", width = "750px", height = "650px"),
      fluidRow(
        column(2, numericInput("niche_plot_width", label = "Plot download width (mm)", value = 250)),
        column(2, numericInput("niche_plot_height", label = "Plot download height (mm)", value = 200)),
        column(2, downloadButton('niche_plot_download', label = "Download overlay plot (pdf)")),
        column(2, downloadButton('niche_plot_download_png', label = "Download overlay plot (png)"))
      ),
      plotOutput("molecule_plot_print", width = "750px", height = "650px"),
      fluidRow(
        column(2, numericInput("molecule_plot_width", label = "Plot download width (mm)", value = 250)),
        column(2, numericInput("molecule_plot_height", label = "Plot download height (mm)", value = 200)),
        column(2, downloadButton('molecule_plot_download', label = "Download molecule plot (pdf)")),
        column(2, downloadButton('molecule_plot_download_png', label = "Download molecule plot (png)"))
      ),
      h6("DAPI is blue, PanCK is green, CD45 is red, CD3 is white (first run only)"),
      imageOutput("stained_image", width = "800px", height = "800px"),
      fluidRow(
        column(2, downloadButton('composite_download', label = "Download composite plot (jpg)"))
      ),
      h5("Slide layout guide"),
      tags$iframe(style="height:500px; width:100%", src="pics/Figure S slides.pdf"),
      h5("Sample info"),
      DT::dataTableOutput("table")
    )
  )
)

#Secure the app with shinymanager
ui <- secure_app(ui)

# Shiny server (does the calculations for things that get displayed in the ui) ----
server <- function(input, output, session) {
  
  #check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  output$auth_output <- renderPrint({
    reactiveValuesToList(res_auth)
  })
  
  # Get a lsit of possible FOVs
  get_FOV <- reactive({
    
    req(input$slide)
    
    FOVs <- list.files(paste0(toplevel_dir, input$slide))
    
    FOVs <- gsub("\\.qs|FOV_", "", FOVs)
    
    FOVs[order(as.numeric(FOVs))]
    
  })
  
  # Update the select input with the FOVs that are in the directory
  observe(updateSelectInput(session = session, inputId = 'FOV', choices = get_FOV()))
  
  # Read the FOV and get plotting options
  pre_read_fov <- reactive({
    
    # You need the FOV to plot read the FOV
    req(get_FOV())
    req(input$slide)
    
    print("read FOV")
    
    file_path <- paste0(paste0(toplevel_dir, input$slide, "/FOV_", input$FOV, ".qs"))
    
    print(file_path)
    
    if(file.exists(file_path)){
      print("File Present")
      cached_fov <- qread(file_path)
      cached_fov$Manual_toplevel_pred <- gsub("_", " ", cached_fov$Manual_toplevel_pred)
      return(cached_fov)
    }else{
      print("File Not Present")
      return(NULL)
    }
    
  })
  
  # Use the read in FOV and plot the selected overlay option
  read_FOV <- reactive({
    
    req(pre_read_fov())
    req(input$niche_overlay)
    
    cached_fov <- pre_read_fov()
    
    if( is.numeric(cached_fov@meta.data[,input$niche_overlay])){
      cached_fov@meta.data[,input$niche_overlay] <- replace(cached_fov@meta.data[,input$niche_overlay], is.na(cached_fov@meta.data[,input$niche_overlay]), min(cached_fov@meta.data[,input$niche_overlay], na.rm = T))
    }else{
      cached_fov@meta.data[,input$niche_overlay] <- replace(cached_fov@meta.data[,input$niche_overlay], is.na(cached_fov@meta.data[,input$niche_overlay]), "QC_fail")
    }
    
    return(cached_fov)
    
  })
  
  sample_name <- reactive({
    
    req(pre_read_fov())
    
    print("get sample names")
    
    sample <- unique(pre_read_fov()$Sample)
    
    title <- paste0(input$slide, " ",input$FOV, " ", sample)
    
  })
  
  celltype_plot <- reactive({
    
    req(pre_read_fov())
    req(sample_name())
    
    plt <- ImageDimPlot(pre_read_fov(), fov = "zoom1", cols = mycols,
                        group.by = "Manual_toplevel_pred", boundaries = "segmentation",
                        mols.size = 0.3, border.color = "black", coord.fixed = FALSE)+
      coord_flip()+
      theme(aspect.ratio = 1)+
      labs(fill = "Cell type", title = sample_name())
    
    if(input$legend == F){
      plt <- plt+
        NoLegend()
    }
    
    if(input$title == F){
      plt <- plt+
        labs(title = NULL)
    }
    
    return(plt)
    
  })
  
  # Get a list of possible cell types
  get_cell_types <- reactive({
    
    req(input$slide)
    req(input$FOV)
    req(pre_read_fov())
    
    print("get_cell_types")
    
    unique(pre_read_fov()$Manual_toplevel_pred)
    
  })
  
  # Get a list of possible genes
  get_genes <- reactive({
    
    req(input$slide)
    req(input$FOV)
    req(read_FOV())
    req(input$niche_overlay)
    req(sample_name())
    
    print("get genes")
    
    rownames(read_FOV())
    
  })
  
  # Get the possible cell types from the loaded FOV
  observe(updatePickerInput(session = session, inputId = 'cell_types', choices = get_cell_types(), selected = get_cell_types()))
  
  observe(updateSelectizeInput(session = session, inputId = 'mols', choices = get_genes(), selected = NULL))
  
  molecule_plot <- reactive({
    
    req(input$slide)
    req(input$FOV)
    req(read_FOV())
    req(input$cell_types)
    req(sample_name())
    req(input$mols)
    
    print("make molecule plot")
    
    cells_plot <- colnames(read_FOV())[read_FOV()$Manual_toplevel_pred %in% input$cell_types] 
    
    # Plot the FOV
    plt <- ImageDimPlot(read_FOV(), fov = "zoom1", 
                        cols = mycols,
                        group.by = "Manual_toplevel_pred", 
                        boundaries = "segmentation",
                        cells = cells_plot, 
                        molecules = input$mols,
                        mols.size = 0.2, 
                        border.color = "black",
                        dark.background = T,
                        coord.fixed = FALSE,
                        nmols = 60000)+
      coord_flip()+
      theme(aspect.ratio = 1)+
      ggtitle(sample_name())+
      labs(fill = "Cell type")
    
    if(input$legend == F){
      
      plt <- plt+NoLegend()
    }
    
    if(input$title == F){
      plt <- plt+
        labs(title = NULL)
    }
    return(plt)
  })
  
  niche_plot <- reactive({
    
    req(input$niche_overlay)
    req(input$cell_types)
    req(input$FOV)

    ct <- isolate({input$cell_types})
    no <- isolate({input$niche_overlay})
    sn <- isolate({sample_name()})
    fv <- isolate({read_FOV()})
    
    cells_plot <- colnames(fv)[fv$Manual_toplevel_pred %in% ct] 
    
    if( is.numeric(fv@meta.data[,no])){
      
      # Make a numeric plot of the feature
      plt <- ImageFeaturePlot(fv, fov = "zoom1", 
                              features = no, 
                              cells = cells_plot, 
                              max.cutoff = "q95", size = 1.5, 
                              dark.background = T,
                              border.color = "black")+
        coord_flip()+
        theme(aspect.ratio = 1,
              plot.title = element_text(hjust = 0))+
        labs(title = paste0(sn, " ", no), fill = "")
    }else{
      
      # Plot the FOV
      plt <-  ImageDimPlot(fv, fov = "zoom1",cells = cells_plot, 
                           group.by = no, 
                           size = 1.5, 
                           cols = "glasbey",
                           dark.background = T,
                           border.color = "black")+
        coord_flip()+
        theme(aspect.ratio = 1)+
        labs(title = paste0(sn, " ", no), fill = "")
      
    }
    
    # Use the niche colours if plotting niches
    if(no == "niches"){
      plt <- plt+ scale_fill_manual(values = nichecols)
    }
    
    if(input$legend == F){
      plt <- plt+NoLegend()
    }
    
    if(input$title == F){
      plt <- plt+
        labs(title = NULL)
    }
    
    return(plt)

  })
  
  # Show the matching protein image
  output$stained_image <- renderImage({
    
    req(input$slide)
    req(input$FOV)
    
    outfile <<- list.files (paste0("/oldvol/shiny/composite/", input$slide, "/CellComposite/"), full.names = T)
    
    fov <- gsub(".*F0|_composite_autocontrast", "", outfile)
    fov <- gsub("^0|\\.jpg", "", fov)
    
    outfile <- outfile[input$FOV == fov]
    
    # Return a list containing the filename
    list(src = outfile,
         width =800,
         height = 800,
         contentType = 'image/jpg',
         alt = "This is alternate text")
  }, deleteFile = F
  )
  
  output$celltype_plot_print = renderPlot({
    
    req(celltype_plot())
    
    # Print the plot
    print(celltype_plot())
    
  })
  
  output$molecule_plot_print = renderPlot({
    
    req(molecule_plot())
    
    # Print the plot
    print(molecule_plot())
    
  })
  
  output$niche_plot_print = renderPlot({
    
    req(niche_plot())
    # Print the plot
    print(niche_plot())
    
  })
  
  # Download the cell type plot
  output$celltype_plot_download =  downloadHandler(
    
    filename = 'Celltype plot.pdf',
    content = function(file) {
      
      device <- function(..., width, height) {
        grDevices::pdf(file = file,
                       width = width,
                       height = height)
      }
      
      ggsave(
        file,
        plot = celltype_plot(),
        device = device,
        width = input$celltype_plot_width,
        units = "mm",
        height = input$celltype_plot_height,
        limitsize = FALSE
      )
    }
  )
  
  # Download the tumour plot on click
  output$celltype_plot_download_png = downloadHandler(
    filename = "Celltype plot.png",
    content = function(file) {
      png(file, height = input$celltype_plot_height, 
          width = input$celltype_plot_width, unit = "mm", res = 600)
      print(celltype_plot())
      dev.off()
    }
  )
  
  # Download the niche plot
  output$niche_plot_download =  downloadHandler(
    
    filename = 'Overlay plot.pdf',
    content = function(file) {
      
      device <- function(..., width, height) {
        grDevices::pdf(file = file,
                       width = width,
                       height = height)
      }
      
      ggsave(
        file,
        plot = niche_plot(),
        device = device,
        width = input$niche_plot_width,
        units = "mm",
        height = input$niche_plot_height,
        limitsize = FALSE
      )
    }
  )
  
  # Download the tumour plot on click
  output$niche_plot_download_png = downloadHandler(
    filename = "Overlay plot.png",
    content = function(file) {
      png(file, height = input$niche_plot_height, 
          width = input$niche_plot_width, unit = "mm", res = 600)
      print(niche_plot())
      dev.off()
    }
  )
  
  # Download the molecule plot
  output$molecule_plot_download =  downloadHandler(
    
    filename = 'Molecule plot.pdf',
    content = function(file) {
      
      device <- function(..., width, height) {
        grDevices::pdf(file = file,
                       width = width,
                       height = height)
      }
      
      ggsave(
        file,
        plot = molecule_plot(),
        device = device,
        width = input$molecule_plot_width,
        units = "mm",
        height = input$molecule_plot_height,
        limitsize = FALSE
      )
    }
  )
  
  # Download the tumour plot on click
  output$molecule_plot_download_png = downloadHandler(
    filename = "Molecule plot.png",
    content = function(file) {
      png(file, height = input$molecule_plot_height, 
          width = input$molecule_plot_width, unit = "mm", res = 600)
      print(molecule_plot())
      dev.off()
    }
  )
  
  # Display the donor info
  output$table <- DT::renderDataTable({
    DT::datatable(donor_stats, options = list(lengthMenu = c(5, 30, 50), pageLength = 5),  filter = "top")
  })
  
  output$composite_download <- downloadHandler(
    filename <- function() {
      "Composite.jpg"
    },
    
    content <- function(file) {
      
      req(input$slide)
      req(input$FOV)
      
      outfile <- list.files (paste0("/oldvol/shiny/composite/", input$slide, "/CellComposite/"), full.names = T)
      
      fov <- gsub(".*F0|_composite_autocontrast", "", outfile)
      fov <- gsub("^0|\\.jpg", "", fov)
      
      outfile <- outfile[input$FOV == fov]

      file.copy(outfile, file)
    },
    contentType = "image/jpg"
  )
  
}

# Run the app ----
shinyApp(ui, server)
