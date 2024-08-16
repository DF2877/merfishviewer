## Should come after spatial-seurat-prep.R ##

## Something about the data loading is wonky ##
## After loading both files, you will need to relaunch the app ##

library(shiny)
library(shinydashboard)
library(Seurat)
library(shinycssloaders)
library(shinyFiles)
library(bslib)
library(scattermore)
library(ggplot2)
library(Biobase)
library(shinybusy)
library(plotrix)
library(shinyWidgets)
library(dplyr)
library(gtools)

options(bitmapType='cairo')
options(spinner.type=sample(1:6,1),spinner.color.background = "#38023B")

identifier <- list() # Global list to store selected positive gene identifiers
antiID <- list() # Global list to store selected negative gene identifiers
if(exists("vizgen.obj")) { # If data is already loaded in the environment, set the axis range for plotting
  axisrange <- c(min(floor(min(vizgen.obj@meta.data[["x"]], na.rm = TRUE)),
                     floor(min(vizgen.obj@meta.data[["y"]], na.rm = TRUE))),
                 max(ceiling(max(vizgen.obj@meta.data[["x"]], na.rm = TRUE)),
                     ceiling(max(vizgen.obj@meta.data[["y"]], na.rm = TRUE))))
} else {
  axisrange <- c(0,0)
}
brainLabel <- if(exists("vizgen.obj")) {"Brain Loaded"} else {"Load Brain Data"}
sectionLabel <- if(exists("sections")) {"Sections Loaded"} else {"Load Section Data"}
recentBrain <- "" # Stores the most recent filepath used to load brain data
recentSection <- "" # Stores the most recent filepath used to load section data
dataCache <- list() # Cache for precalculated plotting data
lVizCache <- list() # Cache for precalculated subsets, used by DE and violin plots
ade <- reactiveValues(a=FALSE) # Reactive tracker for buttons
markersDL <- NULL # Global var used to export DE to the download handler
interestDL <- NULL # Global var used to export DE to the download handler
minify <- reactiveValues(m=FALSE) # Reactive tracker for the minimize button
subclasses <- if(exists("vizgen.obj")) {c(unique(vizgen.obj@meta.data[["subclass"]]),unique(vizgen.obj@meta.data[["supertype"]]))} else {NA} # If data is already loaded, set the subclasses for selectize
if(file.exists("D:/df/targets_and_families.csv")) {targetIDs <- read.csv("D:/df/targets_and_families.csv", skip=1)}
cr <- colorRamp(c("white", "red")) # Colors!
palette( # More colors!
  c(
    rgb(t(col2rgb("white")), max = 255, alpha = 20),
    c(
      rgb(t(col2rgb("red")), max = 255, alpha = 200),
      rgb(t(col2rgb("cyan")), max = 255, alpha = 200),
      rgb(t(col2rgb("green")), max = 255, alpha = 200),
      rgb(t(col2rgb("violet")), max = 255, alpha = 200),
      rgb(t(col2rgb("yellow")), max = 255, alpha = 200),
      rgb(t(col2rgb("magenta")), max = 255, alpha = 200),
      rgb(t(col2rgb("orange")), max = 255, alpha = 200),
      rgb(t(col2rgb("pink")), max = 255, alpha = 200)
    )
  )
)

spinners <- c("double-bounce", "circle", "bounce", "folding-cube", 
              "rotating-plane", "cube-grid", "fading-circle", "dots", "cube", 
              "flower", "pixel", "hollow-dots", "intersecting-circles", 
              "orbit", "radar", "scaling-squares", "half-circle", 
              "trinity-rings", "fulfilling-square", "circles-to-rhombuses", 
              "semipolar", "self-building-square", "swapping-squares", 
              "fulfilling-bouncing-circle", "fingerprint", "spring", "atom", 
              "looping-rhombuses", "breeding-rhombus")
spinners2 <- sample(c("standard", "hourglass", "circle", "arrows", "dots", "pulse"), 1)

ui <- dashboardPage(
  dashboardHeader(title = "Data Viewer"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu( id = "tabs",
      menuItem("Section View", tabName = "sectionView", icon = icon("magnifying-glass")),
      menuItem("View All", tabName = "overview", icon = icon("grip")),
      menuItem("Differential Expression",tabName = "de",icon = icon("compact-disc"))
    ),
    materialSwitch(inputId = "classView", label = "Class View", status = "info", right = TRUE),
    selectizeInput(
      "subclass",
      "Subclass",
      choices = NULL,
      selected = "",
      multiple = TRUE
    ),
    selectizeInput(
      "symbol",
      "Gene (+)",
      # Try to use Whole brain or section data, if either is available
      tryCatch(vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]],silent = TRUE,error = function(cond){
        if (exists("sections")) {
          sections[[1]]@assays[["Vizgen"]]@meta.data[["gene_symbol"]]
        } else {NA}
      }),
      selected = "Chat",
      multiple = TRUE
    ),
    selectizeInput(
      "symbolMinus",
      "Gene (-)",
      # Try to use Whole brain or section data, if either is available
      tryCatch(vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]],silent = TRUE,error = function(cond){
        if (exists("sections")) {
          sections[[1]]@assays[["Vizgen"]]@meta.data[["gene_symbol"]]
        } else {NA}
      }),
      selected = NULL,
      multiple = TRUE
    ),
    materialSwitch(inputId = "strictMarker", label = "Strict DE Selection?", status = "info", right = TRUE),
    textOutput("test"),
    selectInput(
      "section",
      "Section",
      tryCatch(mixedsort(names(sections)),silent = TRUE,error = function(cond){NA})
    ),
    actionBttn("deBtn","DE for Selection"),
    sliderInput("plotSize","Plot Size (Bigger is slower)", 360,1440, 720, step = 360),
    shinyFilesButton("getBrain", brainLabel,
                     title = "Please select a file (.rds):", multiple = FALSE,
                     buttonType = "default", class = NULL),
    
    shinyFilesButton("getSections", sectionLabel,
                     title = "Please select a file (.rds):", multiple = FALSE,
                     buttonType = "default", class = NULL)
  ),
  ## Body content
  dashboardBody(
    tabItems(
      tabItem(tabName = "sectionView",
              box(
                splitLayout(sliderInput("size","Area", 0, 1, 0.7, step = 0.01),
                sliderInput("target", "Target", 0, 0.7, 0.35, step = 0.01),
                sliderInput("x", "X", axisrange[1],axisrange[2],round(mean(axisrange)), step = 0.05),
                sliderInput("y","Y",axisrange[1],axisrange[2],round(mean(axisrange)), step = 0.05),
                cellArgs = list(style = "padding: 12px")),
              width = "auto"),
              textOutput("mini"),
              uiOutput("main")
      ),
      tabItem(tabName = "overview",
              box(
                withSpinner(
                  uiOutput("plots", width = "100%")
                ),width = "auto"
              )
      ),
      tabItem(tabName = "de",
              splitLayout(downloadButton("downloadM","Export Full Table"),
                          downloadButton("downloadI","Export GPCR/Ion Channel Only")),
              box(
                splitLayout(tableOutput("markers"),
                            tableOutput("interest")),
                width = "auto"))
    ),
    use_busy_spinner(spin = sample(spinners,1), 180, width = 180, position = "top-left")
  )
)

server <- function(input, output, session) {
  
  # output$test <- renderText({input$subclass})
  
  # Initialize the subclass selectize
  # This one runs on the server-side because because it has to search many more choices
  updateSelectizeInput(
    inputId = "subclass",
    choices = tryCatch(c(unique(vizgen.obj@meta.data[["subclass"]]),unique(vizgen.obj@meta.data[["supertype"]])),silent = TRUE,error = function(cond){NA}),
    server = TRUE
  )
  
  # Render main display
  output$main <- renderUI({
    if (!minify$m) {
      fluidRow(
        box(
          uiOutput("sectionBox"),width = 6, height = "1200"
        ),
        box(
          splitLayout(
            actionBttn("vlnBtn","Make Violin Plots"),
            actionBttn("sbcBtn","ID Subclasses"),
            textOutput("vlnID"),
            cellArgs = list(style = "padding: 12px"),
            actionBttn("mini","Minimize")
          ),
          h4("Target Selection"),
          plotOutput("vlnTarget"),
          h4("Non-target Area"),
          plotOutput("vlnArea"),
          width=6, height = "1200")
      )
    } else {
      fluidRow(
        box(
          uiOutput("sectionBox"),width = 11, height = "1560"
        ),
        box(
          actionBttn("mini","<"),
          width=1, height = "auto")
      )
    }
  })
  
  # Reactive observer for minimize
  observeEvent(input$mini, {
    minify$m <<- !minify$m
  })
  
  # Render UI container for main plot
  output$sectionBox <- renderUI({
    withSpinner(
      plotOutput("sectionPlot",width = input$plotSize, height = input$plotSize), hide.ui = FALSE
    )
  })
  
  # File upload handler (Whole brain)
  volumes = getVolumes()
  bindEvent(observe(input$getBrain), ignoreNULL = TRUE, {
    shinyFileChoose(input, "getBrain", roots = volumes, session = session)
    if(!is.null(input$getBrain)){
      show_spinner()
      block(
        "getBrain",
        type = spinners2
      )
      file_selected<-parseFilePaths(volumes, input$getBrain)
      try({
        # Since file upload is slow, prevent reading the same file twice in a row
        if(file_selected$datapath != recentBrain) {
          temp <- readRDS(file_selected$datapath)
          if(class(temp) == "Seurat") {
            # Vars need to be set globally
            recentBrain <<- file_selected$datapath
            vizgen.obj <<- temp
            axisrange <<- c(min(floor(min(vizgen.obj@meta.data[["x"]], na.rm = TRUE)),
                                floor(min(vizgen.obj@meta.data[["y"]], na.rm = TRUE))),
                            max(ceiling(max(vizgen.obj@meta.data[["x"]], na.rm = TRUE)),
                                ceiling(max(vizgen.obj@meta.data[["y"]], na.rm = TRUE))))
            # Update relevant UI elements
            updateSliderInput(inputId = "x", min=axisrange[1],max=axisrange[2],value=round(mean(axisrange)), step = 0.05)
            updateSliderInput(inputId = "y",min=axisrange[1],max=axisrange[2],value=round(mean(axisrange)), step = 0.05)
            updateSelectizeInput(
              inputId = "subclass",
              choices = tryCatch(c(unique(vizgen.obj@meta.data[["subclass"]]),unique(vizgen.obj@meta.data[["supertype"]])),silent = TRUE,error = function(cond){NA}),
              server = TRUE
            )
            updateSelectizeInput(
              inputId = "symbol",
              choices = tryCatch(vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]],silent = TRUE,error = function(cond){NA})
            )
            updateSelectizeInput(
              inputId = "symbolMinus",
              choices = tryCatch(vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]],silent = TRUE,error = function(cond){NA})
            )
            updateSelectInput(
              inputId = "section",
              choices=tryCatch(mixedsort(na.omit(unique(vizgen.obj@meta.data[["z"]]))),silent = TRUE,error = function(cond){NA})
            )
            brainLabel <<- "Brain Loaded"
            try(identifier <<- vizgen.obj@assays[["Vizgen"]]@meta.data[which(vizgen.obj@assays[["Vizgen"]]@meta.data[,2]==input$symbol),1])
            # message(identifier)
            mkPlt(input$section, identifier, antiID, input$subclass)
          } else {notify_failure("Invalid file", timeout = 30000, position = "center-bottom")}
        }
      }, silent = TRUE)
      hide_spinner()
      unblock("getBrain")
    }
  })
  
  # File upload handler (Sections)
  bindEvent(observe(input$getSections), {
    shinyFileChoose(input, "getSections", roots = volumes, session = session)
    if(!is.null(input$getSections)){
      show_spinner()
      block(
        "getSections",
        type = spinners2
      )
      file_selected2<-parseFilePaths(volumes, input$getSections)
      try({
        # Since file upload is slow, prevent reading the same file twice in a row
        if(file_selected2$datapath != recentSection) {
          temp <- readRDS(file_selected2$datapath)
          if(class(temp) == "list") {
            recentSection <<- file_selected2$datapath
            sections <<- temp
            updateSelectInput(
              inputId = "section",
              choices=tryCatch(mixedsort(names(sections)),silent = TRUE,error = function(cond){NA})
            )
            updateSelectizeInput(
              inputId = "symbol",
              choices = tryCatch(sections[[1]]@assays[["Vizgen"]]@meta.data[["gene_symbol"]],silent = TRUE,error = function(cond){NA})
            )
            updateSelectizeInput(
              inputId = "symbolMinus",
              choices = tryCatch(vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]],silent = TRUE,error = function(cond){NA})
            )
            sectionLabel <<- "Sections Loaded"
            try(identifier <<- sections[[1]]@assays[["Vizgen"]]@meta.data[which(sections[[1]]@assays[["Vizgen"]]@meta.data[,2]==input$symbol),1], silent = TRUE)
            # message(identifier)
            mkPlt(input$section, identifier, antiID, input$subclass)
          } else {notify_failure("Invalid file", timeout = 30000, position = "center-bottom")}
        }
      },silent = TRUE)
      hide_spinner()
      unblock("getSections")
    }
  })
  
  # Retrieve spatial expression data
  sectionData <- function(sec, id, antiID, sbc) {
    # message("Calculating")
    # If the desired data is already cached, retrieve and return it
    if (paste0(paste0(as.character(sec),unlist(id),collapse = "+"),"-",paste0(unlist(antiID),collapse = "-"),paste0(unlist(sbc),collapse = "+")) %in% names(dataCache)) {
      dataCache[[paste0(paste0(as.character(sec),unlist(id),collapse = "+"),"-",paste0(unlist(antiID),collapse = "-"),paste0(unlist(sbc),collapse = "+"))]]
    } else { # Otherwise, calculate and save new
      df <- FetchData(sections[[as.character(sec)]],c("x","y","subclass","supertype",unlist(id),unlist(antiID)),clean = "all")
      df <- cbind(aggr = 1, df)
      # Collapse expression values from each gene down to one number
      # I've done this by multiplying positive markers and dividing negative markers
      # This method strongly selects for cells with complete matching to the criteria
      # Partial matches are not able to get high values
      for (i in id) {
        df[1] <- df[1] * (df[i]+0.0000001) # Add small amount to avoid div by zero
      }
      for (j in antiID) {
        df[1] <- df[1] / (df[j]+0.0000001) # Must add the same amount as above
      }
      # Identify subclasses, and store alpha values
      subclass <- if (length(sbc)>0) {sapply(sapply(df[4], function(l){l%in%sbc})|sapply(df[5], function(l){l%in%sbc}), function(k){if(k){200}else{20}})} else {200}
      # Reconstruct dataframe without individual gene expression, only aggregate
      DF <- data.frame(x=as.numeric(df[[2]]),y=as.numeric(df[[3]]),feat=as.numeric(df[[1]]),subclass=subclass)
      # Create a copy with rounded values
      DF2 <- data.frame(x=round(as.numeric(df[[2]]),2),y=round(as.numeric(df[[3]]),2),feat=round(as.numeric(df[[1]]),1),subclass=subclass)
      # Remove cells with too much overlap
      # This only removes cells if ALL rounded columns are too similar
      # Speeds up plotting with minor information loss
      DF2 <- DF[!duplicated(DF2),]
      colnames(DF2) <- c("x","y",paste0(paste0(unlist(input$symbol),collapse = "+"),"-",paste0(unlist(input$symbolMinus),collapse = "-")),"subclass")
      # Clear old cached values if cache size is over an arbitrary limit
      if (length(dataCache) >= 512) {dataCache <- dataCache[-1:-255]}
      # Store processed data in the cache
      # The formula for calculating cache names must be unique for any change in parameters
      # The formula must be the same when storing and retrieving
      dataCache[[paste0(paste0(as.character(sec),unlist(id),collapse = "+"),"-",paste0(unlist(antiID),collapse = "-"),paste0(unlist(sbc),collapse = "+"))]] <<- DF2
      DF2
    }
  }
  
  # Retrieve cell class data and generate colors
  classData <- function(sec, sbc) {
    # message("Calculating")
    # If the desired data is already cached, retrieve and return it
    if (paste0(paste0(as.character(sec),collapse = "+"),"-",paste0(unlist(sbc),collapse = "+")) %in% names(dataCache)) {
      dataCache[[paste0(paste0(as.character(sec),collapse = "+"),"-",paste0(unlist(sbc),collapse = "+"))]]
    } else { # Otherwise, calculate and save new
      df <- FetchData(sections[[as.character(sec)]],c("x","y","subclass","supertype"),clean = "all")
      # Identify subclasses, and store alpha values
      if (length(sbc)>0) {
        subclass <- sapply(1:nrow(df), function(m) {if (df[m,4] %in% sbc) {
          match(df[m,4],sbc)+1
          } else if (df[m,3] %in% sbc) {
            match(df[m,3],sbc)+1
            } else {1}})
        } else {subclass <- 1}
      # Reconstruct dataframe without individual gene expression, only aggregate
      DF <- data.frame(x=as.numeric(df[[1]]),y=as.numeric(df[[2]]),subclass=subclass)
      # Create a copy with rounded values
      DF2 <- data.frame(x=round(as.numeric(df[[1]]),2),y=round(as.numeric(df[[2]]),2))
      # Remove cells with too much overlap
      # This only removes cells if ALL rounded columns are too similar
      # Speeds up plotting with minor information loss
      DF2 <- DF[!duplicated(DF2),]
      colnames(DF2) <- c("x","y","subclass")
      # Clear old cached values if cache size is over an arbitrary limit
      if (length(dataCache) >= 512) {dataCache <- dataCache[-1:-255]}
      # Store processed data in the cache
      # The formula for calculating cache names must be unique for any change in parameters
      # The formula must be the same when storing and retrieving
      dataCache[[paste0(paste0(as.character(sec),collapse = "+"),"-",paste0(unlist(sbc),collapse = "+"))]] <<- DF2
      DF2
    }
  }
  
  # Render main plot
  mkPlt <- function(sectionn, id, antiID, sbc) {
    # start.time <- Sys.time()
    tryCatch({
      if (input$classView) {
        DF2 <- classData(sectionn, sbc)
        #
        output$sectionPlot <- renderPlot({
          par(bg="black")
          scattermoreplot(DF2$x,DF2$y,pch=".", col = unlist(DF2[3]),
                          size = c(input$plotSize,input$plotSize),
                          cex=sqrt(input$plotSize)/13, ylim = rev(axisrange), 
                          xlim = axisrange)
          draw.circle(input$x, input$y, input$size, border = "steelblue", lwd = 2)
          draw.circle(input$x, input$y, input$target, border = "forestgreen", lwd = 2)
          legend("topleft", legend=c("",unlist(input$subclass)), pch=16, col=palette(), text.col = "white")
        })
        # palette("default")
      } else {
        if (!is.null(sections[[as.character(sectionn)]])) {
          DF2 <- sectionData(sectionn,id,antiID, sbc)
          output$sectionPlot <- renderPlot({
            par(bg="black")
            scattermoreplot(DF2$x,DF2$y,pch=".", 
                            col = rgb(cr(unlist(DF2[3]) / max(DF2[3])),max=255, alpha=unlist(DF2[4])),
                            size = c(input$plotSize,input$plotSize),
                            cex=sqrt(input$plotSize)/13, ylim = rev(axisrange), 
                            xlim = axisrange)
            draw.circle(input$x, input$y, input$size, border = "steelblue", lwd = 2)
            draw.circle(input$x, input$y, input$target, border = "forestgreen", lwd = 2)
          })
          # end.time <- Sys.time()
          # time.taken <- end.time - start.time
          # message(time.taken)
        } else {
          # message(sectionn,id)
          output$sectionPlot <- renderPlot(plot.new())}
      }
    },error = function(cond) {
      # message("wheee")
      message(conditionMessage(cond))
      output$sectionPlot <- renderPlot(plot.new())})
  }
  
  # Reactive observer for positive gene selector
  bindEvent(observe(input$symbol), {
    # message("updatesy")
    # Do not recalculate violin plots (They will still re-render)
    ade$v <<- FALSE
    # Convert the list of gene symbols to Ensmble identifiers
    # The ABC merfish data uses Ensmble ids. If your data uses symbols, you could change this to simply return the input without any calculation
    try(
      if(length(input$symbol>0)){
        identifier <<- list()
        for (i in 1:length(input$symbol)) {
          identifier[[i]] <<- vizgen.obj@assays[["Vizgen"]]@meta.data[which(vizgen.obj@assays[["Vizgen"]]@meta.data[,2]==input$symbol[[i]]),1]
        }
      } else {identifier <<- list()}
    )
    mkPlt(input$section, identifier, antiID, input$subclass)
  })
  
  # Reactive observer for negative gene selector
  bindEvent(observe(input$symbolMinus), {
    # Do not recalculate violin plots (They will still re-render)
    ade$v <<- FALSE
    # Convert the list of gene symbols to Ensmble identifiers
    # The ABC merfish data uses Ensmble ids. If your data uses symbols, you could change this to simply return the input without any calculation
    try(
      if(length(input$symbolMinus)>0){
        antiID <<- list()
        for (j in 1:length(input$symbolMinus)) {
          antiID[[j]] <<- vizgen.obj@assays[["Vizgen"]]@meta.data[which(vizgen.obj@assays[["Vizgen"]]@meta.data[,2]==input$symbolMinus[[j]]),1]
        }
      } else {antiID <<- list()}
    )
    # Trigger plot generation. It is only important during initialization. 
    # Once created, the plot will re-render automatically. It needs to be in a reactive object to initialize, but which one is arbitrary
    mkPlt(input$section, identifier, antiID, input$subclass)
  })
  
  # Reactive observer for area slider
  bindEvent(observe(input$size), {
    # Force violin plots to unrender by overwriting them with blank plots
    # Otherwise they will automatically recalculate, which is slow and unnecessary
    output$vlnTarget <- renderPlot(plot.new())
    output$vlnArea <- renderPlot(plot.new())
    ade$v <<- FALSE
    # Update the max value of the target slider to match
    updateSliderInput(inputId = "target",max = input$size, step = 0.01)
    # mkPlt(input$section, identifier, antiID, input$subclass)
  })
  
  # Reactive observer for x slider
  bindEvent(observe(input$x), {
    # Force violin plots to unrender by overwriting them with blank plots
    # Otherwise they will automatically recalculate, which is slow and unnecessary
    output$vlnTarget <- renderPlot(plot.new())
    output$vlnArea <- renderPlot(plot.new())
    ade$v <<- FALSE
  #   mkPlt(input$section, identifier, antiID, input$subclass)
  })
  
  # Reactive observer for y slider
  bindEvent(observe(input$y), {
    # Force violin plots to unrender by overwriting them with blank plots
    # Otherwise they will automatically recalculate, which is slow and unnecessary
    output$vlnTarget <- renderPlot(plot.new())
    output$vlnArea <- renderPlot(plot.new())
    ade$v <<- FALSE
  #   mkPlt(input$section, identifier, antiID, input$subclass)
  })
  
  # Reactive observer for section selector
  bindEvent(observe(input$section), {
    # Force violin plots to unrender by overwriting them with blank plots
    # Otherwise they will automatically recalculate, which is slow and unnecessary
    output$vlnTarget <- renderPlot(plot.new())
    output$vlnArea <- renderPlot(plot.new())
    ade$v <<- FALSE
  #   mkPlt(input$section, identifier, antiID, input$subclass)
  })
  
  # Render UI container for multi-section plots
  output$plots <- renderUI({
    # Insert the right number of plot output objects into the web page
    plot_output_list <- lapply(rev(mixedsort(names(sections))), function(i) {
      plotname <- paste("plot", i, sep="")
      withSpinner(
        plotOutput(plotname, height = 240)
      )
    })
    
    withProgress({
      show_spinner()
      maxCol <- 0
      # Calculate plot data, and find maximum expression value
      lapply(rev(mixedsort(names(sections))), function(i) {
        incProgress(0.5/length(sections))
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          plotname <- paste("plot", my_i, sep="")
          if (length(unlist(identifier))>0) {
            sbc <- input$subclass
            # If the desired data is already cached, retrieve and return it
            if (paste0("s",paste0(as.character(my_i),unlist(identifier),
                                  collapse = "+"),"-",paste0(unlist(antiID),
                                                             collapse = "-"),
                       paste0(unlist(sbc),
                              collapse = "+")) %in% names(dataCache)) {
              DF2 <- dataCache[[paste0("s",paste0(as.character(my_i),
                                                  unlist(identifier),
                                                  collapse = "+"),"-",
                                       paste0(unlist(antiID),collapse = "-"),
                                       paste0(unlist(sbc),collapse = "+"))]]
            } else {
              
              if (input$classView) {
                df <- FetchData(sections[[as.character(my_i)]],c("x","y","subclass","supertype"),clean = "all")
                df <- cbind(aggr = 1, df)
                if (length(sbc)>0) {subclass <- sapply(1:nrow(df), function(m) {
                  if (df[m,5] %in% sbc) {match(df[m,5],sbc)+1} else if (df[m,4] %in% sbc) {
                    match(df[m,4],sbc)+1
                    } else {1}
                  })} else {subclass <- 1}
              } else {
                df <- FetchData(sections[[as.character(my_i)]],c("x","y","subclass","supertype",unlist(identifier),unlist(antiID)),clean = "all")
                df <- cbind(aggr = 1, df)
                # Collapse expression values from each gene down to one number
                # I've done this by multiplying positive markers and dividing negative markers
                # This method strongly selects for cells with complete matching to the criteria
                # Partial matches are not able to get high values
                for (i in identifier) {
                  df[1] <- df[1] * (df[i]+0.0000001)
                }
                for (j in antiID) {
                  df[1] <- df[1] / (df[j]+0.0000001)
                }
                subclass <- if (length(sbc)>0) {sapply(sapply(df[4], function(l){l%in%sbc})|sapply(df[5], function(l){l%in%sbc}), function(k){if(k){200}else{10}})} else {200}
              }
              # Reconstruct dataframe without individual gene expression, only aggregate
              DF <- data.frame(x=as.numeric(df[[2]]),y=as.numeric(df[[3]]),feat=as.numeric(df[[1]]),subclass=subclass)
              # Create a copy with rounded values
              DF2 <- data.frame(x=round(as.numeric(df[[2]]),2),y=round(as.numeric(df[[3]]),2),feat=round(as.numeric(df[[1]]),1))
              # Remove cells with too much overlap
              # This only removes cells if ALL rounded columns are too similar
              # Speeds up plotting with minor information loss
              DF2 <- DF[!duplicated(DF2),]
              colnames(DF2) <- c("x","y",paste0(paste0(unlist(input$symbol),collapse = "+"),"-",paste0(unlist(input$symbolMinus),collapse = "-")),"subclass")
              # Clear old cached values if cache size is over an arbitrary limit
              if (length(dataCache) >= 512) {dataCache <- dataCache[-1:-255]}
              # Store processed data in the cache
              # The formula for calculating cache names must be unique for any change in parameters
              # The formula must be the same when storing and retrieving
              dataCache[[paste0("s",paste0(as.character(my_i),
                                           unlist(identifier),collapse = "+"),
                                "-",paste0(unlist(antiID),collapse = "-"),
                                paste0(unlist(sbc),collapse = "+"))]] <<- DF2
            }
            # Determine the maximum expression value
            if (max(DF2[3]) > maxCol) {
              maxCol <<- max(DF2[3])
            }
          }
        })
      })
      # Render plots
      lapply(rev(mixedsort(names(sections))), function(i) {
        incProgress(0.5/length(sections))
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
          my_i <- i
          plotname <- paste("plot", my_i, sep="")
          # Call renderPlot for each one. Plots are only actually generated when they
          # are visible on the web page.
          output[[plotname]] <- renderPlot({
            if (length(unlist(identifier))>0) {
              # Retrieve plot data (because it was  calculated in a different local env)
              DF2 <- dataCache[[paste0("s",paste0(as.character(my_i),
                                                  unlist(identifier),
                                                  collapse = "+"),"-",
                                       paste0(unlist(antiID),collapse = "-"),
                                       paste0(unlist(input$subclass),
                                              collapse = "+"))]]
              par(bg="black")
              if (input$classView) {
                scattermoreplot(DF2$x,DF2$y,pch=".", 
                                main = paste(paste0(unlist(input$subclass),collapse = "+"),"\n",as.character(my_i)), 
                                col.main = "white", col = unlist(DF2[4]), 
                                size = c(240,240), ylim = rev(axisrange), 
                                xlim = axisrange)
                # legend("topleft", legend=c("",unlist(input$subclass)), pch=16, col=palette(), text.col = "white")
              } else {
                scattermoreplot(DF2$x,DF2$y,pch=".", 
                                main = paste(paste0(unlist(input$symbol),collapse = "+"),"\n",as.character(my_i)), 
                                col.main = "white", 
                                col = rgb(cr(unlist(DF2[3]) / max(DF2[3])),max=255, alpha=unlist(DF2[4])), 
                                size = c(240,240), ylim = rev(axisrange), 
                                xlim = axisrange)
              }
            }
          }, height = 240, width = 240)
        })
      })
      hide_spinner()
    }, message = "Plotting...", value = 0)
    do.call(flowLayout, plot_output_list)
  })
  
  # Reactive observer for DE button
  # This doesn't do anything
  # It was supposed to switch the user's view to the DE tab
  # Now it just triggers the tracker
  observeEvent(input$deBtn, {
    show_spinner()
    updateTabItems(session, "tabs","de")
    ade$a <<- TRUE
    hide_spinner()
  },ignoreInit = TRUE)
  
  # Reactive observer for DE button tracker
  observeEvent(ade$a, {
    if (ade$a) {
      ade$a <<- FALSE
      doDE()
    }
  })
  
  # Calculate differentially expressed genes
  doDE <- function() {
    # message("function")
    withProgress({
      localVizgen <- getLocalSubset()
      incProgress(0.25)
      #Find marker genes for each cluster
      #Default is Wilcoxon rank sum test
      markers <- FindMarkers(localVizgen, ident.1 = TRUE, group.by = "io")
      incProgress(0.25)
      markers <- markers[markers$p_val_adj < 0.05,]
      incProgress(0.25)
      if (exists("targetIDs")) {
        targets <- targetIDs$Type[match(rownames(markers),targetIDs$MGI.symbol)]
        # mgiSymbols <- localVizgen@assays[["Vizgen"]]@meta.data$gene_symbol[match(rownames(markers),localVizgen@assays[["Vizgen"]]@meta.data$gene_identifier)]
        markers <- cbind(markers, targets)
        # rownames(markers) <- make.names(mgiSymbols, unique = TRUE)
      }
      interest <- markers[markers$targets %in% c("gpcr","vgic","lgic","other_ic"),]
      # rownames(markers) <- vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_symbol"]][match(rownames(markers),vizgen.obj@assays[["Vizgen"]]@meta.data[["gene_identifier"]])]
      output$markers <- renderTable(markers, 
                                    striped = TRUE,
                                    hover = TRUE,
                                    rownames = TRUE)
      output$interest <- renderTable(interest,
                                     striped = TRUE,
                                     hover = TRUE,
                                     rownames = TRUE)
      # Send results to global var so that the download handler can access it
      markersDL <<- markers
      interestDL <<- interest
    }, message = "Analysing...", value = 0.25)
  }
  
  # Subset target cells
  # It's pretty slow but it makes subsequent operations much faster
  getLocalSubset <- function() {
    # Generate cache name addon if strict selection is being used
    symbolString <- if (input$strictMarker) {
      paste0(input$strictMarker,paste0(unlist(input$symbol),unlist(input$symbolMinus),unlist(input$subclass),collapse = "-"))
    } else {""}
    isolate(
      # If subset has already been done, retrieve and return it
      if (paste(sep="-",input$target,input$size,input$x,input$y,input$section,symbolString,collapse = "-") %in% names(lVizCache)) {
        lVizCache[[paste(sep="-",input$target,input$size,input$x,input$y,input$section,symbolString,collapse = "-")]]
      } else {
        withProgress({
          # Find cells near a point
          distances <- list()
          dat <- list()
          for (z in names(sections)) {
            distances[[z]] <- matchpt(as.matrix(sections[[z]][[c("x","y")]]),
                                      as.matrix(data.frame(input$x,input$y)))
            distances[[z]]$distance <- sqrt(((as.numeric(z)-as.numeric(input$section))**2)+(distances[[z]]$distance**2))
            dat[[z]] <- distances[[z]][rownames(distances[[z]][distances[[z]]$distance<input$size,]),]
            # Mark whether cells are within the inner or outer area
            dat[[z]]$io <- dat[[z]]$distance<input$target
            incProgress(0.4/length(names(sections)))
          }
          # Combine individual lists of cells per section into one list
          datdat <- do.call(rbind,unname(dat))
          # Additional filtering if strict selection is enabled
          if (input$strictMarker) {
            DF2 <- data.frame()
            datLen <- vector()
            # Find lengths of each list of cells
            # The lists would be empty if that section did not intersect the target area
            for(i in 1:length(dat)) {
              datLen[names(sections)[i]] <- nrow(dat[[i]])
            }
            # Get expression data for each relevant section
            for (aSection in names(which(datLen>0))) {
              temp <- sectionData(aSection,identifier,antiID, input$subclass)
              # Combine them all into one big dataframe as we go, no need for them to be separate
              DF2 <- rbind(DF2, temp)
            }
            # DF2 <- sectionData(input$section,identifier,antiID)
            feat <- paste0(paste0(unlist(input$symbol),collapse = "+"),"-",paste0(unlist(input$symbolMinus),collapse = "-"))
            # Check each cell. If it does not express the desired combination of genes, 
            # or is not part of the desired subclass, mark it as part of the outer group, regardless of spatial location
            for (aRow in rownames(datdat)) {
              tryCatch({
                if (DF2[aRow, feat] < 0.000001 || DF2[aRow, "subclass"] < 200) {
                  datdat[aRow,"io"] <- FALSE
                }
              }, error = function(cond) {message(aRow)})
              incProgress(0.4/length(rownames(datdat)))
            }
          }
          datCells <- rownames(datdat)
          # Create a subset of cells within the target area
          # This is the slow part
          localVizgen <- suppressWarnings(subset(vizgen.obj, subset = cell_label %in% datCells))
          localVizgen <- AddMetaData(localVizgen, datdat)
          rownames(localVizgen) <- localVizgen@assays[["Vizgen"]]@meta.data[["gene_symbol"]]
          incProgress(0.2)
        }, message = "Subsetting...", value = 0)
        # Clear old cached values if cache size is over an arbitrary limit
        if (length(lVizCache) >= 128) {lVizCache <- lVizCache[-1:-63]}
        # Store processed data in the cache
        # The formula for calculating cache names must be unique for any change in parameters
        # The formula must be the same when storing and retrieving
        lVizCache[[paste(sep="-",input$target,input$size,input$x,input$y,input$section,symbolString,collapse = "-")]] <<- localVizgen
        localVizgen
      }
    )
  }
  
  # Reactive observer for the violin plot button
  # This doesn't do anything
  # It was supposed to switch the user's view to the DE tab
  # Now it just triggers the tracker
  observeEvent(input$vlnBtn, {
    show_spinner()
    ade$v <<- TRUE
    hide_spinner()
  },ignoreInit = TRUE)
  
  # Reactive observer for the violin plot button tracker
  observeEvent(ade$v, {
    # Render the violin plots
    if (ade$v) {
      ade$v <<- FALSE
      output$vlnID <- renderText({
        isolate(paste("Area:",input$size,"Target:",input$target,"X:",input$x,"Y:",input$y,"Section:",input$section))
      })
      output$vlnTarget <- renderPlot({
        VlnPlot(subset(getLocalSubset(), subset = io == TRUE), c(unlist(input$symbol),unlist(input$symbolMinus)), ncol = 5, same.y.lims = TRUE)
      })
      output$vlnArea <- renderPlot({
        VlnPlot(subset(getLocalSubset(), subset = io == FALSE), c(unlist(input$symbol),unlist(input$symbolMinus)), ncol = 5, same.y.lims = TRUE)
      })
    }
  })
  
  # Download handler
  output$downloadM <- downloadHandler(
    filename = function() {paste0(paste(sep="-",input$target,input$size,input$x,
                                        input$y,input$section,
                                        paste(sep="-",input$strictMarker,
                                              paste0(unlist(input$symbol),
                                                     unlist(input$symbolMinus),
                                                     unlist(input$subclass),
                                                     collapse = "-")),
                                        collapse = "-"),".csv")},
    content = function(fname) {
      write.csv(markersDL,fname)
    }
  )
  
  # Download handler
  output$downloadI <- downloadHandler(
    filename = function() {paste0(paste(sep="-",input$target,input$size,input$x,
                                        input$y,input$section,
                                        paste(sep="-",input$strictMarker,
                                              paste0(unlist(input$symbol),
                                                     unlist(input$symbolMinus),
                                                     unlist(input$subclass),
                                                     collapse = "-")),
                                        collapse = "-"),"-gpcr-ic.csv")},
    content = function(fname) {
      write.csv(interestDL,fname)
    }
  )
  
  observeEvent(input$sbcBtn, {
    show_spinner()
    ade$s <<- TRUE
    hide_spinner()
  },ignoreInit = TRUE)
  
  observeEvent(ade$s, {
    if (ade$s) {
      ade$s <<- FALSE
      output$vlnID <- renderText({
        isolate(paste("Area:",input$size,"Target:","\n",input$target,"X:",input$x,"Y:",input$y,"\n","Section:",input$section))
      })
      localVizgen <- getLocalSubset()
      output$vlnTarget <- renderPlot({
        iLViz <- subset(localVizgen, subset = io == TRUE)
        sbcCounts <- slice_max(data.frame(table(iLViz@meta.data$subclass)), order_by = table(iLViz@meta.data$subclass), n = 15)
        sptpCounts <- slice_max(data.frame(table(iLViz@meta.data$supertype)), order_by = table(iLViz@meta.data$supertype), n = 15)
        iTopSbc <- slice_max(rbind(sbcCounts,sptpCounts), order_by = rbind(sbcCounts,sptpCounts)$Freq, n = 15)
        rownames(iTopSbc) <- iTopSbc$Var1
        par("mar" = c(5.1, max(nchar(unfactor(iTopSbc$Var1)))/1.5, 4.1, 2.1))
        barplot(t(iTopSbc[2]),horiz = T, las = 1)
        # par("mar" = c(5.1, 4.1, 4.1, 2.1))
      })
      output$vlnArea <- renderPlot({
        oLViz <- subset(localVizgen, subset = io == FALSE)
        sbcCounts <- slice_max(data.frame(table(oLViz@meta.data$subclass)), order_by = table(oLViz@meta.data$subclass), n = 15)
        sptpCounts <- slice_max(data.frame(table(oLViz@meta.data$supertype)), order_by = table(oLViz@meta.data$supertype), n = 15)
        oTopSbc <- slice_max(rbind(sbcCounts,sptpCounts), order_by = rbind(sbcCounts,sptpCounts)$Freq, n = 15)
        rownames(oTopSbc) <- oTopSbc$Var1
        par("mar" = c(5.1, max(nchar(unfactor(oTopSbc$Var1)))/1.5, 4.1, 2.1))
        barplot(t(oTopSbc[2]),horiz = T, las = 1)
        # par("mar" = c(5.1, 4.1, 4.1, 2.1))
      })
    }
  })
  
  # Clean up memory when app is closed
  # Not technically necessary, since RStudio does this automatically
  # But this app can use a lot of memory, and RStudio only does gc when it gets full
  # No point in continuing to claim the memory once we're done
  # It doesn't remove anything, just frees up memory
  session$onSessionEnded(function() {
    gc()
  })
}

shinyApp(ui, server)
