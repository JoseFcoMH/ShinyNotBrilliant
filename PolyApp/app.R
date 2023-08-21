#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

lapply(c('shiny','tidyverse', 'ggthemes', 'pracma', 'ggprism', 'data.table', 'polyprof', 'readxl', 'pzfx', 'openxlsx'), require, character.only = TRUE)

norm_methods <- c('None', 'AUC', 'AUC_Rib', '80S')

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Polysome profiles"),
    
    # Upload files
    fileInput('profiles', 'Upload .csv files straight from the Biocomp fractionator here', multiple = TRUE, accept = c('.csv')),
    tableOutput('files'),
    
    # Layout
    sidebarLayout(
      position = 'left',
      sidebarPanel(
        h4('Profile visualization'),
                   
        # Sliders for peak selection
        sliderInput('PeakStart', label = 'Aligning peak search - start', value = 35, min = 6, max = 66),
        sliderInput('PeakEnd', label = 'Aligning peak search - end', value = 46, min = 6, max = 66),
        
        # Data processing options
        checkboxInput('zbaseline', 'Set all baselines to 0', value = TRUE),
        checkboxInput('smoothen', 'Smooth graph', value = TRUE),
        selectInput('norm', 'AUC normalization method', norm_methods),
        numericInput('maxAbs', 'Maximum absorbance value', value = 2),
        
        actionButton('calcgraph', 'Graph'),
                   
        h4('Polysome/Monosome ratios'), 
                   
        # Sliders for monosome and polysome areas
        sliderInput('MonosomeStart', label = 'Monosome area start', value = 26, min = 6, max = 66),
        sliderInput('MonosomeEnd', label = 'Monosome area end', value = 37, min = 6, max = 66),
        sliderInput('PolysomeStart', label = 'Polysome area start', value = 37, min = 6, max = 66),
        sliderInput('PolysomeEnd', label = 'Polysome area end', value = 60, min = 6, max = 70),
        
        actionButton('calcpm', 'Calculate'),
                   
      ),
      
      # Main panel with stuff
      mainPanel(
        # Plots
        plotOutput("rawPlot"),
        plotOutput("ProfPlot"),
        # Downloads
        downloadButton('postPlotImg', 'Download High-Res image'),
        downloadButton('postPlotDtf', 'Download processed data (default format)'),
        downloadButton('postPlotExcel', 'Download processed data (ready to GraphPad Prism copy-paste, .xlsx) (Recommended)'),
        downloadButton('postPlotPrism', 'Download processed data (GraphPad Prism format, .pzfx) (Not recommended for plots with >100 data points)')
      )
    )
)

# Server serving
server <- function(input, output) {
  
    output$files <- renderTable(input$profiles, rownames = TRUE)
    
    profils1 <- eventReactive(input$profiles, {
      profs <- list()
      n = 1
      for(nr in 1:nrow(input$profiles)){
        f <- input$profiles[nr, 'datapath']
        f <- read_csv(f, skip = 47)
        f$Sample_ID = n
        profs <- c(profs, list(f))
        n = n+1
      }
      return(profs)
    })
    
    profils2 <- eventReactive(input$calcgraph, {
      temp <- lapply(profils1(), align, ref = profils1()[[1]], by_peaks = TRUE, npeaks = 1, minPeakPos = input$PeakStart, maxPeakPos = input$PeakEnd)
      temp <- lapply(temp, normalize, to = input$norm, max_abs = input$maxAbs, max_jump = 0.5, pos_start = 10, pos_end = 70, smoothen = input$smoothen, zero_baseline = input$zbaseline)
      return(temp)
    })
    
    output$rawPlot <- renderPlot({
      ggplot(data = rbindlist(profils1()), aes(`Position(mm)`, Absorbance, color = as.factor(Sample_ID))) +
        geom_line() +
        labs(title = 'Raw plot', color = 'ID') +
        theme_bw()
    })
    
    plotInput <- function(){
      ggplot(data = rbindlist(profils2()), aes(`Position(mm)`, Absorbance, color = as.factor(Sample_ID))) +
        geom_line() +
        labs(title = 'Processed plot', color = 'ID') +
        theme_bw()
    }
    
    output$ProfPlot <- renderPlot({
      plotInput()
    })
    
    output$postPlotImg <- downloadHandler(filename = 'ProcessedProf.tiff', 
                                          content = function(file) {
                                            ggsave(file, plot = plotInput(), width = 6, height = 4, dpi = 600, device = 'tiff')
                                            }
                                          )
    
    output$postPlotDtf <- downloadHandler(filename = 'ProcessedProf.csv', 
                                          content = function(file) {
                                            write.csv(profils2(), file)
                                            }
                                          )
    
    output$postPlotPrism <- downloadHandler(filename = 'ProcessedProfPrism.pzfx', 
                                            content = function(file) {
                                              write_pzfx(PrismExport2(profils2(), mode = 'percentage'),
                                                        file, x_col = 1)
                                              }
                                            )
    
    output$postPlotExcel <- downloadHandler(filename = 'ProcessedProfPrism.xlsx', 
                                            content = function(file) {
                                              write.xlsx(PrismExport2(profils2(), mode = 'percentage'),
                                                         file)
                                            }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
