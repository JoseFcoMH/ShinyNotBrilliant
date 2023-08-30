#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#


lapply(c('shiny', 'DT', 'tidyverse', 'ggthemes', 'pracma', 'ggprism', 'data.table', 'polyprof', 'readxl', 'pzfx', 'openxlsx'), require, character.only = TRUE)

norm_methods <- c('None', 'AUC', 'AUC_Rib', '80S')

IDer <- function(da, iden=c()){
          iden = unlist(str_split(iden, ','))
          if(length(iden) == nrow(da)){
            mer <- da %>% mutate(Sample_ID = iden)
          } else {
            mer <- da %>% mutate(Sample_ID = as.factor(row_number()))
          }
          mer
        }

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Polysome profiles"),
    
    # Upload files
    fileInput('profiles', 'Upload .csv files straight from the Biocomp fractionator here', multiple = TRUE, accept = c('.csv')),
    textInput('IDs', 'Enter your data labels separated by a comma (,)', placeholder = '1,abc,s 1,'),
    DTOutput('files'),

    # Layout
    sidebarLayout(
      position = 'left',
      sidebarPanel(
        h4('Profile visualization'),
                   
        # Sliders for peak selection
        checkboxInput('cutoffsP', 'Show cutoffs for peak selection in Raw plot', value = FALSE),
        sliderInput('PeakStart', label = 'Aligning peak search - start', value = 35, min = 6, max = 66),
        sliderInput('PeakEnd', label = 'Aligning peak search - end', value = 46, min = 6, max = 66),
        
        # Data processing options
        checkboxInput('zbaseline', 'Set all baselines to 0', value = TRUE),
        checkboxInput('smoothen', 'Smooth graph', value = TRUE),
        selectInput('norm', 'AUC normalization method', norm_methods),
        numericInput('maxAbs', 'Maximum absorbance value', value = 2),
        
        # PM Ratios
        h4('Polysome/Monosome ratios'), 
                   
        # Sliders for monosome and polysome areas
        checkboxInput('cutoffsPM', 'Show cutoffs for PM ratio in Processed plot', value = FALSE),
        sliderInput('MonosomeStart', label = 'Monosome area start', value = 26, min = 6, max = 66),
        sliderInput('MonosomeEnd', label = 'Monosome area end', value = 37, min = 6, max = 66),
        sliderInput('PolysomeStart', label = 'Polysome area start', value = 37, min = 6, max = 66),
        sliderInput('PolysomeEnd', label = 'Polysome area end', value = 60, min = 6, max = 70),
        
      ),
      
      # Main panel with stuff
      mainPanel(
        # Plots
        plotOutput("rawPlot"),
        plotOutput("ProfPlot"),
        # Downloads
        downloadButton('postPlotImg', 'Download High-Res image'),
        downloadButton('postPlotDtf', 'Download processed data (.csv format)'),
        downloadButton('postPlotExcel', 'Download processed data (ready to GraphPad Prism copy-paste, .xlsx)'),
        downloadButton('postPlotPrism', 'Download processed data (.pzfx) (Not recommended, slow)'),
        
        # Poly/Mono ratios
        DTOutput('PMtable'),
        downloadButton('PMcalcs', 'Download PM ratio table'),
      )
    )
)

# Server serving
server <- function(input, output) {
  
    ftable <- reactive({ 
      req(input$profiles)
      IDer(input$profiles, input$IDs)
    })
  
    output$files <- renderDT(select(ftable(), name, Sample_ID), editable = FALSE)
    
    profils1 <- reactive({
      req(input$profiles)
      profs <- list()
      for(nr in 1:nrow(ftable())){
        f <- read_csv(ftable()[nr, 'datapath'], skip = 47)
        f$Sample_ID <- ftable()[nr, 'Sample_ID']
        f <- f %>% select(`Position(mm)`, Absorbance, Sample_ID)
        profs <- c(profs, list(f))
      }
      return(profs)
    })
    
    profils2 <- reactive({
      temp <- lapply(profils1(), align, ref = profils1()[[1]], by_peaks = TRUE, npeaks = 1, minPeakPos = input$PeakStart, maxPeakPos = input$PeakEnd)
      temp <- lapply(temp, normalize, to = input$norm, max_abs = input$maxAbs, max_jump = 0.5, pos_start = 10, pos_end = 70, smoothen = input$smoothen, zero_baseline = input$zbaseline)
      return(temp)
    })
    
    plotRaw <- function(){
      p1 <- ggplot(data = rbindlist(profils1()), aes(`Position(mm)`, Absorbance, color = Sample_ID)) +
        geom_line() +
        labs(title = 'Raw plot') +
        theme_bw() +
        theme(legend.title = element_blank())
      
      if(input$cutoffsP) {p1 <- p1 +
        geom_vline(xintercept = input$PeakStart, color = 'gray30') + 
        geom_vline(xintercept = input$PeakEnd, color = 'gray20')}
      p1
    }
    
    output$rawPlot <- renderPlot({
      plotRaw()
    })
    
    plotInput <- function(){
      
      p2 <- ggplot(data = rbindlist(profils2()), aes(`Position(mm)`, Absorbance, color = Sample_ID)) +
        geom_line() + 
        labs(title = 'Processed plot') +
        theme_bw() +
        theme(legend.title = element_blank())
      
      if(input$cutoffsPM) {p2 <- p2 +
        geom_vline(xintercept = input$MonosomeStart, color = 'blue') + 
        geom_vline(xintercept = input$MonosomeEnd, color = 'blue4') +
        geom_vline(xintercept = input$PolysomeStart, color = 'red') +
        geom_vline(xintercept = input$PolysomeEnd, color = 'red3')}
     
       p2
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
    
    PMrat <- reactive(lapply(profils2(), PolyMonoRatio, mono_start = input$MonosomeStart, mono_end = input$MonosomeEnd, poly_start = input$PolysomeStart, poly_end = input$PolysomeEnd))
    
    PMtab <- reactive({
      pmt <- tibble(unique(rbindlist(profils2())$Sample_ID), unlist(PMrat()))    
      colnames(pmt) = c('ID', 'Poly/Mono ratio')
      pmt
      })
    
    output$PMtable <- renderDT(PMtab())
    
    output$PMcalcs <- downloadHandler(filename = 'PMratios.xlsx', 
                                            content = function(file) {
                                              write.xlsx(PMtab(),
                                                         file)
                                            }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
