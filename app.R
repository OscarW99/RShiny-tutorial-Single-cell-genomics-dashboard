
source('global.R')

ui <- dashboardPage(
    dashboardHeader(title = "scRNAseq Analysis"),
    dashboardSidebar(
        tags$head(
        tags$style(HTML(".skin-blue .main-header .sidebar-toggle {display: none;}"))
        ),
        sidebarMenu(id='tab',
            useShinyjs(),
            menuItem("Home Page", tabName = "home", icon = icon("list")),
            menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),
            conditionalPanel(condition = "input.tab == 'input'",
                div(
                    fileInput("file", "Upload File", multiple=TRUE, accept=c('.rds')),
                    actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                    actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                    )
                )
        )
    ), 
    dashboardBody(
            tabItems(
                tabItem(tabName = "input", # tabItem refers to tab in sidebar (not main panel)
                tabsetPanel(id = 'main_tabs',
                    tabPanel("Instructions",
                            includeMarkdown("./markdown/instructions.md")
                            )
                        )
                    ),
                tabItem(tabName = "home",
                tags$h1(HTML("<u>Welcome to The scRNAseq Suerat analysis RShiny app</u>")),
                )
                )
            )         
)

server <- function(input, output, session) {
    options(shiny.maxRequestSize=300*1024^2)

    values <- reactiveValues()

    # Disable Run by default
    shinyjs::disable("run")

    observe({
    if(is.null(input$file) != TRUE) {
        shinyjs::enable("run")
    } else {
        shinyjs::disable("run")
    }
    })

    observeEvent(input$run, {
        shinyjs::disable("run")

        # Clear tabs before 'Run' is ran another time
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")

        show_modal_spinner(text = "Preparing plots...")

        obj <- load_seurat_obj(input$file$datapath)
        if (is.vector(obj)){
            showModal(modalDialog(
                title = "Error with file",
                HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
                    paste(unlist(obj), collapse = "<br><br>"))
            ))
            shinyjs::enable("run")

        } else {
            
            output$umap <- renderPlot({
                if (!is.null(input$metadata_col)) {
                    create_metadata_UMAP(obj, input$metadata_col)
                }
            })

            output$featurePlot <- renderPlot({
                if (!is.null(input$gene)) {
                    create_feature_plot(obj, input$gene)
                }
            })

            output$downloadFeaturePlot <- downloadHandler(
                filename = function(){
                    paste0(input$gene, '_feature_plot', '.png')
                },
                content = function(file){
                    plot <- create_feature_plot(obj, input$gene)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )
            output$download_umap <- downloadHandler(
                filename = function(){
                    paste0(input$metadata_col, '_UMAP', '.png')
                },
                content = function(file){
                    plot <- create_metadata_UMAP(obj, input$metadata_col)
                    ggsave(filename=file, width = 10, height = 5, type = "cairo")
                }
            )

            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "UMAP",
                    fluidRow(
                    column(
                        width = 8,
                        plotOutput(outputId = 'umap'),
                        downloadButton("download_umap", "Download UMAP")
                    ),
                    column(
                        width = 4,
                        selectizeInput("metadata_col", 
                            "Metadata Column", 
                            colnames(obj@meta.data)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                ),
                select = TRUE
            )
            insertTab(
                inputId = "main_tabs",
                tabPanel(
                    "Gene Expression",
                    fluidRow(
                    column(
                        width = 8,
                        plotOutput(outputId = 'featurePlot'),
                        downloadButton("downloadFeaturePlot", "Download Feature Plot")
                    ),
                    column(
                        width = 4,
                        selectizeInput("gene", 
                            "Genes", 
                            rownames(obj)
                        )
                    )
                    ),
                    style = "height: 90%; width: 95%; padding-top: 5%;"
                )
            )

            remove_modal_spinner()
            shinyjs::enable("run")
            
        }
    })

    # Clear all sidebar inputs when 'Reset' button is clicked
    observeEvent(input$reset, {
        shinyjs::reset("file")
        removeTab("main_tabs", "UMAP")
        removeTab("main_tabs", "Gene Expression")
        shinyjs::disable("run")
    })

}

shinyApp(ui, server)


# By default, Shiny limits file uploads to 5MB per file. You can modify this limit by using the shiny.maxRequestSize option. For example, adding options(shiny.maxRequestSize=30*1024^2) to the top of server.R would increase the limit to 30MB.

# shinyjs lets you perform common useful JavaScript operations in Shiny applications without having to know any JavaScript