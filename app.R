
source('global.R')
library(dplyr)
library(shinyjs)
library(shiny)
library(DT)
library(shinydashboard)
library(shinydashboardPlus)
library(ggplot2)
library(shinybusy)
library(glue)

ui <- dashboardPage(
    dashboardHeader(title = "Sngle Cell RNAseq Analysis"),
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
                fileInput("file", "Upload File", multiple=TRUE, accept=c('.Rda')),
                actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%"),
                )
            )
    )
    ),
    dashboardBody(
            tabItems(
                tabItem(tabName = "input", # tabItem refers to tab in sidebar (not main panel)
                h2("Output"),
                tabsetPanel(id = 'main_tabs',
                    tabPanel("Instructions", br(), h2("Instructions - Inputs"), br(),p(style="font-size:20px;",
                            "- Upload a Seurat object and click Submit!", br(), 
                            "- Press the 'Reset' button to clear file inputs.", br(), 
                            ))
                        )
                    ),
                tabItem(tabName = "home", # tabItem refers to tab in sidebar (not main panel)
                tags$h1(HTML("<u>Welcome to The scRNAseq Suerat analysis RShiny app</u>")),
                )
                )
            )         
            
)

server <- function(input, output, session) {

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
}

shinyApp(ui, server)