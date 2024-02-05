#' @importFrom shiny tabPanel
#' @importFrom shiny plotOutput
report_strudel <- tabPanel("Strudel plot",
   plotOutput("strudel_plot")
)

#' @importFrom shiny tabPanel
report_hot_gene <- tabPanel("Hot insert genes",
        # fluidRow(
        #     plotOutput("hot_gene_host_plot")), 
        #         style = "font-size: 70%;"),
        #     shinydashboard::box(title = "hot_gene_host_table",
        #         width = 4,
        #         collapsible = T,
        #         shinycssloaders::withSpinner(DT::dataTableOutput("hot_gene_host_table")), 
        #         style = "font-size: 70%;")),
        # fluidRow(
        #     shinydashboard::box(
        #         width = 8,
        #         collapsible = T,
        #         shinycssloaders::withSpinner(plotOutput("hot_gene_virus_plot")), 
        #         style = "font-size: 70%;"),
        #     shinydashboard::box(title = "hot_gene_virus_table",
        #         width = 4,
        #         collapsible = T,
        #         shinycssloaders::withSpinner(DT::dataTableOutput("hot_gene_virus_table")), 
        #         style = "font-size: 70%;")
        # )
        plotOutput("hot_gene_host_plot"),
        DT::dataTableOutput("hot_gene_host_table"),
        plotOutput("hot_gene_virus_plot"),
        DT::dataTableOutput("hot_gene_virus_table")
)

box1 <- box(width = 9, 
            tabBox(
            width=12,
            title = "",
            selected = "Strudel plot",
            report_strudel,
            report_hot_gene))


box2 <-  box(width = 3,
    h4("parameters for strudel plot"),
    # shinyWidgets::colorPickr("virus_color", label="virus rect color", "#EAFEFF",width=6),
    # shinyWidgets::colorPickr("host_color", label="host rect color", "#EAFEFF", width=6),
    fluidRow(
        column(width = 6, shinyWidgets::colorPickr("virus_color", label="virus rect color", "#EAFEFF")),
        column(width = 6, shinyWidgets::colorPickr("host_color", label="host rect color", "#EAFEFF"))
    ),
    fluidRow(
        column(width = 6, textInput("label_virus", label="virus label", value = "HPV16")),
        column(width = 6, textInput("label_host", label="host label", value = "Host"))
    ),
    numericInput("hot_gene",label="number of host genes", value = 5, step = 1),
    numericInput("size_gene",label="size of gene labels", value = 6, step = 1),
    numericInput("size_label",label="size of label_virus and label_host", value = 6, step = 1),
    h4("parameters for hot gene plot"),
    fluidRow(
        column(width = 6, textInput("tssRegion_left","tssRegion left:", value = "3000")),
        column(width = 6, textInput("tssRegion_right","tssRegion right:", value = "3000"))
    ),
    numericInput("hot_gene_host",label="number of host hot genes", value = 5,step = 1),
    numericInput("hot_gene_virus",label="number of virus hot genes", value = 5,step = 1)
)

#' @importFrom shinydashboard dashboardSidebar
sidebar <- dashboardSidebar(   
    textInput("accession_number","Accession number: ", value = "NC_001526.2"),
    helpText("Tip: Download the genome annotation information for the corresponding virus from NCBI."),
    # fileInput('virus_info', 'Upload virus genome annotation file',
    #                         accept=c('text/csv', 'text/comma-separated-values, text/plain')),
    # helpText("Tip: Or manually upload the virus genome annotation file. It contains at least three columns,
    #     The first column is the gene name, the second column is the start site,
    #     and the third column is the stop site."),
    # radioButtons("virus_info_select","Virus genom annotation select:",
	#                              list("NCBI" = "NCBI",
	# 			                "Upload file" = "upload"),
    #                             selected = "NCBI"),
    fileInput('insert_info', 'Upload virus inserts result file',
                            accept=c('text/csv', 'text/comma-separated-values, text/plain')),
    helpText("Tip: The virus inserts result file contains at least four columns. 
             The first column is the chromosome,
             the second column is the host insertion site,
             the third column is the virus break site, and the fourth column is the number of reads."),
    br(),
    actionButton("submit", label = "Submit",
                 style="background:#6fa6d6;color:white;
                 border: none;text-align: center;font-size: 16px;
                 font-family: 'Times New Roman', Times, serif;")
)

#' @importFrom shiny fluidRow
tab_normal <- fluidRow(
    # 改成一个tabpanel
    box1,
    box2
)

#' @importFrom shinydashboard dashboardBody
body <- dashboardBody(
    tab_normal 
)

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @noRd
app_ui <- function(request) {
    dashboardPage(skin = "green",
        dashboardHeader(title = "virusPlot"),
        sidebar,
        body
    )
}