

#' @import shiny
box1_strudel <- box(width = 9,
    uiOutput("strudel_plot_ui"),
    fluidRow(
        column(width = 3, selectInput("format1","Format", list("pdf", "jpg", "png", "tiff"), selected = "pdf")),
        column(width = 3, numericInput("dpi1","Dpi", value = 300,step = 10)),
        column(width = 3, numericInput("w1","Width", value = 800,step = 10)),
        column(width = 3, numericInput("h1","Height", value = 450,step = 10))
    ),
    fluidRow(    
        column(width = 6, downloadButton("down1","Download"))
    )
)

#' @import shiny
box1_hotgene <- box(width = 9,
    uiOutput("hot_gene_host_plot_ui"),
    fluidRow(
        column(width = 3, selectInput("format2","Format", list("pdf", "jpg", "png", "tiff"), selected = "pdf")),
        column(width = 3, numericInput("dpi2","Dpi", value = 300,step = 10)),
        column(width = 3, numericInput("w2","Width", value = 800,step = 10)),
        column(width = 3, numericInput("h2","Height", value = 450,step = 10))
    ),
    fluidRow(    
        column(width = 6, downloadButton("down2","Download"))
    ),
    # plotOutput("hot_gene_host_plot"),
    DT::dataTableOutput("hot_gene_host_table"),
    uiOutput("hot_gene_virus_plot_ui"),
    fluidRow(
        column(width = 3, selectInput("format3","Format", list("pdf", "jpg", "png", "tiff"), selected = "pdf")),
        column(width = 3, numericInput("dpi3","Dpi", value = 300,step = 10)),
        column(width = 3, numericInput("w3","Width", value = 800,step = 10)),
        column(width = 3, numericInput("h3","Height", value = 450,step = 10))
    ),
    fluidRow(    
        column(width = 6, downloadButton("down3","Download"))
    ),
    # plotOutput("hot_gene_virus_plot"),
    dataTableOutput("hot_gene_virus_table")
)


box2_strudel <- box(width = 3,
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
    numericInput("size_label",label="size of label_virus and label_host", value = 6, step = 1)
)

box2_hotgene <- box(width = 3,
    h4("parameters for hot gene plot"),
    fluidRow(
        column(width = 6, shinyWidgets::colorPickr("observed_color", label="Observed color", "#4d4d4d")),
        column(width = 6, shinyWidgets::colorPickr("expected_color", label="Expected color", "#999999"))
    ),
    fluidRow(
        column(width = 6, textInput("tssRegion_left","tssRegion left:", value = "3000")),
        column(width = 6, textInput("tssRegion_right","tssRegion right:", value = "3000"))
    ),
    numericInput("hot_gene_host",label="number of host hot genes", value = 5,step = 1),
    numericInput("hot_gene_virus",label="number of virus hot genes", value = 8,step = 1)
)


#' @import shiny
report_strudel <- tabPanel("Strudel plot",
    fluidRow(
        box1_strudel,
        box2_strudel
    )
)

#' @importFrom shiny tabPanel
report_hot_gene <- tabPanel("Hot insert genes",
    fluidRow(
        box1_hotgene,
        box2_hotgene
    )
)

body <- box(width = 12, 
            tabBox(
            width=12,
            title = "",
            selected = "Strudel plot",
            report_strudel,
            report_hot_gene))


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
        tags$div(
      style = "height: 5cm; overflow-y: auto; overflow-x: auto; background-color: white;",
      span(verbatimTextOutput("fileContent"), style="color:black")
    ),
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


#' @importFrom shinydashboard dashboardBody
body <- dashboardBody(
    body 
)

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shinydashboard
#' @noRd
app_ui <- function(request) {
    dashboardPage(skin = "green",
        dashboardHeader(title = "virusPlot"),
        sidebar,
        body
    )
}