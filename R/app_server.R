#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @importFrom shiny observeEvent
#' @importFrom shiny renderPlot
#' @importFrom shiny withProgress
# @importFrom plotly renderPlotly
#' @importFrom DT dataTableOutput
#' @importFrom DT renderDataTable
#' @importFrom plotly plotlyOutput
#' @noRd
app_server <- function(input, output, session) {
    observeEvent(input$submit, {  
        withProgress({
            incProgress(message = "Downloading virus annotation data")
            accession_number <- input$accession_number
            gene_features <- get_virus_annotation(accession_number = accession_number,
                    email = "13766876214@163.com")
            virus_info <- deal_virus_annotation(gene_features)
        })
        # strudel plot      
        # virus_info_select <- input$virus_info_select
        # if (virus_info_select == "NCBI") {
        #     accession_number <- input$accession_number
        #     gene_features <- get_virus_annotation(accession_number = accession_number,
        #         email = "13766876214@163.com")
        #     virus_info <- deal_virus_annotation(gene_features)
        # } else {
        #     virus_info <- read.table(input$virus_info$datapath[1], sep = "\t", header = TRUE)      
        # }
        withProgress({
            incProgress(message = "Ploting strudel plot")
            insert_info <- read.table(input$insert_info$datapath[1], sep = "\t", header = TRUE)
            virus_color <- input$virus_color
            host_color <- input$host_color
            label_virus <- input$label_virus
            label_host <- input$label_host
            hot_gene  <- as.numeric(input$hot_gene)
            size_gene <- as.numeric(input$size_gene)
            size_label <- as.numeric(input$size_label)
            
            p_strudel <- strudel_plot(virus_info = virus_info, insert_info = insert_info,
                                  virus_color = virus_color, host_color = host_color,
                                  label_virus = label_virus, label_host = label_host, 
                                  hot_gene = hot_gene, size_gene = size_gene, 
                                  size_label = size_label) 
            output$strudel_plot <- renderPlot({
                p_strudel 
            }) 
        })


        # hot genes
        withProgress({
            incProgress(message = "Ploting hot genes")
            tssRegion_left <- input$tssRegion_left |> as.numeric()
            tssRegion_right <- input$tssRegion_right |> as.numeric()
            hot_gene_host <- input$hot_gene_host
            hot_gene_virus <- input$hot_gene_virus
            tssRegion = c(-tssRegion_left, tssRegion_right)
            hot_gene <- get_hot_gene(virus_info, insert_info, tssRegion = tssRegion)
            insert_plot <- hot_gene_plot(hot_gene, hot_gene_host = hot_gene_host, hot_gene_virus = hot_gene_virus)
            p_host <- insert_plot$p_host
            p_virus <- insert_plot$p_virus
            t_host <- hot_gene$host
            t_virus <- hot_gene$virus
            output$hot_gene_host_plot <- renderPlot({
                p_host
            })
            output$hot_gene_virus_plot <- renderPlot({
                p_virus
            })
    
            output$hot_gene_host_table <- renderDataTable({
                t_host
            })
    
            output$hot_gene_virus_table <- renderDataTable({
                t_virus
            })
        })

    })
}

