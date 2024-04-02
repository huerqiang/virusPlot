#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @importFrom shiny observeEvent
#' @importFrom shiny renderPlot
#' @importFrom shiny withProgress
# @importFrom plotly renderPlotly
# @importFrom DT dataTableOutput
# @importFrom DT renderDataTable
#' @importFrom plotly plotlyOutput
#' @noRd
app_server <- function(input, output, session) {
    output$fileContent <- renderText({
        if (is.null(input$insert_info)) {
          content <- readLines(system.file("extdata", "insert_info.txt", package = "virusPlot"))
          paste(content, collapse = "\n")
        } else {
          inFile <- input$insert_info
          content <- readLines(inFile$datapath)
          paste(content, collapse = "\n")
        }
    })
    
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
            # insert_info <- read.table(input$insert_info$datapath[1], sep = "\t", header = TRUE)
            if (is.null(input$insert_info)) {
              data(insert_info, package = "virusPlot")
              insert_info2 <- insert_info
            } else {
              insert_info2 <- read.table(input$insert_info$datapath[1], sep = "\t", header = TRUE)
            }
            virus_color <- input$virus_color
            host_color <- input$host_color
            label_virus <- input$label_virus
            label_host <- input$label_host
            hot_gene  <- as.numeric(input$hot_gene)
            size_gene <- as.numeric(input$size_gene)
            size_label <- as.numeric(input$size_label)
            
            p_strudel <- strudel_plot(virus_info = virus_info, insert_info = insert_info2,
                                  virus_color = virus_color, host_color = host_color,
                                  label_virus = label_virus, label_host = label_host, 
                                  hot_gene = hot_gene, size_gene = size_gene, 
                                  size_label = size_label) 
            output$strudel_plot <- renderPlot({
                p_strudel 
            }) 

            output$strudel_plot_ui <- renderUI({
                ns <- session$ns
                plotOutput("strudel_plot",width = paste0(input$w1, "px"),
                           height = paste0(input$h1, "px"))
            })
        })


        # hot genes
        withProgress({
            incProgress(message = "Ploting hot genes")
            tssRegion_left <- input$tssRegion_left |> as.numeric()
            tssRegion_right <- input$tssRegion_right |> as.numeric()
            hot_gene_host <- input$hot_gene_host
            hot_gene_virus <- input$hot_gene_virus
            observed_color <- input$observed_color
            expected_color <- input$expected_color
            tssRegion = c(-tssRegion_left, tssRegion_right)
            hot_gene <- get_hot_gene(virus_info, insert_info2, tssRegion = tssRegion)
            insert_plot <- hot_gene_plot(hot_gene, hot_gene_host = hot_gene_host, 
                hot_gene_virus = hot_gene_virus, observed_color = observed_color, 
                expected_color = expected_color)
            p_host <- insert_plot$p_host
            p_virus <- insert_plot$p_virus
            t_host <- hot_gene$host
            t_virus <- hot_gene$virus
 
            output$hot_gene_host_plot <- renderPlot({
                p_host
            })

            output$hot_gene_host_plot_ui <- renderUI({
                ns <- session$ns
                plotOutput("hot_gene_host_plot",width = paste0(input$w2, "px"),
                           height = paste0(input$h2, "px"))
            })
            output$hot_gene_virus_plot <- renderPlot({
                p_virus
            })
            output$hot_gene_virus_plot_ui <- renderUI({
                ns <- session$ns
                plotOutput("hot_gene_virus_plot",width = paste0(input$w3, "px"),
                           height = paste0(input$h3, "px"))
            })
            t_host$gene <- sapply(t_host$gene, function(gene) {  
                   url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gene, "&keywords=", gene)  
                   HTML(paste0("<a href='", url, "' target='_blank'>", gene, "</a>"))  
            }) 
            #######
            # save(t_host, file = "E:\\linshi\\2024_3_31\\test_shiny\\t_host.Rdata")
            ########
            output$hot_gene_host_table <- DT::renderDataTable({
                # t_host
                DT::datatable(t_host,   
                      extensions = 'Responsive', # 可选，使表格响应式  
                      options = list(dom = 'Bfrtip', # 可选，定义表格控件的布局  
                                     columnDefs = list(list(visible = FALSE, targets = 0))), # 隐藏第一列（原始的Gene列）  
                      escape = FALSE) # 允许HTML渲染  
            })
    
            output$hot_gene_virus_table <- renderDataTable({
                t_virus
            })
        })

        output$down1 <- downloadHandler(
            filename = function(){
              paste0("strudel_plot_",Sys.Date(),".",input$format1)
            },
            content = function(file){
                ggplot2::ggsave(plot = p_strudel, file = file, width = input$w1/72,
                                height = input$h1/72, dpi = input$dpi1)
                
            }
        )
        output$down2 <- downloadHandler(
            filename = function(){
              paste0("strudel_plot_",Sys.Date(),".",input$format2)
            },
            content = function(file){
                ggplot2::ggsave(plot = p_host, file = file, width = input$w2/72,
                                height = input$h2/72, dpi = input$dpi2)
                
            }
        )
        output$down3 <- downloadHandler(
            filename = function(){
              paste0("strudel_plot_",Sys.Date(),".",input$format3)
            },
            content = function(file){
                ggplot2::ggsave(plot = p_virus, file = file, width = input$w3/72,
                                height = input$h3/72, dpi = input$dpi3)   
            }
        )       
    })
}

