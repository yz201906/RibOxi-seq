#!/usr/bin/env Rscript
# Title     : Visualization and Analysis of RibOxi-seq Count Data
# Objective : Shiny app that allows visualization based on gene symbol and enables various statistics.
# Created by: yinzh
# Created on: 5/4/2020
require(tidyverse)
require(DT)
require(Gviz)
require(GenomicRanges)
require(biomaRt)
require(shiny)
require(shinycssloaders)
require(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)

raw_data <- readRDS("raw_data_melt.rds")
### UI side ------------------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
    shiny::tags$h1("RibOxi-seq Counts Visualization"),
    # Sidebar+Main ------------------------------------------------------------------------------------------------------------------------------------------
    sidebarLayout(
        sidebarPanel(
            selectizeInput('samples', 'Sample:', choices = NULL, multiple = TRUE),
            numericInput('min_counts', 'Single Base Minimum counts:', 50),
#            actionButton('limit_sample', 'Show/Update Table'),
            selectizeInput('genes', 'Gene Symbol:', choices = NULL),
            selectInput(
                'genome',
                label = 'Genome:',
                choices = c("mm10", "hg19", "hg38"),
                selected = "mm10"
            ),
            actionButton('gene_sample', 'Plot'),
            width = 3
        ),
        mainPanel(
            dataTableOutput("usr_selected") %>% withSpinner(color="#0dc5c1"),
            downloadButton("download_full_table", "Download table"),
            width = 9
        )
    ),
    # Other layouts ------------------------------------------------------------------------------------------------------------------------------------------
    plotOutput("model_counts") %>% withSpinner(color="#0dc5c1"),

    fluidRow(column(1, {
        downloadButton("download_pdf1", "Download plot")
    }),
    column(1, offset = 2, {
        downloadButton("download_individual_table", "Download table for this gene")
    })),
    shiny::tags$h1(
        " ____________________________________________________________________________________________________________________________________________________________________________"
    ),

    shiny::tags$h1("                                 "),
    sidebarLayout(
        sidebarPanel(
            uiOutput("my_base_list"),
            actionButton('plot_zoomed', 'Plot'),
            downloadButton("download_pdf2", "Download plot"),
            width = 3
        ),
        mainPanel(plotOutput('zoomed_in')%>% withSpinner(color="#0dc5c1"),
                  width = 9)
    )
)

### Server side--------------------------------------------------------------------------------------------------------------------------------------

server <- function (input, output, session) {
    updateSelectizeInput(session,
                         'genes',
                         choices = unique(raw_data$gene),
                         server = TRUE)
    updateSelectizeInput(session,
                         'samples',
                         choices = unique(raw_data$sample_id),
                         server = TRUE)

    my_sample_filtered <- reactive({
        sample_filtered <- filter(raw_data, sample_id %in% input$samples)
    })

    # Table ------------------------------------------------------------------------------------------------------------------------------------------
    full_table <- reactive({
        counts_filtered <-
            subset(my_sample_filtered(), counts >= input$min_counts)
        filtered <-
            counts_filtered[order(-counts_filtered$counts), ]
        filtered
    })

    gene_table <- reactive({
        gene_range <-
            data.frame(seq(from = my_gene_to_plot_anno()@start, to = my_gene_to_plot_anno()@end))
        colnames(gene_range) <- 'base'
        gene_download <-
            left_join(gene_range, my_gene_to_plot(), by = 'base')
        gene_download <-
            gene_download %>% replace_na(list(counts = 0))
        gene_download$chr <- my_chr()
        gene_download$gene <- input$genes
        gene_download$sample_id <- input[['samples']][1]
        gene_download
    })

    # Plot components------------------------------------------------------------------------------------------------------------------------------------------
    my_ensembl <- reactive({
        if (input$genome == "mm10") {
            ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        }
        else {
            ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        }
    })

    my_gene_to_plot <- eventReactive(input$gene_sample, {
        gene_to_plot <-
            filter(my_sample_filtered(),
                   gene == input$genes &
                       sample_id == input[['samples']][1])
    })

    my_gr <- eventReactive(input$gene_sample, {
        gr <- makeGRangesFromDataFrame(
            my_gene_to_plot(),
            keep.extra.columns = TRUE,
            ignore.strand = TRUE,
            starts.in.df.are.0based = TRUE,
            seqnames.field = "chr",
            start.field = "base",
            end.field = "base"
        )
    })
    my_chr <- eventReactive(input$gene_sample, {
        chr <- as.character(unique(seqnames(my_gr())))
    })

    my_gen <- eventReactive(input$gene_sample, {
        gen <- input$genome
    })

    my_dtrack <- eventReactive(input$gene_sample, {
        dtrack <-
            DataTrack(
                my_gr(),
                name = input[['samples']][1],
                background.title = "blue",
                fontsize = 15
            )
    })

    my_gtrack <- eventReactive(input$gene_sample, {
        gtrack <- GenomeAxisTrack()
    })

    my_itrack <- eventReactive(input$gene_sample, {
        itrack <- IdeogramTrack(genome = my_gen(), chromosome = my_chr())
    })

    my_gene_to_plot_anno <- eventReactive(input$gene_sample, {
        gene_to_plot_anno <-
            BiomartGeneRegionTrack(
                biomart = my_ensembl(),
                genome = input$genome,
                symbol = input$genes
            )
    })

    my_position <- eventReactive(input$plot_zoomed, {
        position <- input[['selected_site']]
        position <- as.data.frame(position)
        position <-
            separate(position, 1, c('base', 'b', 'c'), sep = ':')
        position$base <- as.numeric(position$base)
        base_position <-
            my_gene_to_plot() %>% filter(base == position$base)
    })

    my_strack <- eventReactive(input$gene_sample, {
        if (input$genome == "mm10") {
            strack <- SequenceTrack(Mmusculus, chromosome = my_chr())
        }
        else {
            strack <- SequenceTrack(Hsapiens, chromosome = my_chr())
        }
    })

    my_grtrack <- eventReactive(input$gene_sample, {
        grtrack <-
            GeneRegionTrack(
                rstarts = start(my_gene_to_plot_anno()),
                rends = end(my_gene_to_plot_anno()),
                strand = strand(my_gene_to_plot_anno()),
                feature = feature(my_gene_to_plot_anno()),
                exon = exon(my_gene_to_plot_anno()),
                transcript = transcript(my_gene_to_plot_anno()),
                gene = gene(my_gene_to_plot_anno()),
                symbol = symbol(my_gene_to_plot_anno()),
                genome = my_gen(),
                chromosome = my_chr(),
                name = "Gene Model",
                transcriptAnnotation = "symbol",
                fontsize = 15,
                fontsize.group = 15,
                fill = 'black',
                lwd = 0.7,
                lex = 10,
                background.title = "brown"
            )
    })

    # Outputs ------------------------------------------------------------------------------------------------------------------------------------------
    output$my_base_list <- renderUI({
        for_list <-
            my_gene_to_plot() %>% unite(my_list, base, sample_id, counts, sep = ":")
        selectInput(
            'selected_site',
            'Select a site:',
            choices =  for_list$my_list,

            multiple = TRUE,
            selectize = FALSE
        )
    })

    output$usr_selected <-
        renderDT(full_table(), class = "display nowrap compact", filter = "top")

    output$model_counts <-
        renderPlot({
            plotTracks(
                list(my_itrack(), my_gtrack(), my_dtrack(), my_grtrack()),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = "histogram"
            )
        })

    output$zoomed_in <-
        renderPlot({
            plotTracks(
                list(
                    my_itrack(),
                    my_gtrack(),
                    my_dtrack(),
                    my_grtrack(),
                    my_strack()
                ),
                from = my_position()$base - 15,
                to = my_position()$base + 15,
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = "histogram"
            )
        })

    # Downloads ------------------------------------------------------------------------------------------------------------------------------------------
    output$download_full_table <- downloadHandler(
        filename = 'Filtered_data_table.csv',
        content = function(file) {
            write.csv(full_table(), file, row.names = FALSE)
        }
    )
    output$download_individual_table <- downloadHandler(
        filename = function() {
            paste(input$genes, '_data_table.csv', sep = '')
        },
        content = function(file) {
            write.csv(gene_table(), file, row.names = FALSE)
        }
    )

    output$download_pdf1 <- downloadHandler(
        filename = function() {
            paste(input$genes, '.pdf', sep = '')
        },
        content = function(file) {
            pdf(file, width = 14, height = 6)
            plotTracks(
                list(my_itrack(), my_gtrack(), my_dtrack(), my_grtrack()),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = "histogram"
            )
            dev.off()
        }
    )

    output$download_pdf2 <- downloadHandler(
        filename = function() {
            paste(input$genes, '_zoomed_in.pdf', sep = '')
        },
        content = function(file) {
            pdf(file, width = 12, height = 6)
            plotTracks(
                list(
                    my_itrack(),
                    my_gtrack(),
                    my_dtrack(),
                    my_grtrack(),
                    my_strack()
                ),
                from = my_position()$base - 15,
                to = my_position()$base + 15,
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = "histogram"
            )
            dev.off()
        }
    )
}


shinyApp(ui, server)
