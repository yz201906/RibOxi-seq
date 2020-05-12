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
ensemblm <-useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
ensemblh <-useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
gtrack <- GenomeAxisTrack()
### UI side ------------------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
    shiny::tags$h1("RibOxi-seq Counts Visualization"),
    # Sidebar+Main ------------------------------------------------------------------------------------------------------------------------------------------
    sidebarLayout(
        sidebarPanel(
            selectizeInput('samples', 'Sample:', choices = NULL, multiple = TRUE),
            numericInput('min_counts', 'Single Base Minimum counts:', 50),
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
            dataTableOutput("usr_selected") %>% withSpinner(color = "#0dc5c1"),
            downloadButton("download_full_table", "Download table"),
            width = 9
        )
    ),
    # Other layouts ------------------------------------------------------------------------------------------------------------------------------------------
    plotOutput("model_counts") %>% withSpinner(color = "#0dc5c1"),

    fluidRow(column(1, {
        downloadButton("download_pdf1", "Download plot")
    }),
    ),
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
        mainPanel(
            plotOutput('zoomed_in') %>% withSpinner(color = "#0dc5c1"),
            width = 9
        )
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

    # Table ------------------------------------------------------------------------------------------------------------------------------------------
    full_table <- reactive({
        counts_filtered <-
            subset(raw_data, counts >= input$min_counts & sample_id %in% input$samples)
        filtered <-
            counts_filtered[order(-counts_filtered$counts),]
        filtered
    })


    # Plot components------------------------------------------------------------------------------------------------------------------------------------------
    my_ensembl <- reactive({
        if (my_gen() == "mm10") {
            ensemblm
        }
        else {
            ensemblh
        }
    })

    my_gene_to_plot <- eventReactive(input$gene_sample, {
        gene_to_plot <-
            filter(raw_data,
                   gene == input$genes & sample_id %in% input$samples)
        gene_to_plot <- gene_to_plot%>%spread(sample_id, counts)
        gene_to_plot$gene <- NULL
        gene_to_plot$seq <- NULL
        gene_to_plot

    })

    my_gr <- reactive({
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
                name = "Nm Counts",
                background.title = "blue",
                fontsize = 15
            )
    })


    my_itrack <- eventReactive(input$gene_sample, {
        itrack <- IdeogramTrack(genome = my_gen(), chromosome = my_chr())
    })

    my_gene_to_plot_anno <- eventReactive(input$gene_sample, {
        gene_to_plot_anno <-
            BiomartGeneRegionTrack(
                biomart = my_ensembl(),
                genome = my_gen(),
                symbol = input$genes
            )
    })

    my_position <- eventReactive(input$plot_zoomed, {
        base_position <- my_gene_to_plot() %>% filter(base == input$selected_site)
    })

    my_strack <- eventReactive(input$gene_sample, {
        if (my_gen() == "mm10") {
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
        selectInput(
            'selected_site',
            'Select a site:',
            choices =  my_gene_to_plot()$base,
            multiple = TRUE,
            selectize = FALSE
        )
    })

    output$usr_selected <-
        renderDT(full_table(), class = "display nowrap compact", filter = "top")

    output$model_counts <-
        renderPlot({
            plotTracks(
                list(my_itrack(), gtrack, my_dtrack(), my_grtrack()),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("confint", "p"),
                groups = rep(c("WT", "KO")), aggregateGroups = TRUE, aggregation = "mean"
            )
        })

    output$zoomed_in <-
        renderPlot({
            plotTracks(
                list(
                    my_itrack(),
                    gtrack,
                    my_dtrack(),
                    my_grtrack(),
                    my_strack()
                ),
                from = my_position()$base - 15,
                to = my_position()$base + 15,
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("confint", "p"),
                groups = rep(c("WT", "KO")), aggregation = "mean"
            )
        })

    # Downloads ------------------------------------------------------------------------------------------------------------------------------------------
    output$download_full_table <- downloadHandler(
        filename = 'Filtered_data_table.csv',
        content = function(file) {
            write.csv(full_table(), file, row.names = FALSE)
        }
    )

    output$download_pdf1 <- downloadHandler(
        filename = function() {
            paste(input$genes, '.pdf', sep = '')
        },
        content = function(file) {
            pdf(file, width = 14, height = 6)
            plotTracks(
                list(my_itrack(), gtrack, my_dtrack(), my_grtrack()),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("confint", "p"),
                groups = rep(c("WT", "KO")), aggregation = "mean"
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
                    gtrack,
                    my_dtrack(),
                    my_grtrack(),
                    my_strack()
                ),
                from = my_position()$base - 15,
                to = my_position()$base + 15,
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("confint", "p"),
                groups = rep(c("WT", "KO")), aggregation = "mean"
            )
            dev.off()
        }
    )
}


shinyApp(ui, server)
