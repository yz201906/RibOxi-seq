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

shinyOptions(cache = diskCache(file.path(dirname(tempdir()), "riboxi_shiny_cache")))
load("raw_data.RData")

if (grepl('mm', args[3])) {
    connection<-try(ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
    if ("try-error" %in% class(connection)) ensembl <- useMart("ensembl", host="http://uswest.ensembl.org/", dataset = "mmusculus_gene_ensembl")
} else {
    connection<-try(ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl"))
    if ("try-error" %in% class(connection)) ensembl <- useMart("ensembl", host="http://uswest.ensembl.org/", dataset = "hsapiens_gene_ensembl")
}

gtrack <- GenomeAxisTrack()
my_gen <- args[3]

if (grepl('mm', args[3])) {
    my_strack <- SequenceTrack(Mmusculus)
} else {
    my_strack <- SequenceTrack(Hsapiens)
}


### UI side ------------------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
    shiny::tags$h1("RibOxi-seq Counts Visualization"),
    # Sidebar+Main ------------------------------------------------------------------------------------------------------------------------------------------
    sidebarLayout(
        sidebarPanel(
            selectizeInput('samples', 'Sample:', choices = NULL, multiple = TRUE),
            numericInput('min_counts', 'Base position min mean counts:', 50),
            selectizeInput('genes', 'Gene Symbol:', choices = NULL),
            uiOutput("my_chr_list"),
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
            sliderInput('window_size', "Select window size:", 10, 200, 30, step = 20),
            actionButton('plot_zoomed', 'Plot'),
            downloadButton("download_pdf2", "Download plot"),
            width = 2
        ),
        mainPanel(
            plotOutput('zoomed_in') %>% withSpinner(color = "#0dc5c1"),
            width = 10
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
                         choices = list(colnames(dplyr::select(raw_data, starts_with(args[1])|starts_with(args[2])))),
                         server = TRUE)

    # Table ------------------------------------------------------------------------------------------------------------------------------------------
    selected_samples <- reactive({
        sample_filtered <- dplyr::select(raw_data, chr, base, gene, seq, any_of(input$samples))
        sample_filtered
    })

    full_table <- reactive({
        sample_filtered <- mutate(selected_samples(), counts_mean=rowMeans(dplyr::select(selected_samples(), contains(args[1])|contains(args[2]))))
        counts_filtered <- filter(sample_filtered, counts_mean >= input$min_counts)
        filtered <-
            counts_filtered[order(-counts_filtered$counts_mean),]
    })


    # Plot components------------------------------------------------------------------------------------------------------------------------------------------

    my_grouping <- eventReactive(input$gene_sample, {
        sample_list<-data.frame(str_split(input$samples, ' '))
        Control_num<-dplyr::select(sample_list, contains(args[1]))
        Treatment_num<-dplyr::select(sample_list, contains(args[2]))
        grouping<-factor(c(rep(args[1], length(Control_num)), rep(args[2], length(Treatment_num))), levels = c(args[1], args[2]))
        grouping
    })

    my_gene_to_plot <- reactive({
        gene_to_plot <- filter(selected_samples(),gene == input$genes)
        gene_to_plot$gene <- NULL
        gene_to_plot$seq <- NULL
        gene_to_plot
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
        itrack <- IdeogramTrack(genome = my_gen, chromosome = input$my_chr)
    })

    my_gene_to_plot_anno <- eventReactive(input$genes, {
        gene_to_plot_anno <- BiomartGeneRegionTrack(biomart = ensembl,genome = my_gen,symbol = input$genes,)
        ranges(gene_to_plot_anno) <- subset(ranges(gene_to_plot_anno), symbol == input$genes)
        gene_to_plot_anno
    })

    my_site_list <- eventReactive(input$gene_sample, {
        site_list<-mutate(my_gene_to_plot(), counts_mean=rowMeans(dplyr::select(my_gene_to_plot(), contains(args[1])|contains(args[2]))))
        site_list<-site_list[order(-site_list$counts_mean),]
    })

    my_position <- eventReactive(input$plot_zoomed, {
        base_position <- my_gene_to_plot() %>% filter(base == input$selected_site)
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
                genome = my_gen,
                chromosome = input$my_chr,
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

    my_window_size <- reactive({input$window_size/2})

    # Outputs ------------------------------------------------------------------------------------------------------------------------------------------
    output$my_chr_list <- renderUI({
        chr_list <- as.character(unique(my_gene_to_plot()$chr))
        if (length (str_split(chr_list, ' '))>1) {
            chr_list <- str_split(chr_list, ' ')
            chr_list
        }
        selectInput('my_chr', "Chromosome:", choices = chr_list, multiple = FALSE, selectize = FALSE)
    })


    output$my_base_list <- renderUI({
        my_site_list()
        selectInput(
            'selected_site',
            'Select a site:',
            choices =  my_site_list()$base,
            multiple = TRUE,
            selectize = FALSE
        )
    })

    output$usr_selected <-
        renderDT(full_table(), class = "display nowrap compact", filter = "top")

    output$model_counts <-
        renderCachedPlot({
            plotTracks(
                list(my_itrack(), gtrack, my_dtrack(), my_grtrack()),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("p","boxplot"),
                groups = my_grouping(),
                add53 = TRUE
            )
        },
        cacheKeyExpr = {list(my_itrack(), my_dtrack(), my_grtrack(), my_grouping())}
        )

    output$zoomed_in <-
        renderCachedPlot({
            plotTracks(
                list(
                    my_itrack(),
                    gtrack,
                    my_dtrack(),
                    my_grtrack(),
                    my_strack
                ),
                from = my_position()$base - my_window_size(),
                to = my_position()$base + my_window_size(),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("p","boxplot"),
                groups = my_grouping(),
                add53 = TRUE
            )
        },
        cacheKeyExpr = {list(my_itrack(), my_dtrack(), my_grtrack(), my_position(), my_grouping(), my_window_size())}
        )

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
                type = c("p","boxplot"),
                groups = my_grouping(),
                add53 = TRUE
            )
            dev.off()
        }
    )

    output$download_pdf2 <- downloadHandler(
        filename = function() {
            paste(input$genes, '_zoomed_in.pdf', sep = '')
        },
        content = function(file) {
            pdf(file, width = 12, height = 8)
            plotTracks(
                list(
                    my_itrack(),
                    gtrack,
                    my_dtrack(),
                    my_grtrack(),
                    my_strack
                ),
                from = my_position()$base - my_window_size(),
                to = my_position()$base + my_window_size(),
                cex = 1.5,
                background.panel = "#FFFEDB",
                type = c("p","boxplot"),
                groups = my_grouping(),
                add53 = TRUE
            )
            dev.off()
        }
    )
}


shinyApp(ui, server)
