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
library(shiny)
library(shinycssloaders)
require(BSgenome)
library(shinydashboard)

shinyOptions(cache = diskCache(file.path(dirname(tempdir()), "riboxi_shiny_cache")))
load("raw_data.RData")

if (grepl('mm', args[3])) {
    library(BSgenome.Mmusculus.UCSC.mm10)
    connection<-try(ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
    if ("try-error" %in% class(connection)) ensembl <- useMart("ensembl", host="http://uswest.ensembl.org/", dataset = "mmusculus_gene_ensembl")
} else {
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(BSgenome.Hsapiens.UCSC.hg38)
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

gene_plot <- function(itrack, dtrack, grtrack, groups, position = NULL, window_size = NULL,  font_size = 1.5) {
    if (!is.null(position) & !is.null(window_size)) {
        plotTracks(
            list(itrack, gtrack, dtrack, grtrack, my_strack),
            cex = font_size,
            background.panel = "#FFFEDB",
            type = c("p","boxplot"),
            groups = groups,
            add53 = TRUE,
            from = position$base - window_size,
            to = position$base + window_size)
    } else {
        plotTracks(
            list(itrack, gtrack, dtrack, grtrack),
            cex = font_size,
            background.panel = "#FFFEDB",
            type = c("p","boxplot"),
            groups = groups,
            add53 = TRUE,
        )
    }
}


### UI side ------------------------------------------------------------------------------------------------------------------------------------------
ui <- dashboardPage(
    dashboardHeader(title = "RibOxi-Seq End Counts Visualization", titleWidth = 400),
    # Sidebar ----
    dashboardSidebar(
        selectizeInput('samples', 'Sample:', choices = NULL, multiple = TRUE),
        numericInput('min_counts', 'Base position min mean counts:', 50),
        selectizeInput('genes', 'Gene Symbol:', choices = NULL),
        uiOutput("my_chr_list"),
        actionButton('gene_sample', 'Plot'),
        conditionalPanel(
            condition = "input.gene_sample>=1",
            h3(" "),
            downloadButton("download_pdf1", "Download plot"),
            h3("Generate zoomed-in plot"),
            uiOutput("my_base_list"),
            sliderInput('window_size', "Select window size:", 10, 200, 30, step = 20),
            actionButton('plot_zoomed', 'Plot'),
            h3(" "),
            downloadButton("download_pdf2", "Download plot")
        )
    ),
    dashboardBody(
        tabsetPanel(
            type = "tabs",
            id = "inTabset",
            tabPanel(
                "Data",
                shiny::tags$a(
                    href = "https://www.youtube.com/watch?v=UFB993xufUU",
                    "Click here for a demonstration on how counts are normalized (Median of ratios).",
                    style = "color:green;font-size:20px;"
                ),
                h1("Click any row below to select a gene and fields on the left will be filled automatically.", style = "font-size:20px;"),
                dataTableOutput("usr_selected") %>% withSpinner(color = "#0dc5c1"),
                downloadButton("download_full_table", "Download table")
            ),
            tabPanel(
                "Visualization",
                plotOutput("model_counts") %>% withSpinner(color = "#0dc5c1"),
                plotOutput('zoomed_in') %>% withSpinner(color = "#0dc5c1")
            ),
            tabPanel(
                "Summary",
                plotOutput("pca_plot"),
                downloadButton("download_pdf3", "Download plot")
                )
        )
    )
)



### Server side--------------------------------------------------------------------------------------------------------------------------------------

server <- function (input, output, session) {
    observe({
        updateSelectizeInput(session,
                             'genes',
                             choices = unique(raw_data$gene),
                             selected = clicked_gene()$gene,
                             server = TRUE)
    })
    updateSelectizeInput(session,
                         'samples',
                         choices = list(colnames(dplyr::select(raw_data, starts_with(args[1])|starts_with(args[2])))),
                         server = TRUE)
    clicked_gene <- eventReactive(input$usr_selected_rows_selected,{
        selected_gene <- full_table() %>% dplyr::slice(input$usr_selected_rows_selected)
    })
    
    observeEvent(input$usr_selected_rows_selected, {
        updateSelectInput(session, "my_chr_list", selected = clicked_gene()$chr)
    })
    
    # Table ------------------------------------------------------------------------------------------------------------------------------------------
    # selected_samples <- reactive({
    #     sample_filtered <- dplyr::select(raw_data, chr, base, gene, seq, any_of(input$samples))
    #     sample_filtered
    # })

    selected_samples <- reactive({
        sample_filtered <- dplyr::select(raw_data, chr, base, gene, any_of(input$samples))
        sample_filtered
    })
    
    full_table <- reactive({
        sample_filtered <- mutate(selected_samples(), 
                                  counts_mean=rowMeans(dplyr::select(selected_samples(), 
                                                                     contains(args[1])|contains(args[2]))))
        counts_filtered <- filter(sample_filtered, counts_mean >= input$min_counts)
        filtered <- counts_filtered[order(-counts_filtered$counts_mean),]
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
#        gene_to_plot$seq <- NULL
        gene_to_plot
    })

    my_gr <- eventReactive(input$gene_sample, {
        my_gene <- filter(my_gene_to_plot(), chr==input$my_chr)
        gr <- makeGRangesFromDataFrame(
            my_gene,
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
        gene_to_plot_anno <- BiomartGeneRegionTrack(biomart = ensembl,genome = my_gen,symbol = input$genes)
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
                background.title = "brown",
                collapseTranscripts="meta"
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
            selected = (my_site_list()$base)[1],
            selectize = FALSE
        )
    })

    output$usr_selected <-
        renderDT(full_table(), class = "display nowrap compact",
                 caption = "Please select a row.",
                 selection = 'single',
                 filter = "top")
    
    
    observeEvent(input$gene_sample, {
        output$model_counts <- renderCachedPlot({
            gene_plot(isolate(my_itrack()), isolate(my_dtrack()), 
                      isolate(my_grtrack()), isolate(my_grouping()))},
            cacheKeyExpr = {list(isolate(my_itrack()), 
                                 isolate(my_dtrack()), isolate(my_grtrack()), 
                                 isolate(my_grouping()))}
            )
    })
    observeEvent(input$plot_zoomed, {
        req(my_position())
        req(my_window_size())
        output$zoomed_in <- renderCachedPlot({
            gene_plot(isolate(my_itrack()), isolate(my_dtrack()), 
                      isolate(my_grtrack()), isolate(my_grouping()), 
                      isolate(my_position()), isolate(my_window_size()))},
            cacheKeyExpr = {list(isolate(my_itrack()), isolate(my_dtrack()), 
                                 isolate(my_grtrack()), isolate(my_position()), 
                                 isolate(my_grouping()), isolate(my_window_size()))}
            )
    })
    
    output$pca_plot <-
        renderPlot({
            pca_plot_rl_dist
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
            gene_plot(my_itrack(), my_dtrack(), my_grtrack(), my_grouping())
            dev.off()
        }
    )

    output$download_pdf2 <- downloadHandler(
        filename = function() {
            paste(input$genes, '_zoomed_in.pdf', sep = '')
        },
        content = function(file) {
            pdf(file, width = 12, height = 8)
            gene_plot(my_itrack(), my_dtrack(), my_grtrack(), my_grouping(), my_position(), my_window_size())
            dev.off()
        }
    )
    
    output$download_pdf3 <- downloadHandler(
        filename = function() {
            paste('PCA_plot.pdf', sep = '')
        },
        content = function(file) {
            ggsave(file, plot = pca_plot_rl_dist, device = "pdf")
        }
    )
    
    observeEvent(input$gene_sample, {
        updateTabsetPanel(session, "inTabset",
                          selected = "Visualization")
    })
}


shinyApp(ui, server)
