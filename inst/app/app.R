library(shiny)
library(xsmut)
library(data.table)

options(shiny.maxRequestSize = 2000 * 1024^2)  # COSMIC files are large

ui <- fluidPage(
  titlePanel("xsmut — species mutations vs human COSMIC"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. Genes"),
      textInput("genes", "Gene symbol(s), comma-separated", "TP53"),
      h4("2. Species"),
      textInput("species", "Species (any name Ensembl knows)", "macaque"),
      textInput("assembly", "Assembly (checked against Ensembl)", "Mmul_10"),
      selectInput("human_asm", "Human assembly (must match COSMIC file)", c("GRCh38", "GRCh37")),
      h4("3. COSMIC files (local)"),
      helpText("Download from cancer.sanger.ac.uk/cosmic/download (login required)."),
      fileInput("cosmic", "Mutation export (.tsv / .tsv.gz)"),
      fileInput("census", "Cancer Gene Census (optional)"),
      h4("4. Species mutations"),
      fileInput("muts", "CSV/TSV with columns: gene, aa_change or aa_pos, [sample_id, consequence]"),
      checkboxInput("compress", "Compress introns", TRUE),
      numericInput("min_label", "Label COSMIC hotspots with ≥ N samples", 5, min = 1),
      actionButton("go", "Run", class = "btn-primary"),
      br(), br(),
      downloadButton("dl_plot", "Download plot (PDF)"),
      downloadButton("dl_table", "Download mapped mutations (CSV)")
    ),
    mainPanel(
      width = 9,
      uiOutput("gene_tabs")
    )
  )
)

server <- function(input, output, session) {

  cosmic_dt <- reactive({
    req(input$cosmic)
    genes <- trimws(strsplit(input$genes, ",")[[1]])
    withProgress(message = "Reading COSMIC…", cosmic_load(input$cosmic$datapath, genes = genes))
  })

  user_muts <- reactive({
    req(input$muts)
    fread(input$muts$datapath)
  })

  results <- eventReactive(input$go, {
    genes <- trimws(strsplit(input$genes, ",")[[1]])
    cos <- cosmic_dt(); mu <- user_muts()
    census <- if (!is.null(input$census)) input$census$datapath else NULL
    withProgress(message = "Querying Ensembl & mapping…", value = 0, {
      out <- list()
      for (g in genes) {
        incProgress(1 / length(genes), detail = g)
        out[[g]] <- tryCatch(
          xsmut_pipeline(g, mu, input$species, input$assembly, cos, census,
                         human_assembly = input$human_asm,
                         compress_introns = input$compress, min_label_n = input$min_label),
          error = function(e) list(error = conditionMessage(e)))
      }
      out
    })
  })

  output$gene_tabs <- renderUI({
    res <- results()
    tabs <- lapply(names(res), function(g) {
      r <- res[[g]]
      if (!is.null(r$error)) return(tabPanel(g, div(class = "alert alert-danger", r$error)))
      chk <- r$gene_check
      tabPanel(g,
        div(class = if (isTRUE(chk$in_census)) "alert alert-success" else "alert alert-warning",
            sprintf("%s — %s. COSMIC mutations in file: %d. Ortholog: %s (%s%% identity, %s).",
                    g,
                    if (isTRUE(chk$in_census)) paste0("in Cancer Gene Census", if (!is.na(chk$tier)) paste0(" (tier ", chk$tier, ")") else "")
                    else if (is.na(chk$in_census)) "census file not supplied" else "NOT in Cancer Gene Census",
                    chk$n_cosmic_mutations,
                    attr(r$species_mutations, "ortholog")$target$protein_id,
                    round(as.numeric(attr(r$species_mutations, "ortholog")$target$perc_id), 1),
                    attr(r$species_mutations, "ortholog")$target$type)),
        plotOutput(paste0("plot_", g), height = "650px"),
        h4("Species mutations mapped to human"),
        DT::dataTableOutput(paste0("tbl_sp_", g)),
        h4("Human COSMIC mutations"),
        DT::dataTableOutput(paste0("tbl_cos_", g))
      )
    })
    do.call(tabsetPanel, c(tabs, id = "tabs"))
  })

  observe({
    res <- results()
    for (g in names(res)) local({
      gg <- g; r <- res[[gg]]
      if (is.null(r$error)) {
        output[[paste0("plot_", gg)]] <- renderPlot(print(r$plot))
        output[[paste0("tbl_sp_", gg)]] <- DT::renderDataTable(r$species_mutations, options = list(pageLength = 10, scrollX = TRUE))
        output[[paste0("tbl_cos_", gg)]] <- DT::renderDataTable(
          r$cosmic_mutations[, .(aa_change, cds_change, consequence, chrom, start, primary_site, sample_id)],
          options = list(pageLength = 10, scrollX = TRUE))
      }
    })
  })

  output$dl_plot <- downloadHandler(
    filename = function() paste0("xsmut_", gsub("[^A-Za-z0-9]", "_", input$genes), ".pdf"),
    content = function(file) {
      res <- results(); ok <- Filter(function(r) is.null(r$error), res)
      pdf(file, width = 12, height = 8)
      for (r in ok) print(r$plot)
      dev.off()
    })

  output$dl_table <- downloadHandler(
    filename = function() "xsmut_mapped_mutations.csv",
    content = function(file) {
      res <- results(); ok <- Filter(function(r) is.null(r$error), res)
      fwrite(rbindlist(lapply(ok, `[[`, "species_mutations"), fill = TRUE), file)
    })
}

shinyApp(ui, server)
