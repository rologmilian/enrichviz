# app.R
# Heatmap, Bar/Bubble Plot, Chord Diagram, Boxplot/Violin visualization
# v.1.3.0  — adds a molecule frequency per pathway plot
# Copyright 2026 RGM
# MIT License — see README.md for full license text
# ─────────────────────────────────────────────────────────

library(shiny)
library(pheatmap)
library(tidyverse)
library(circlize)

options(shiny.maxRequestSize = 500 * 1024^2)

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  
  tags$style(HTML("
    html, body { height: 100%; overflow: hidden; }

    .col-sm-3 {
      position : fixed; top: 60px; left: 0;
      height: calc(100vh - 60px);
      overflow-y: auto; overflow-x: hidden;
      padding: 10px 15px; z-index: 100;
    }
    .col-sm-9 {
      position: fixed; top: 60px; left: 25%; width: 75%;
      height: calc(100vh - 60px);
      overflow-y: auto; overflow-x: hidden;
      padding: 10px 15px;
    }

    .container-fluid > h2 {
      position: fixed; top: 0; left: 0; width: 100%;
      background: white; z-index: 200;
      padding: 10px 20px 5px 20px; margin: 0;
      border-bottom: 1px solid #ddd;
    }

    #app_blurb{
      position: fixed;
      top: 80px;
      left: 0;
      width: 100%;
      background: white;
      z-index: 1000;
      padding: 8px 20px;
      border-bottom: 1px solid #eee;
      color: #444;
      font-size: 14px;
      line-height: 1.35;
    }

    .col-sm-3{
      top: 170px !important;
      height: calc(100vh - 170px) !important;
    }
    .col-sm-9{
      top: 170px !important;
      height: calc(100vh - 170px) !important;
    }

    .col-sm-3::-webkit-scrollbar,
    .col-sm-9::-webkit-scrollbar { width: 6px; }
    .col-sm-3::-webkit-scrollbar-thumb,
    .col-sm-9::-webkit-scrollbar-thumb {
      background: #cccccc; border-radius: 3px;
    }

    /* filter panel highlight */
    #filter_panel {
      background: #f0f7ff;
      border: 1px solid #b8d4f0;
      border-radius: 6px;
      padding: 10px 12px 6px 12px;
      margin-bottom: 6px;
    }
    #filter_panel h4 { margin-top: 4px; color: #1a4f8a; }
  ")),
  
  titlePanel("Heatmap, Bar Plot and Chord Diagram of Normalized Data by Pathway or Function"),
  
  div(
    id = "app_blurb",
    HTML(
      "This app integrates normalized expression/abundance data with sample metadata and pathway/function enrichment results to generate bar/bubble plots, chord diagrams, heatmaps, and box/violin plots.<br>
       <b>Files needed (CSV):</b> (1) Normalized data (rows = genes/proteins; columns = samples), (2) metadata/annotation (sample IDs + group/condition), (3) enrichment results (pathway/function name, p-value/FDR or -log10 value, and member molecules/genes)."
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Upload Data Files"),
      fileInput("file_norm_prot",  "Normalized Data (.csv)",               accept = ".csv"),
      fileInput("file_annotation", "Metadata or Sample Annotation (.csv)", accept = ".csv"),
      fileInput("file_ipa_funct",  "Enriched Pathway or Functions (.csv)", accept = ".csv"),
      
      hr(),
      
      #  Metadata slice filter 
      div(
        id = "filter_panel",
        h4("\U0001f50d Filter Samples by Metadata"),
        helpText("Select a metadata column and keep only the desired values.
                  All plots will use only the matching samples."),
        uiOutput("ui_filter_col"),
        uiOutput("ui_filter_vals"),
        uiOutput("ui_filter_summary")
      ),
      
      hr(),
      
      h4("Select Pathway / Functions \u2014 Column"),
      helpText("Category = col 1. Molecules = last column by default."),
      uiOutput("ui_col_category"),
      uiOutput("ui_col_molecules"),
      
      radioButtons("mol_separator",
                   label    = "Gene/molecule separator in pathway file",
                   choices  = c("Comma  ( , )"     = ",",
                                "Slash  ( / )"     = "/",
                                "Semicolon  ( ; )" = ";",
                                "Pipe  ( | )"      = "|"),
                   selected = ",",
                   inline   = TRUE),
      
      hr(),
      
      h4("Normalized counts or protein abundance \u2014 Column"),
      helpText("Gene/Protein ID column (e.g. 'gene_symbol')."),
      uiOutput("ui_col_gene"),
      
      hr(),
      
      h4("Sample Annotation \u2014 Column"),
      uiOutput("ui_col_sample_id"),
      uiOutput("ui_col_group"),
      
      hr(),
      
      h4("Bar / Bubble Plot Settings"),
      
      radioButtons("barplot_type",
                   label    = "Plot type",
                   choices  = c("Bar plot" = "bar", "Bubble plot" = "bubble"),
                   selected = "bar",
                   inline   = TRUE),
      
      uiOutput("ui_col_pvalue"),
      uiOutput("ui_pvalue_transform_info"),
      numericInput("bar_top_n", "Show top N pathways", value = 20, min = 5, max = 200),
      
      conditionalPanel(
        condition = "input.barplot_type == 'bar'",
        textInput("bar_color", "Bar fill colour (hex or name)", value = "steelblue")
      ),
      
      conditionalPanel(
        condition = "input.barplot_type == 'bubble'",
        uiOutput("ui_bubble_size_col"),
        uiOutput("ui_bubble_color_col"),
        sliderInput("bubble_size_range",
                    label = "Bubble size range (min / max)",
                    min = 1, max = 30, value = c(3, 15), step = 1),
        selectInput("bubble_palette",
                    label   = "Colour palette",
                    choices = c("Blues", "Reds", "Purples", "Greens",
                                "OrRd", "YlOrRd", "RdYlBu", "viridis",
                                "magma", "plasma"),
                    selected = "Blues")
      ),
      
      numericInput("bar_plot_height", "Plot height (px)", value = 600, min = 300, max = 2000),
      downloadButton("download_barplot", "Download Plot (.png)"),
      
      hr(),
      
      h4("Chord Diagram Settings"),
      numericInput("chord_top_n",           "Top N pathways to include",            value = 10,   min = 2,   max = 50),
      numericInput("chord_min_freq",        "Minimum times a molecule must appear", value = 2,    min = 1,   max = 20),
      helpText("Increase 'minimum appearances' to reduce clutter."),
      sliderInput("chord_inner_radius",     "Inner circle size",                    min = 0.1, max = 0.7, value = 0.4, step = 0.05),
      helpText("Reduce to shrink central hole and give more space to labels."),
      numericInput("chord_label_size_path", "Pathway label font size",              value = 0.55, min = 0.3, max = 1.5, step = 0.05),
      numericInput("chord_label_size_mol",  "Molecule label font size",             value = 0.70, min = 0.3, max = 1.5, step = 0.05),
      textInput("chord_pathway_color",      "Pathway sector colour (hex or name)",  value = "tomato"),
      textInput("chord_molecule_color",     "Molecule sector colour (hex or name)", value = "steelblue"),
      numericInput("chord_plot_height",     "Chord plot height (px)",               value = 750,  min = 400, max = 2000),
      downloadButton("download_chord", "Download Chord Diagram (.png)"),
      
      hr(),
      
      h4("Gene Frequency Plot Settings"),
      radioButtons(
        "freq_pathway_scope",
        "Pathways used to compute gene frequency",
        choices  = c("All pathways" = "all", "Top N pathways (by significance)" = "top"),
        selected = "all",
        inline   = FALSE
      ),
      conditionalPanel(
        condition = "input.freq_pathway_scope == 'top'",
        numericInput("freq_top_n_pathways",
                     "Top N pathways",
                     value = 20, min = 2, max = 200)
      ),
      numericInput("freq_min_pathways",
                   "Minimum # pathways a gene must appear in",
                   value = 5, min = 1, max = 100),
      numericInput("freq_top_n_genes",
                   "Show top N genes (after filtering)",
                   value = 50, min = 5, max = 500),
      numericInput("freq_plot_height",
                   "Frequency plot height (px)",
                   value = 600, min = 300, max = 2000),
      downloadButton("download_gene_freq", "Download Frequency Plot (.png)"),
      
      hr(),
      
      h4("Heatmap Settings"),
      uiOutput("category_selector"),
      hr(),
      uiOutput("ui_group_colors"),
      numericInput("plot_height", "Heatmap height (px)", value = 700, min = 300, max = 2000),
      
      h4("Download Heatmaps"),
      downloadButton("download_single", "Download Selected Heatmap (.png)"),
      br(), br(),
      downloadButton("download_all",    "Download ALL Heatmaps (.zip)"),
      
      hr(),
      
      h4("Boxplot / Violin Settings"),
      radioButtons("plot_type",
                   label    = "Plot type",
                   choices  = c("Boxplot" = "box", "Violin plot" = "violin"),
                   selected = "box",
                   inline   = TRUE),
      uiOutput("ui_boxplot_protein"),
      helpText("Group colours are shared with Heatmap Settings above."),
      numericInput("boxplot_height", "Plot height (px)", value = 500, min = 300, max = 2000),
      downloadButton("download_boxplot_sidebar", "Download Plot (.png)")
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("Bar Plot",
                 br(),
                 verbatimTextOutput("barplot_status"),
                 plotOutput("barplot", height = "600px"),
                 br(),
                 downloadButton("download_barplot_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("Chord Diagram",
                 br(),
                 verbatimTextOutput("chord_status"),
                 plotOutput("chord_plot", height = "750px"),
                 br(),
                 downloadButton("download_chord_main", "\u2b07 Download Chord Diagram (.png)")
        ),
        
        tabPanel("Gene Frequency",
                 br(),
                 verbatimTextOutput("gene_freq_status"),
                 plotOutput("gene_freq_plot", height = "600px"),
                 br(),
                 downloadButton("download_gene_freq_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("Heatmap",
                 br(),
                 verbatimTextOutput("status"),
                 plotOutput("heatmap_plot", height = "700px"),
                 br(),
                 downloadButton("download_single_main", "\u2b07 Download Selected Heatmap (.png)"),
                 br(), br(),
                 downloadButton("download_all_main",    "\u2b07 Download ALL Heatmaps (.zip)")
        ),
        
        tabPanel("Boxplot",
                 br(),
                 verbatimTextOutput("boxplot_status"),
                 plotOutput("boxplot_plot", height = "500px"),
                 br(),
                 downloadButton("download_boxplot_main", "\u2b07 Download Plot (.png)")
        )
      )
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  # ── 1. Load uploaded data ─────────────────────────────────────────────────
  norm_prot  <- reactive({ req(input$file_norm_prot);  read.csv(input$file_norm_prot$datapath,  check.names = FALSE) })
  ipa_funct  <- reactive({ req(input$file_ipa_funct);  read.csv(input$file_ipa_funct$datapath,  check.names = FALSE) })
  annotation_raw <- reactive({ req(input$file_annotation); read.csv(input$file_annotation$datapath, check.names = FALSE) })
  
  # ── 2. Metadata filter widgets ────────────────────────────────────────────
  
  # Column selector for filtering (excludes the sample-ID column)
  output$ui_filter_col <- renderUI({
    req(annotation_raw())
    cols <- colnames(annotation_raw())
    selectInput("filter_col",
                label   = "Filter column",
                choices = c("(none)" = "__none__", cols),
                selected = "__none__")
  })
  
  # Value checkboxes — populated from the chosen filter column
  output$ui_filter_vals <- renderUI({
    req(annotation_raw(), input$filter_col)
    if (is.null(input$filter_col) || input$filter_col == "__none__") return(NULL)
    
    col    <- input$filter_col
    vals   <- sort(unique(as.character(annotation_raw()[[col]])))
    # Remove NA string representations
    vals   <- vals[!vals %in% c("NA", "NaN", "")]
    
    tagList(
      checkboxGroupInput(
        inputId  = "filter_vals",
        label    = paste0("Keep values in '", col, "'"),
        choices  = vals,
        selected = vals          # all selected = no filtering by default
      ),
      actionButton("filter_select_all",  "Select all",  class = "btn-xs btn-default"),
      actionButton("filter_select_none", "Clear all",   class = "btn-xs btn-default")
    )
  })
  
  # Select-all / Clear-all buttons
  observeEvent(input$filter_select_all, {
    req(annotation_raw(), input$filter_col)
    if (input$filter_col == "__none__") return()
    vals <- sort(unique(as.character(annotation_raw()[[input$filter_col]])))
    vals <- vals[!vals %in% c("NA", "NaN", "")]
    updateCheckboxGroupInput(session, "filter_vals", selected = vals)
  })
  
  observeEvent(input$filter_select_none, {
    updateCheckboxGroupInput(session, "filter_vals", selected = character(0))
  })
  
  # The filtered annotation that all plots use
  # ── KEY reactive ─────────────────────────────────────────────────────────
  annotation <- reactive({
    req(annotation_raw())
    ann <- annotation_raw()
    
    # No filter column chosen → return everything
    if (is.null(input$filter_col) || input$filter_col == "__none__")
      return(ann)
    
    # Filter column chosen but no values selected → warn but return empty
    req(input$filter_vals)
    
    col  <- input$filter_col
    keep <- as.character(ann[[col]]) %in% input$filter_vals
    ann[keep, , drop = FALSE]
  })
  
  # Summary badge shown in the filter panel
  output$ui_filter_summary <- renderUI({
    req(annotation_raw())
    n_total    <- nrow(annotation_raw())
    n_filtered <- nrow(annotation())
    if (is.null(input$filter_col) || input$filter_col == "__none__") {
      helpText(paste0("\u2139\ufe0f  Using all ", n_total, " samples (no filter active)."))
    } else {
      col <- input$filter_col
      helpText(paste0(
        "\u2705  ", n_filtered, " / ", n_total, " samples retained",
        " \u2014 ", col, ": ",
        paste(input$filter_vals, collapse = ", ")
      ))
    }
  })
  
  # ── 3. Column-selector widgets ────────────────────────────────────────────
  output$ui_col_category <- renderUI({
    req(ipa_funct())
    cols <- colnames(ipa_funct())
    selectInput("col_category", label = paste0("Category column (detected: '", cols[1], "')"),
                choices = cols, selected = cols[1])
  })
  
  output$ui_col_molecules <- renderUI({
    req(ipa_funct())
    cols     <- colnames(ipa_funct())
    keywords <- c("geneid","gene_id","gene","molecules","mol","members","genes")
    detected <- cols[sapply(tolower(cols), function(x)
      any(sapply(keywords, function(k) grepl(k, x, fixed = TRUE))))]
    default  <- if (length(detected) > 0) detected[1] else cols[length(cols)]
    selectInput("col_molecules", label = paste0("Molecules column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  output$ui_col_gene <- renderUI({
    req(norm_prot())
    cols     <- colnames(norm_prot())
    keywords <- c("gene_symbol", "gene.symbol", "genesymbol",
                  "symbol", "gene_name", "genename",
                  "gene", "protein", "uniprot", "id", "name")
    detected <- cols[sapply(tolower(cols), function(x)
      any(sapply(keywords, function(k) grepl(k, x, fixed = TRUE))))]
    default  <- if (length(detected) > 0) detected[1] else cols[1]
    selectInput("col_gene",
                label   = paste0("Gene/Protein ID column (detected: '", default, "')"),
                choices = cols,
                selected = default)
  })
  
  output$ui_col_sample_id <- renderUI({
    req(annotation())
    cols <- colnames(annotation())
    selectInput("col_sample_id", label = paste0("Sample ID column (detected: '", cols[1], "')"),
                choices = cols, selected = cols[1])
  })
  
  output$ui_col_group <- renderUI({
    req(annotation())
    cols    <- colnames(annotation())
    default <- if (length(cols) >= 2) cols[2] else cols[1]
    selectInput("col_group", label = paste0("Group column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  # ── 4. P-value column selector ────────────────────────────────────────────
  output$ui_col_pvalue <- renderUI({
    req(ipa_funct())
    cols     <- colnames(ipa_funct())
    keywords <- c("fdr","adj","p\\.val","pval","p-val","p_val",
                  "p\\.adjust","padj","log\\(p","-log","bh","bonf")
    detected <- cols[sapply(cols, function(x)
      any(sapply(keywords, function(k) grepl(k, x, ignore.case = TRUE, perl = TRUE))))]
    default  <- if (length(detected) > 0) detected[1] else cols[2]
    selectInput("col_pvalue", label = "P-value or -log10(p-value) column (X-axis)",
                choices = cols, selected = default)
  })
  
  # ── 5. Transformation info ────────────────────────────────────────────────
  pvalue_transform <- reactive({
    req(ipa_funct(), input$col_pvalue)
    raw <- suppressWarnings(as.numeric(ipa_funct()[[input$col_pvalue]]))
    raw <- raw[!is.na(raw)]
    if (length(raw) == 0) return("unknown")
    if (all(raw >= 0 & raw <= 1)) "raw" else "log10"
  })
  
  output$ui_pvalue_transform_info <- renderUI({
    req(pvalue_transform())
    if (pvalue_transform() == "raw")
      helpText("\u2714 Raw p-values detected (0-1). Will apply -log10 transformation.")
    else if (pvalue_transform() == "log10")
      helpText("\u2714 Already -log10 transformed values detected. Will use as-is.")
    else
      helpText("\u26a0 Could not determine transformation. Check column selection.")
  })
  
  # ── 6. Dynamic colour inputs: one per group (uses filtered annotation) ────
  output$ui_group_colors <- renderUI({
    req(annotation(), input$col_group)
    groups <- unique(as.character(annotation()[[input$col_group]]))
    default_colors <- c("purple","darkorange","steelblue","forestgreen",
                        "firebrick","gold","hotpink","cyan4","sienna","gray40")
    tagList(lapply(seq_along(groups), function(i) {
      textInput(inputId = paste0("color_group_", i),
                label   = paste0("Colour \u2014 ", groups[i], " (hex or name)"),
                value   = default_colors[i])
    }))
  })
  
  # ── 7. Sample columns resolved (uses filtered annotation) ─────────────────
  sample_cols_resolved <- reactive({
    req(norm_prot(), annotation(), input$col_sample_id)
    ann_ids   <- as.character(annotation()[[input$col_sample_id]])
    data_cols <- colnames(norm_prot())
    data_cols[data_cols %in% ann_ids]
  })
  
  # ── 8. Category dropdown (heatmap) ───────────────────────────────────────
  output$category_selector <- renderUI({
    req(ipa_funct(), input$col_category)
    choices <- unique(ipa_funct()[[input$col_category]])
    selectInput("selected_category", label = paste0("Select: ", input$col_category),
                choices = choices, selected = choices[1])
  })
  
  # ── 9. Gene list from selected category ──────────────────────────────────
  genes_interest <- reactive({
    req(ipa_funct(), input$selected_category, input$col_category,
        input$col_molecules, input$mol_separator)
    ipa_funct() %>%
      filter(.data[[input$col_category]] == input$selected_category) %>%
      pull(.data[[input$col_molecules]]) %>%
      strsplit(input$mol_separator, fixed = TRUE) %>%
      unlist() %>% trimws() %>% unique()
  })
  
  # ── 10. Heatmap status ────────────────────────────────────────────────────
  output$status <- renderText({
    req(genes_interest(), norm_prot(), input$col_gene, sample_cols_resolved())
    matched <- sum(as.character(norm_prot()[[input$col_gene]]) %in% genes_interest())
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter  : ", input$filter_col, " = ",
             paste(input$filter_vals, collapse = ", "), "\n")
    else ""
    paste0(
      filter_line,
      "Samples used     : ", length(sample_cols_resolved()),               "\n",
      "Category         : ", input$selected_category,                      "\n",
      "Genes in list    : ", length(genes_interest()),                      "\n",
      "Proteins matched : ", matched,                                       "\n",
      "Sample columns   : ", paste(sample_cols_resolved(), collapse = ", ")
    )
  })
  
  # ══════════════════════════════════════════════════════════════════════════
  # BAR / BUBBLE PLOT
  # ══════════════════════════════════════════════════════════════════════════
  output$ui_bubble_size_col <- renderUI({
    req(ipa_funct())
    cols     <- colnames(ipa_funct())
    keywords <- c("count","size","gene","overlap","number","n_gene","ngene","ratio")
    detected <- cols[sapply(tolower(cols), function(x)
      any(sapply(keywords, function(k) grepl(k, x, fixed = TRUE))))]
    default  <- if (length(detected) > 0) detected[1] else cols[1]
    selectInput("bubble_size_col",
                label   = paste0("Bubble SIZE column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  output$ui_bubble_color_col <- renderUI({
    req(ipa_funct())
    cols     <- colnames(ipa_funct())
    keywords <- c("fdr","adj","p\\.val","pval","p-val","p_val",
                  "p\\.adjust","padj","log\\(p","-log","bh","bonf")
    detected <- cols[sapply(cols, function(x)
      any(sapply(keywords, function(k) grepl(k, x, ignore.case = TRUE, perl = TRUE))))]
    default  <- if (length(detected) > 0) detected[1] else cols[2]
    selectInput("bubble_color_col",
                label   = paste0("Bubble COLOUR intensity column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  barplot_data <- reactive({
    req(ipa_funct(), input$col_category, input$col_pvalue, pvalue_transform())
    df       <- ipa_funct()
    raw_vals <- suppressWarnings(as.numeric(df[[input$col_pvalue]]))
    
    df_out <- df %>%
      mutate(
        .raw       = raw_vals,
        .neg_log10 = if (pvalue_transform() == "raw") -log10(.raw) else .raw
      ) %>%
      filter(!is.na(.neg_log10), is.finite(.neg_log10)) %>%
      arrange(desc(.neg_log10)) %>%
      slice_head(n = input$bar_top_n)
    
    cat_levels <- rev(unique(as.character(df_out[[input$col_category]])))
    df_out %>%
      mutate(.category = factor(as.character(.data[[input$col_category]]),
                                levels = cat_levels))
  })
  
  output$barplot_status <- renderText({
    req(ipa_funct(), input$col_pvalue, pvalue_transform(), input$barplot_type)
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter  : ", input$filter_col, " = ",
             paste(input$filter_vals, collapse = ", "), "\n")
    else ""
    base <- paste0(
      filter_line,
      "Plot type        : ", if (input$barplot_type == "bubble") "Bubble plot" else "Bar plot", "\n",
      "P-value column   : ", input$col_pvalue,                                                  "\n",
      "Values detected  : ", if (pvalue_transform() == "raw")
        "raw p-values (0-1) -> -log10 applied" else "already -log10 transformed -> used as-is", "\n",
      "Total pathways   : ", nrow(ipa_funct()),                                                  "\n",
      "Showing top      : ", min(input$bar_top_n, nrow(ipa_funct()))
    )
    if (input$barplot_type == "bubble" &&
        !is.null(input$bubble_size_col) && !is.null(input$bubble_color_col)) {
      base <- paste0(base,
                     "\nBubble size col  : ", input$bubble_size_col,
                     "\nBubble colour col: ", input$bubble_color_col)
    }
    base
  })
  
  make_barplot <- function() {
    req(barplot_data(), input$barplot_type)
    df      <- barplot_data()
    x_label <- if (pvalue_transform() == "raw")
      paste0("-log10(", input$col_pvalue, ")") else input$col_pvalue
    
    if (input$barplot_type == "bar") {
      bar_color <- if (is.null(input$bar_color) || input$bar_color == "") "steelblue" else input$bar_color
      ggplot(df, aes(x = .neg_log10, y = .category)) +
        geom_col(fill = bar_color, width = 0.7) +
        labs(title = paste0("Top ", nrow(df), " Enriched Pathways / Functions"),
             x = x_label, y = NULL) +
        theme_bw(base_size = 12) +
        theme(axis.text.y         = element_text(size = 10),
              axis.text.x         = element_text(size = 10),
              plot.title          = element_text(face = "bold", size = 13),
              panel.grid.major.y  = element_blank())
      
    } else {
      req(input$bubble_size_col, input$bubble_color_col)
      size_raw  <- suppressWarnings(as.numeric(ipa_funct()[[input$bubble_size_col]]))
      color_raw <- suppressWarnings(as.numeric(ipa_funct()[[input$bubble_color_col]]))
      df$.bub_size  <- size_raw [match(as.character(df[[input$col_category]]),
                                       as.character(ipa_funct()[[input$col_category]]))]
      df$.bub_color <- color_raw[match(as.character(df[[input$col_category]]),
                                       as.character(ipa_funct()[[input$col_category]]))]
      color_needs_transform <- all(df$.bub_color >= 0 & df$.bub_color <= 1, na.rm = TRUE)
      if (color_needs_transform) df$.bub_color <- -log10(df$.bub_color)
      color_label <- if (color_needs_transform)
        paste0("-log10(", input$bubble_color_col, ")") else input$bubble_color_col
      palette     <- if (is.null(input$bubble_palette) || input$bubble_palette == "") "Blues" else input$bubble_palette
      use_viridis <- palette %in% c("viridis", "magma", "plasma")
      p <- ggplot(df, aes(x = .neg_log10, y = .category,
                          size = .bub_size, fill = .bub_color)) +
        geom_point(shape = 21, alpha = 0.85, colour = "grey30") +
        scale_size(name  = input$bubble_size_col, range = input$bubble_size_range) +
        labs(title = paste0("Top ", nrow(df), " Enriched Pathways / Functions"),
             x = x_label, y = NULL, fill = color_label) +
        theme_bw(base_size = 12) +
        theme(axis.text.y        = element_text(size = 10),
              axis.text.x        = element_text(size = 10),
              plot.title         = element_text(face = "bold", size = 13),
              panel.grid.major.y = element_blank(),
              legend.position    = "right")
      if (use_viridis) {
        p <- p + scale_fill_viridis_c(option = palette, name = color_label)
      } else {
        p <- p + scale_fill_distiller(palette = palette, direction = 1, name = color_label)
      }
      p
    }
  }
  
  output$barplot <- renderPlot({ req(barplot_data()); make_barplot() },
                               height = function() input$bar_plot_height)
  
  dl_barplot_handler <- function(file) {
    png(file, width = 4200, height = 2400, res = 300, type = "cairo-png")
    print(make_barplot())
    dev.off()
  }
  dl_barplot_filename <- function() {
    paste0(if (input$barplot_type == "bubble") "Bubbleplot_" else "Barplot_",
           "Enriched_Pathways_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  }
  output$download_barplot      <- downloadHandler(filename = dl_barplot_filename, content = dl_barplot_handler)
  output$download_barplot_main <- downloadHandler(filename = dl_barplot_filename, content = dl_barplot_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # CHORD DIAGRAM
  # ══════════════════════════════════════════════════════════════════════════
  chord_data <- reactive({
    req(ipa_funct(), input$col_category, input$col_molecules,
        input$col_pvalue, pvalue_transform(), input$mol_separator)
    df       <- ipa_funct()
    raw_vals <- suppressWarnings(as.numeric(df[[input$col_pvalue]]))
    df <- df %>%
      mutate(.neg_log10 = if (pvalue_transform() == "raw") -log10(raw_vals) else raw_vals) %>%
      filter(!is.na(.neg_log10), is.finite(.neg_log10)) %>%
      arrange(desc(.neg_log10)) %>%
      slice_head(n = input$chord_top_n)
    edges <- df %>%
      select(pathway = all_of(input$col_category), molecules = all_of(input$col_molecules)) %>%
      mutate(molecules = strsplit(as.character(molecules), input$mol_separator, fixed = TRUE)) %>%
      unnest(molecules) %>%
      mutate(molecules = trimws(molecules)) %>%
      filter(molecules != "", !is.na(molecules))
    mol_freq <- edges %>% count(molecules, name = "freq") %>% filter(freq >= input$chord_min_freq)
    edges %>% filter(molecules %in% mol_freq$molecules)
  })
  
  output$chord_status <- renderText({
    req(chord_data())
    cd <- chord_data()
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter  : ", input$filter_col, " = ",
             paste(input$filter_vals, collapse = ", "), "\n")
    else ""
    paste0(
      filter_line,
      "Pathways shown   : ", n_distinct(cd$pathway),  "\n",
      "Molecules shown  : ", n_distinct(cd$molecules), "\n",
      "Total connections: ", nrow(cd),                 "\n",
      "Min appearances  : ", input$chord_min_freq,     "\n",
      "Inner radius     : ", input$chord_inner_radius
    )
  })
  
  make_chord <- function() {
    cd <- chord_data()
    validate(
      need(nrow(cd) > 0, "No connections to display. Try lowering 'Minimum appearances'."),
      need(n_distinct(cd$pathway) >= 2, "At least 2 pathways needed for a chord diagram.")
    )
    mat    <- as.matrix(table(cd$pathway, cd$molecules))
    n_path <- nrow(mat); n_mol <- ncol(mat)
    pathway_color  <- if (is.null(input$chord_pathway_color)  || input$chord_pathway_color  == "") "tomato"    else input$chord_pathway_color
    molecule_color <- if (is.null(input$chord_molecule_color) || input$chord_molecule_color == "") "steelblue" else input$chord_molecule_color
    path_cols   <- colorRampPalette(c(pathway_color,  "lightyellow"))(n_path)
    mol_cols    <- colorRampPalette(c(molecule_color, "lightyellow"))(n_mol)
    grid_colors <- c(setNames(path_cols, rownames(mat)), setNames(mol_cols, colnames(mat)))
    circos.clear()
    par(mar = c(2, 2, 3, 2))
    circos.par("gap.after" = 2, "start.degree" = 90)
    chordDiagram(mat, h = input$chord_inner_radius, grid.col = grid_colors,
                 transparency = 0.4, annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.08), directional = FALSE)
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      sector_name <- get.cell.meta.data("sector.index")
      xlim        <- get.cell.meta.data("xlim")
      ylim        <- get.cell.meta.data("ylim")
      is_pathway  <- sector_name %in% rownames(mat)
      circos.text(mean(xlim), ylim[1] + 0.1, sector_name,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                  cex = if (is_pathway) input$chord_label_size_path else input$chord_label_size_mol)
    }, bg.border = NA)
    title(paste0("Chord Diagram \u2014 Top ", n_distinct(cd$pathway), " Pathways & Molecules"), cex.main = 1.1)
    legend("bottomleft", legend = c("Pathways / Functions", "Molecules"),
           fill = c(pathway_color, molecule_color), border = NA, bty = "n", cex = 0.85)
    circos.clear()
  }
  
  output$chord_plot <- renderPlot({ req(chord_data()); make_chord() },
                                  height = function() input$chord_plot_height)
  dl_chord_handler  <- function(file) { png(file, width = 10, height = 10, units = "in", res = 300); make_chord(); dev.off() }
  output$download_chord      <- downloadHandler(filename = function() "Chord_Diagram_Pathways_Molecules.png", content = dl_chord_handler)
  output$download_chord_main <- downloadHandler(filename = function() "Chord_Diagram_Pathways_Molecules.png", content = dl_chord_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # GENE FREQUENCY BAR PLOT (All pathways OR Top N pathways)
  # ══════════════════════════════════════════════════════════════════════════
  gene_freq_data <- reactive({
    req(ipa_funct(), input$col_category, input$col_molecules,
        input$mol_separator, input$freq_pathway_scope)
    
    df <- ipa_funct()
    
    # Option: restrict to Top N pathways by significance
    if (input$freq_pathway_scope == "top") {
      req(input$col_pvalue, pvalue_transform(), input$freq_top_n_pathways)
      
      raw_vals <- suppressWarnings(as.numeric(df[[input$col_pvalue]]))
      df <- df %>%
        mutate(.neg_log10 = if (pvalue_transform() == "raw") -log10(raw_vals) else raw_vals) %>%
        filter(!is.na(.neg_log10), is.finite(.neg_log10)) %>%
        arrange(desc(.neg_log10)) %>%
        slice_head(n = input$freq_top_n_pathways)
    }
    
    # Expand pathway→gene membership
    edges <- df %>%
      select(pathway = all_of(input$col_category), molecules = all_of(input$col_molecules)) %>%
      mutate(molecules = strsplit(as.character(molecules), input$mol_separator, fixed = TRUE)) %>%
      unnest(molecules) %>%
      mutate(molecules = trimws(molecules)) %>%
      filter(molecules != "", !is.na(molecules)) %>%
      distinct(pathway, molecules)
    
    # Count distinct pathways per molecule
    edges %>%
      count(molecules, name = "pathway_count") %>%
      arrange(desc(pathway_count), molecules)
  })
  
  output$gene_freq_status <- renderText({
    req(gene_freq_data(), input$freq_min_pathways, input$freq_top_n_genes)
    
    gf <- gene_freq_data()
    gf_filt <- gf %>% filter(pathway_count >= input$freq_min_pathways)
    
    scope_line <- if (input$freq_pathway_scope == "all") {
      paste0("Pathways used   : ALL (", n_distinct(ipa_funct()[[input$col_category]]), ")")
    } else {
      paste0("Pathways used   : TOP ", input$freq_top_n_pathways,
             " (ranked by ", input$col_pvalue, ")")
    }
    
    paste0(
      scope_line, "\n",
      "Genes total     : ", nrow(gf), "\n",
      "Genes with >= ", input$freq_min_pathways, " pathways: ", nrow(gf_filt), "\n",
      "Showing top     : ", min(input$freq_top_n_genes, nrow(gf_filt))
    )
  })
  
  make_gene_freq_plot <- function() {
    req(gene_freq_data(), input$freq_min_pathways, input$freq_top_n_genes)
    
    gf <- gene_freq_data() %>%
      filter(pathway_count >= input$freq_min_pathways) %>%
      slice_head(n = input$freq_top_n_genes)
    
    validate(need(nrow(gf) > 0,
                  "No genes meet the frequency cutoff. Lower the cutoff or expand pathway scope."))
    
    ggplot(gf, aes(x = reorder(molecules, -pathway_count), y = pathway_count)) +
      geom_col(fill = "steelblue", color = "white") +
      labs(
        title = paste0("Genes appearing in \u2265 ", input$freq_min_pathways, " pathways"),
        x = "Gene",
        y = "Number of pathways"
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  }
  
  output$gene_freq_plot <- renderPlot({
    make_gene_freq_plot()
  }, height = function() input$freq_plot_height)
  
  dl_gene_freq_filename <- function() {
    paste0("Gene_Frequency_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  }
  dl_gene_freq_handler <- function(file) {
    ggsave(file, plot = make_gene_freq_plot(),
           width = 12, height = 7, dpi = 300, units = "in")
  }
  
  output$download_gene_freq      <- downloadHandler(filename = dl_gene_freq_filename, content = dl_gene_freq_handler)
  output$download_gene_freq_main <- downloadHandler(filename = dl_gene_freq_filename, content = dl_gene_freq_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # HEATMAP  — now uses filtered annotation + filtered sample columns
  # ══════════════════════════════════════════════════════════════════════════
  make_heatmap <- function(gene_filter, title_text) {
    req(input$col_gene, input$col_sample_id, input$col_group)
    data <- norm_prot()
    ann  <- annotation()                       # ← filtered
    
    validate(need(nrow(ann) > 0,
                  "No samples remain after filtering. Please broaden your metadata filter selection."))
    
    ann_ordered <- ann[order(as.character(ann[[input$col_group]])), ]
    ann_ids_ord <- as.character(ann_ordered[[input$col_sample_id]])
    data_cols   <- colnames(data)
    sample_cols <- data_cols[data_cols %in% ann_ids_ord]
    sample_cols <- sample_cols[order(match(sample_cols, ann_ids_ord))]
    
    validate(need(length(sample_cols) > 0,
                  "No sample columns in the normalized data matched the filtered annotation IDs."))
    
    if (!is.null(gene_filter)) {
      data <- data[as.character(data[[input$col_gene]]) %in% as.character(gene_filter), ]
    }
    
    validate(
      need(nrow(data) > 1,
           paste0("Fewer than 2 genes matched in the normalized data.\n",
                  "Check that the correct Gene/Protein ID column is selected.\n",
                  "Genes searched: ", paste(head(gene_filter, 5), collapse = ", "), "...\n",
                  "Gene column used: '", input$col_gene, "'"))
    )
    
    mat           <- as.matrix(data[, sample_cols, drop = FALSE])
    rownames(mat) <- as.character(data[[input$col_gene]])
    
    row_vars <- apply(mat, 1, var, na.rm = TRUE)
    mat      <- mat[!is.na(row_vars) & row_vars > 0, , drop = FALSE]
    
    validate(need(nrow(mat) > 1,
                  "All matched genes have zero variance across samples \u2014 cannot draw heatmap."))
    
    mat_scaled <- t(scale(t(mat)))
    col_annotation           <- data.frame(as.character(ann_ordered[[input$col_group]]),
                                           row.names   = as.character(ann_ordered[[input$col_sample_id]]),
                                           check.names = FALSE)
    colnames(col_annotation) <- input$col_group
    col_annotation           <- col_annotation[colnames(mat), , drop = FALSE]
    groups   <- unique(as.character(ann[[input$col_group]]))
    grp_cols <- setNames(
      sapply(seq_along(groups), function(i) {
        val <- input[[paste0("color_group_", i)]]
        if (is.null(val) || val == "") "grey" else val
      }),
      groups
    )
    ann_colors        <- list(grp_cols)
    names(ann_colors) <- input$col_group
    
    # Build a subtitle that shows the active filter
    filter_subtitle <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("\nFilter: ", input$filter_col, " = ", paste(input$filter_vals, collapse = ", "))
    else ""
    
    pheatmap(
      mat_scaled,
      annotation_col    = col_annotation,
      annotation_colors = ann_colors,
      cluster_rows      = TRUE, cluster_cols = FALSE,
      show_rownames     = nrow(mat) <= 100, show_colnames = TRUE,
      scale             = "none",
      color             = colorRampPalette(c("darkblue", "white", "darkred"))(100),
      border_color      = NA, fontsize_row = 8, fontsize_col = 10,
      main              = paste0(title_text, filter_subtitle)
    )
  }
  
  output$heatmap_plot <- renderPlot({
    req(norm_prot(), annotation(), genes_interest(), input$selected_category)
    make_heatmap(gene_filter = genes_interest(),
                 title_text  = paste0(input$selected_category, "\n(Z-score, by Group)"))
  }, height = function() input$plot_height)
  
  dl_single_filename <- function()
    paste0("Heatmap_", gsub("[^A-Za-z0-9]", "_", input$selected_category), ".png")
  dl_single_handler  <- function(file) {
    png(file, width = 8, height = 10, units = "in", res = 300)
    make_heatmap(gene_filter = genes_interest(),
                 title_text  = paste0(input$selected_category, "\n(Z-score, by Group)"))
    dev.off()
  }
  output$download_single      <- downloadHandler(filename = dl_single_filename, content = dl_single_handler)
  output$download_single_main <- downloadHandler(filename = dl_single_filename, content = dl_single_handler)
  
  download_all_handler <- function(file) {
    req(input$col_category, input$col_molecules, input$col_gene, input$mol_separator)
    tmp_dir <- tempdir()
    saved   <- c()
    for (cat_val in unique(ipa_funct()[[input$col_category]])) {
      genes <- ipa_funct() %>%
        filter(.data[[input$col_category]] == cat_val) %>%
        pull(.data[[input$col_molecules]]) %>%
        strsplit(input$mol_separator, fixed = TRUE) %>%
        unlist() %>% trimws() %>% unique()
      if (length(genes) == 0) next
      data <- norm_prot()[as.character(norm_prot()[[input$col_gene]]) %in% genes, ]
      if (nrow(data) <= 1) next
      clean    <- gsub("[^A-Za-z0-9]", "_", cat_val)
      png_path <- file.path(tmp_dir, paste0("Heatmap_", clean, ".png"))
      png(png_path, width = 8, height = 10, units = "in", res = 300)
      make_heatmap(gene_filter = genes, title_text = paste0(cat_val, "\n(Z-score, by Group)"))
      dev.off()
      saved <- c(saved, png_path)
    }
    zip(zipfile = file, files = saved, flags = "-j")
  }
  output$download_all      <- downloadHandler(filename = function() "All_Heatmaps.zip", content = download_all_handler)
  output$download_all_main <- downloadHandler(filename = function() "All_Heatmaps.zip", content = download_all_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # BOXPLOT / VIOLIN PLOT  — uses filtered annotation + filtered samples
  # ══════════════════════════════════════════════════════════════════════════
  all_molecules <- reactive({
    req(ipa_funct(), input$col_molecules, input$mol_separator)
    ipa_funct()[[input$col_molecules]] %>%
      strsplit(input$mol_separator, fixed = TRUE) %>%
      unlist() %>% trimws() %>% unique() %>% sort()
  })
  
  output$ui_boxplot_protein <- renderUI({
    req(all_molecules(), norm_prot(), input$col_gene)
    gene_col     <- as.character(norm_prot()[[input$col_gene]])
    mols_present <- all_molecules()[all_molecules() %in% gene_col]
    validate(need(length(mols_present) > 0,
                  "No molecules from the pathway file found in the normalized data."))
    selectizeInput(
      inputId  = "boxplot_protein",
      label    = paste0("Search protein / gene (", length(mols_present), " matched)"),
      choices  = mols_present,
      selected = mols_present[1],
      options  = list(placeholder = "Type to search...", maxOptions = 2000)
    )
  })
  
  boxplot_data <- reactive({
    req(norm_prot(), annotation(), input$col_gene, input$col_sample_id,
        input$col_group, input$boxplot_protein, sample_cols_resolved())
    protein     <- input$boxplot_protein
    data        <- norm_prot()
    ann         <- annotation()                 # ← filtered
    sample_cols <- sample_cols_resolved()       # ← already filtered via annotation()
    
    validate(need(nrow(ann) > 0,
                  "No samples remain after filtering. Please broaden your metadata filter selection."))
    validate(need(length(sample_cols) > 0,
                  "No sample columns matched the filtered annotation IDs."))
    
    gene_col    <- as.character(data[[input$col_gene]])
    row_idx     <- which(gene_col == as.character(protein))
    validate(need(length(row_idx) >= 1,
                  paste0("Protein '", protein, "' not found.")))
    prot_row   <- data[row_idx[1], , drop = FALSE]
    vals       <- as.numeric(unlist(prot_row[, sample_cols, drop = FALSE]))
    ann_ids    <- as.character(ann[[input$col_sample_id]])
    ann_groups <- as.character(ann[[input$col_group]])
    ann_pos    <- match(sample_cols, ann_ids)
    validate(need(!all(is.na(ann_pos)),
                  "None of the sample columns matched the annotation Sample ID column."))
    groups <- ann_groups[ann_pos]
    df <- data.frame(sample = sample_cols, value = vals, group = groups,
                     stringsAsFactors = FALSE, row.names = NULL)
    df <- df[!is.na(df$value) & !is.na(df$group) & df$group != "", ]
    validate(need(nrow(df) > 0, "No valid rows after group assignment."))
    df
  })
  
  boxplot_group_colors <- reactive({
    req(annotation(), input$col_group)
    groups <- unique(as.character(annotation()[[input$col_group]]))
    setNames(
      sapply(seq_along(groups), function(i) {
        val <- input[[paste0("color_group_", i)]]
        if (is.null(val) || val == "") "grey" else val
      }),
      groups
    )
  })
  
  output$boxplot_status <- renderText({
    req(boxplot_data(), input$boxplot_protein, input$plot_type)
    df <- boxplot_data()
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter  : ", input$filter_col, " = ",
             paste(input$filter_vals, collapse = ", "), "\n")
    else ""
    paste0(
      filter_line,
      "Plot type        : ", if (input$plot_type == "violin") "Violin plot" else "Boxplot", "\n",
      "Selected protein : ", input$boxplot_protein,                    "\n",
      "Samples found    : ", nrow(df),                                 "\n",
      "Groups           : ", paste(unique(df$group), collapse = ", "), "\n",
      "Values range     : ", round(min(df$value, na.rm = TRUE), 2),
      " \u2014 ", round(max(df$value, na.rm = TRUE), 2)
    )
  })
  
  make_boxplot <- function() {
    req(boxplot_data(), boxplot_group_colors(), input$plot_type)
    df         <- boxplot_data()
    grp_colors <- boxplot_group_colors()
    protein    <- input$boxplot_protein
    ann_group_order <- unique(as.character(annotation()[[input$col_group]]))
    present_groups  <- intersect(ann_group_order, unique(df$group))
    grp_colors      <- grp_colors[names(grp_colors) %in% present_groups]
    df$group        <- factor(df$group, levels = present_groups)
    
    filter_subtitle <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("\n(", input$filter_col, " = ", paste(input$filter_vals, collapse = ", "), ")")
    else ""
    
    p <- ggplot(df, aes(x = group, y = value, fill = group, colour = group))
    if (input$plot_type == "violin") {
      p <- p +
        geom_violin(alpha = 0.5, trim = FALSE, width = 0.8) +
        geom_jitter(width = 0.12, size = 2.5, alpha = 0.85)
    } else {
      p <- p +
        geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
        geom_jitter(width = 0.15, size = 2.5, alpha = 0.85)
    }
    p +
      scale_fill_manual(values   = grp_colors) +
      scale_colour_manual(values = grp_colors) +
      labs(title    = paste0("Normalized Abundance / counts \u2014 ", protein, filter_subtitle),
           x = NULL, y = "Normalized abundance / counts") +
      theme_bw(base_size = 13) +
      theme(legend.position    = "none",
            plot.title         = element_text(face = "bold", size = 14),
            axis.text.x        = element_text(size = 12),
            axis.text.y        = element_text(size = 11),
            panel.grid.major.x = element_blank())
  }
  
  output$boxplot_plot <- renderPlot({ req(boxplot_data()); make_boxplot() },
                                    height = function() input$boxplot_height)
  
  dl_boxplot_filename <- function()
    paste0(if (input$plot_type == "violin") "Violin_" else "Boxplot_",
           gsub("[^A-Za-z0-9]", "_", input$boxplot_protein), ".png")
  dl_boxplot_handler  <- function(file)
    ggsave(file, plot = make_boxplot(), width = 8, height = 6, dpi = 300, units = "in")
  
  output$download_boxplot_sidebar <- downloadHandler(filename = dl_boxplot_filename, content = dl_boxplot_handler)
  output$download_boxplot_main    <- downloadHandler(filename = dl_boxplot_filename, content = dl_boxplot_handler)
}

# ══════════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)