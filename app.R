# app.R
# Heatmaps, Bar/Bubble Plots, Chord Diagrams, Boxplot/Violin, Volcano and PCA plots visualization
# v.1.4.6  — American spelling throughout the interface (color, centered, gray, ...)
# v.1.4.5  — bubble colors: diverging scale for signed columns (NES, z-score), numeric columns only, -log10 only for p-value columns
# v.1.4.4  — chord labels can wrap; optional prefix removal (e.g. GOBP_) for displayed pathway names
# v.1.4.3  — download width/height per plot and text-size scaling (ragg); "Molecules" wording in heatmap status
# v.1.4.2  — Molecule Frequency (renamed) colored by log2FC; chord label/color fixes; box/violin x label, trimmed violins
# v.1.4.1  — chord diagram no longer fails with many molecules (adaptive sector gaps)
# v.1.4.0  — adds optional differential analysis file: volcano plot, PCA, DE heatmap
# Copyright 2026 RGM
# MIT License — see README.md for full license text
# ─────────────────────────────────────────────────────────

library(shiny)
library(pheatmap)
library(tidyverse)
library(circlize)
library(ggrepel)
library(ragg)      # PNG downloads with adjustable text scaling


options(shiny.maxRequestSize = 500 * 1024^2)

# Default download sizes (inches) for each plot; users can change them per tab
dl_defaults <- list(bar = c(14, 8), chord = c(10, 10), freq = c(12, 7), heat = c(8, 10),
                    box = c(8, 6), volcano = c(9, 7), pca = c(9, 7), deheat = c(8, 10))

# Width / height inputs shown next to each tab's download button
dl_size_ui <- function(key) {
  d <- dl_defaults[[key]]
  fluidRow(
    column(3, numericInput(paste0("dlw_", key), "Download width (in)",  value = d[1], min = 1, max = 30, step = 0.25)),
    column(3, numericInput(paste0("dlh_", key), "Download height (in)", value = d[2], min = 1, max = 30, step = 0.25))
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- fluidPage(
  
  tags$head(tags$title("EnrichViz v1.4.6")),
  
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
    #keyword_filter_panel {
      background: #fff8e1;
      border: 1px solid #ffe082;
      border-radius: 6px;
      padding: 10px 12px 6px 12px;
      margin-bottom: 6px;
    }
    #keyword_filter_panel h4 { margin-top: 4px; color: #7a5800; }
    #de_panel {
      background: #eef8f0;
      border: 1px solid #b5dcbd;
      border-radius: 6px;
      padding: 10px 12px 6px 12px;
      margin-bottom: 6px;
    }
    #de_panel h4 { margin-top: 4px; color: #1e6b33; }
  ")),
  
  tags$div(
    style = paste0(
      "position:fixed; top:0; left:0; width:100%; background:white;",
      "z-index:200; border-bottom:1px solid #ddd;",
      "padding:8px 20px 6px 20px; margin:0;"
    ),
    
    # Row 1 – logo + title side by side
    tags$div(
      style = "display:flex; align-items:center; gap:14px;",
      tags$img(
        src    = "enrichviz.jpg",
        height = "64px",
        style  = "border-radius:6px; flex-shrink:0;"
      ),
      tags$div(
        tags$h2(
          "Heatmap, Bar/Bubble, Chord, Molecule Frequency, Box/Violin, Volcano, and PCA Plots of Pathway and Differential Analysis Results",
          style = "margin:0 0 2px 0; font-size:18px; font-weight:bold;"
        ),
        tags$p(
          style = "margin:0; color:#444; font-size:13px; line-height:1.35;",
          HTML(
            "This app integrates normalized expression/abundance data with sample metadata
             and pathway/function enrichment results to generate bar/bubble plots, chord
             diagrams, heatmaps, molecule frequency, and box/violin plots.
             <b>Files needed (CSV):</b> (1) Normalized data (rows = genes/proteins; columns = samples),
             (2) metadata/annotation (sample IDs + group/condition),
             (3) enrichment results (pathway/function name, p-value/FDR or -log10 value,
             and member molecules/genes), and optionally (4) differential analysis results
             (gene ID, log2 fold change, adjusted p-value) for volcano, PCA and DE heatmap."
          )
        )
      )
    ),
    
    # Row 2 – citation
    tags$div(
      style = paste0(
        "margin-top:6px; padding:5px 10px; border-radius:4px;",
        "background:#f0f7ff; border:1px solid #b8d4f0;",
        "font-size:12px; color:#1a4f8a; line-height:1.4;"
      ),
      HTML(
        "<b>Citation:</b> Garcia-Milian R. EnrichViz: An Interactive R Shiny Application
         for Visualization of Pathway Enrichment Results from Omics Data.
         <em>bioRxiv</em> [Preprint]. 2026 Jun 24:2026.06.19.733398.
         doi: <a href='https://doi.org/10.64898/2026.06.19.733398'
                 target='_blank'>10.64898/2026.06.19.733398</a>.
         PMID:&nbsp;42395511; PMCID:&nbsp;PMC13320841."
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Upload Data Files"),
      fileInput("file_norm_prot",  "Normalized Data (.csv)",               accept = ".csv"),
      fileInput("file_annotation", "Metadata or Sample Annotation (.csv)", accept = ".csv"),
      fileInput("file_ipa_funct",  "Enriched Pathway or Functions (.csv)", accept = ".csv"),
      fileInput("file_de",         "Differential Analysis Results (.csv) \u2014 optional", accept = ".csv"),
      numericInput("dl_text_scale", "Text size in downloaded plots (×)", value = 1, min = 0.3, max = 3, step = 0.1),
      helpText("Downloads are 300 dpi PNG files of the width and height set below each plot.",
               "For a journal panel, set them to the final printed size (e.g. 3.5 in for one column,",
               "7.5 in for full width); text then prints at its real point size. Use this factor to",
               "enlarge (> 1) or shrink (< 1) all text and lines in the downloads."),
      
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
      helpText("Columns and separator are detected automatically from names and contents;",
               "change them here if needed."),
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
      div(
        id = "keyword_filter_panel",
        h4("\U0001f50e Keyword Filter \u2014 Pathways"),
        helpText(tags$b("Filter by pathway name:"),
                 " comma-separated keywords (partial text match), OR logic, case-insensitive."),
        textInput("pathway_keyword_filter",
                  label       = NULL,
                  placeholder = "e.g. immune, signaling, metabolism"),
        helpText(tags$b("Filter by molecule / gene symbol:"),
                 " comma-separated symbols, OR logic, case-insensitive. Exact symbol match;",
                 " add * as a wildcard (e.g. SYN* matches SYN1, SYN2, SYN3)."),
        textInput("molecule_keyword_filter",
                  label       = NULL,
                  placeholder = "e.g. MAPK1, TP53, SYN*"),
        uiOutput("ui_keyword_filter_summary")
      ),
      hr(),
      
      radioButtons("barplot_type",
                   label    = "Plot type",
                   choices  = c("Bar plot" = "bar", "Bubble plot" = "bubble"),
                   selected = "bar",
                   inline   = TRUE),
      
      uiOutput("ui_col_pvalue"),
      uiOutput("ui_pvalue_transform_info"),
      numericInput("bar_top_n", "Show top N pathways", value = 20, min = 5, max = 200),
      fluidRow(
        column(7, selectInput("label_mode", "Long pathway names",
                              choices  = c("Wrap" = "wrap", "Truncate" = "truncate", "Show in full" = "full"),
                              selected = "wrap")),
        column(5, numericInput("label_width", "Max characters", value = 50, min = 15, max = 200, step = 5))
      ),
      helpText("Wrap splits long names over up to 3 lines (at spaces, '_' or '-'); Truncate cuts them",
               "with '\u2026'. Only the axis labels change, not the data."),
      textInput("strip_prefix", "Remove prefix from pathway names (display only)", value = "",
                placeholder = "e.g. GOBP_, KEGG_, REACTOME_"),
      helpText("Comma-separated prefixes removed from the start of pathway names in the bar / bubble,",
               "chord and heatmap labels. Filters and downloads of tables are not affected."),
      
      conditionalPanel(
        condition = "input.barplot_type == 'bar'",
        uiOutput("ui_bar_fill_col"),
        conditionalPanel(
          condition = "input.bar_fill_col == ''",
          textInput("bar_color", "Bar fill color (hex or name)", value = "steelblue")
        ),
        conditionalPanel(
          condition = "input.bar_fill_col != ''",
          radioButtons("bar_fill_scale", "Color scale",
                       choices  = c("Auto" = "auto", "Gradient" = "gradient", "Diverging (centered at 0)" = "diverging"),
                       selected = "auto", inline = TRUE),
          helpText("Auto: diverging (white = 0) when the column has negative values (e.g. z-score, NES),",
                   "otherwise a light-to-dark gradient. P-value columns are shown as -log10.",
                   "Pathways with no value are gray."),
          conditionalPanel(
            condition = "input.bar_fill_scale != 'diverging'",
            selectInput("bar_palette", "Gradient palette",
                        choices = c("Blues", "Reds", "Purples", "Greens",
                                    "OrRd", "YlOrRd", "viridis", "magma", "plasma"),
                        selected = "Blues")
          ),
          conditionalPanel(
            condition = "input.bar_fill_scale != 'gradient'",
            fluidRow(
              column(6, textInput("bar_neg_color", "Negative color", value = "steelblue")),
              column(6, textInput("bar_pos_color", "Positive color", value = "darkorange"))
            )
          )
        )
      ),
      
      conditionalPanel(
        condition = "input.barplot_type == 'bubble'",
        uiOutput("ui_bubble_size_col"),
        uiOutput("ui_bubble_color_col"),
        sliderInput("bubble_size_range",
                    label = "Bubble size range (min / max)",
                    min = 1, max = 30, value = c(3, 15), step = 1),
        radioButtons("bubble_fill_scale", "Color scale",
                     choices  = c("Auto" = "auto", "Gradient" = "gradient", "Diverging (centered at 0)" = "diverging"),
                     selected = "auto", inline = TRUE),
        helpText("Auto: diverging (white = 0) when the color column has negative values (e.g. NES,",
                 "z-score), otherwise a gradient. P-value columns are shown as -log10."),
        conditionalPanel(
          condition = "input.bubble_fill_scale != 'diverging'",
          selectInput("bubble_palette",
                      label   = "Gradient palette",
                      choices = c("Blues", "Reds", "Purples", "Greens",
                                  "OrRd", "YlOrRd", "RdYlBu", "viridis",
                                  "magma", "plasma"),
                      selected = "Blues")
        ),
        conditionalPanel(
          condition = "input.bubble_fill_scale != 'gradient'",
          fluidRow(
            column(6, textInput("bubble_neg_color", "Negative color", value = "steelblue")),
            column(6, textInput("bubble_pos_color", "Positive color", value = "darkorange"))
          )
        )
      ),
      
      numericInput("bar_plot_height", "Plot height (px)", value = 600, min = 300, max = 2000),
      downloadButton("download_barplot", "Download Plot (.png)"),
      
      hr(),
      
      h4("Chord Diagram Settings"),
      numericInput("chord_top_n",           "Top N pathways to include",            value = 10,   min = 2,   max = 50),
      numericInput("chord_min_freq",        "Minimum times a molecule must appear", value = 2,    min = 1,   max = 20),
      helpText("Increase 'minimum appearances' to reduce clutter."),
      sliderInput("chord_inner_radius",     "Link curve depth",                     min = 0.2, max = 0.9, value = 0.4, step = 0.05),
      helpText("How far the links bend toward the center. Lower = links stay near the outer",
               "ring (can look crowded); higher = links arc through the center and are easier to follow.",
               "This changes appearance only; it does not filter pathways or molecules."),
      numericInput("chord_label_size_path", "Pathway label font size",              value = 0.55, min = 0.3, max = 1.5, step = 0.05),
      numericInput("chord_label_size_mol",  "Molecule label font size",             value = 0.70, min = 0.3, max = 1.5, step = 0.05),
      radioButtons("chord_label_mode", "Long pathway labels",
                   choices = c("Truncate" = "truncate", "Wrap (up to 3 lines)" = "wrap"),
                   selected = "truncate", inline = TRUE),
      numericInput("chord_label_width",     "Max characters per pathway label line", value = 35, min = 10, max = 150, step = 5),
      sliderInput("chord_canvas", "Space around the circle for labels", min = 1.2, max = 2.5, value = 1.5, step = 0.05),
      checkboxInput("chord_show_title", "Show plot title", value = TRUE),
      helpText("Truncate cuts long pathway names with '\u2026'; Wrap splits them over up to 3 lines",
               "(at '_', '-' or spaces). Wrapped labels need wide sectors; raise the plot size if they overlap."),
      textInput("chord_pathway_color",      "Pathway sector color (hex or name)",  value = "tomato"),
      textInput("chord_molecule_color",     "Molecule sector color (hex or name)", value = "steelblue"),
      numericInput("chord_plot_height",     "Chord plot height (px)",               value = 750,  min = 400, max = 2000),
      downloadButton("download_chord", "Download Chord Diagram (.png)"),
      
      hr(),
      
      h4("Molecule Frequency Plot Settings"),
      radioButtons(
        "freq_pathway_scope",
        "Pathways used to compute molecule frequency",
        choices  = c("All pathways" = "all", "Top N pathways (by significance)" = "top"),
        selected = "all",
        inline   = FALSE
      ),
      checkboxInput("freq_use_keyword_filter",
                    "Apply the pathway keyword filter (if one is active)", value = TRUE),
      conditionalPanel(
        condition = "input.freq_pathway_scope == 'top'",
        numericInput("freq_top_n_pathways",
                     "Top N pathways",
                     value = 20, min = 2, max = 200)
      ),
      numericInput("freq_min_pathways",
                   "Minimum # pathways a molecule must appear in",
                   value = 5, min = 1, max = 100),
      numericInput("freq_top_n_genes",
                   "Show top N molecules (after filtering)",
                   value = 50, min = 5, max = 500),
      numericInput("freq_plot_height",
                   "Frequency plot height (px)",
                   value = 600, min = 300, max = 2000),
      checkboxInput("freq_color_lfc",
                    "Color bars by log2 fold change (needs the differential analysis file)",
                    value = TRUE),
      downloadButton("download_gene_freq", "Download Frequency Plot (.png)"),
      
      hr(),
      
      h4("Heatmap Settings"),
      # Choices are sent from the server on demand (server-side selectize)
      selectizeInput("selected_category", "Select pathway / function", choices = NULL,
                     options = list(placeholder = "Upload the enrichment file, then type to search...")),

      # clustering toggles 
      checkboxInput("heatmap_cluster_rows", "Cluster rows",    value = TRUE),
      checkboxInput("heatmap_cluster_cols", "Cluster columns", value = FALSE),
      numericInput("heatmap_fontsize_row", "Row label font size (pt)", value = 8, min = 2, max = 20, step = 0.5),
      
      # Color scale limit
      sliderInput("heatmap_scale_limit",
                  "Heatmap color scale limit (± Z-score)",
                  min = 1, max = 4, value = 2, step = 0.5),
      helpText("Z-scores beyond \u00b1 this limit get the darkest color. Lower = darker, higher-contrast",
               "heatmap (small differences stand out); higher = lighter heatmap with more shades",
               "between extremes. Also used by the DE Heatmap."),
      
      hr(),
      uiOutput("ui_group_colors"),
      numericInput("plot_height", "Heatmap height (px)", value = 700, min = 300, max = 2000),
      
      h4("Download Heatmaps"),
      downloadButton("download_single", "Download Selected Heatmap (.png)"),
      br(), br(),
      uiOutput("ui_download_all"),

      hr(),
      
      h4("Boxplot / Violin Settings"),
      radioButtons("plot_type",
                   label    = "Plot type",
                   choices  = c("Boxplot" = "box", "Violin plot" = "violin"),
                   selected = "box",
                   inline   = TRUE),
      selectizeInput("boxplot_protein", "Search protein / gene / metabolite", choices = NULL,
                     options = list(placeholder = "Type to search...")),
      uiOutput("ui_boxplot_protein_msg"),
helpText("Group colors are shared with Heatmap Settings above."),
      numericInput("boxplot_height", "Plot height (px)", value = 500, min = 300, max = 2000),
      downloadButton("download_boxplot_sidebar", "Download Plot (.png)"),
      
      hr(),
      
      #  Differential analysis: volcano / PCA / DE heatmap 
      div(
        id = "de_panel",
        h4("\U0001f30b Differential Analysis"),
        helpText("Upload the optional differential analysis file (e.g. DESeq2, limma, MetaboAnalyst results)
                  to enable the Volcano and DE Heatmap tabs. Works for genes, proteins or metabolites. PCA works without it."),
        uiOutput("ui_de_cols"),
        numericInput("de_lfc_cut",  "|log2 fold change| \u2265", value = 1,    min = 0, max = 10, step = 0.25),
        numericInput("de_padj_cut", "Adjusted p-value <",          value = 0.05, min = 0, max = 1,  step = 0.01),
        uiOutput("ui_de_summary"),
        checkboxInput("de_log_transform",
                      "Apply log2(x + 1) to normalized data (PCA & DE heatmap)", value = TRUE),
        helpText("Untick if the normalized data are already log-transformed (common for proteomics and metabolomics).")
      ),
      
      h4("Volcano Plot Settings"),
      numericInput("volcano_top_labels", "Label top N significant features (by padj)", value = 10, min = 0, max = 200),
      textInput("volcano_label_genes", "Also label these features (comma-separated IDs)",
                placeholder = "e.g. NPAS4, P04637, Glutamate"),
      textInput("volcano_up_color",   "Up-regulated color (hex or name)",   value = "firebrick"),
      textInput("volcano_down_color", "Down-regulated color (hex or name)", value = "steelblue"),
      numericInput("volcano_height", "Plot height (px)", value = 600, min = 300, max = 2000),
      downloadButton("download_volcano", "Download Volcano Plot (.png)"),
      
      hr(),
      
      h4("PCA Settings"),
      radioButtons("pca_gene_set",
                   label    = "Features used for PCA",
                   choices  = c("All features"                          = "all",
                                "Top N most variable features"          = "top",
                                "DE-significant features (needs DE file)" = "de"),
                   selected = "top"),
      conditionalPanel(
        condition = "input.pca_gene_set == 'top'",
        numericInput("pca_top_n", "Top N most variable features", value = 500, min = 10, max = 50000)
      ),
      checkboxInput("pca_scale", "Scale features to unit variance", value = FALSE),
      fluidRow(
        column(6, numericInput("pca_x", "X-axis PC", value = 1, min = 1, max = 20)),
        column(6, numericInput("pca_y", "Y-axis PC", value = 2, min = 1, max = 20))
      ),
      checkboxInput("pca_show_labels", "Show sample labels", value = FALSE),
      helpText("Point colors are shared with Heatmap Settings above."),
      numericInput("pca_height", "Plot height (px)", value = 600, min = 300, max = 2000),
      downloadButton("download_pca", "Download PCA Plot (.png)"),
      
      hr(),
      
      h4("DE Heatmap Settings"),
      radioButtons("de_heatmap_direction",
                   label    = "Features to show",
                   choices  = c("Up & down" = "both", "Up only" = "up", "Down only" = "down"),
                   selected = "both", inline = TRUE),
      numericInput("de_heatmap_top_n", "Max features (top N by padj)", value = 50, min = 2, max = 2000),
      helpText("Clustering toggles, color-scale limit and group colors are shared with Heatmap Settings above."),
      numericInput("de_heatmap_height", "Plot height (px)", value = 800, min = 300, max = 3000),
      downloadButton("download_de_heatmap", "Download DE Heatmap (.png)")
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
                 dl_size_ui("bar"),
                 downloadButton("download_barplot_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("Chord Diagram",
                 br(),
                 verbatimTextOutput("chord_status"),
                 plotOutput("chord_plot", height = "750px"),
                 br(),
                 dl_size_ui("chord"),
                 downloadButton("download_chord_main", "\u2b07 Download Chord Diagram (.png)")
        ),
        
        tabPanel("Molecule Frequency",
                 br(),
                 verbatimTextOutput("gene_freq_status"),
                 plotOutput("gene_freq_plot", height = "600px"),
                 br(),
                 dl_size_ui("freq"),
                 downloadButton("download_gene_freq_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("Heatmap",
                 br(),
                 verbatimTextOutput("status"),
                 plotOutput("heatmap_plot", height = "700px"),
                 br(),
                 dl_size_ui("heat"),
                 downloadButton("download_single_main", "\u2b07 Download Selected Heatmap (.png)"),
                 br(), br(),
                 uiOutput("ui_download_all_main")
),
        
        tabPanel("Boxplot",
                 br(),
                 verbatimTextOutput("boxplot_status"),
                 plotOutput("boxplot_plot", height = "500px"),
                 br(),
                 dl_size_ui("box"),
                 downloadButton("download_boxplot_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("Volcano",
                 br(),
                 verbatimTextOutput("volcano_status"),
                 plotOutput("volcano_plot", height = "600px"),
                 br(),
                 dl_size_ui("volcano"),
                 downloadButton("download_volcano_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("PCA",
                 br(),
                 verbatimTextOutput("pca_status"),
                 plotOutput("pca_plot", height = "600px"),
                 br(),
                 dl_size_ui("pca"),
                 downloadButton("download_pca_main", "\u2b07 Download Plot (.png)")
        ),
        
        tabPanel("DE Heatmap",
                 br(),
                 verbatimTextOutput("de_heatmap_status"),
                 plotOutput("de_heatmap_plot", height = "800px"),
                 br(),
                 dl_size_ui("deheat"),
                 downloadButton("download_de_heatmap_main", "\u2b07 Download Plot (.png)")
        )
      )
    )
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {
  
  # ── 0. Helpers for reading and recognising enrichment files ──────────────
  # Enrichment exports differ by tool (IPA, clusterProfiler, g:Profiler, Enrichr,
  # DAVID, fgsea, ...). These helpers find the right columns by name and content.

  # Reads a CSV, skipping any lines above the real header (e.g. IPA's
  # "(c) QIAGEN. All rights reserved." line) and naming unnamed columns.
  read_enrichment_csv <- function(path) {
    lines <- readLines(path, n = 50, warn = FALSE)
    nf    <- suppressWarnings(count.fields(textConnection(lines), sep = ",", quote = "\"",
                                           blank.lines.skip = FALSE, comment.char = ""))
    nf[is.na(nf)] <- 0
    multi <- nf[nf >= 2]
    skip  <- if (length(multi) > 0) {
      typical <- as.integer(names(which.max(table(multi))))   # most common field count
      which(nf == typical)[1] - 1
    } else 0
    df <- read.csv(path, check.names = FALSE, skip = skip)
    empty <- is.na(names(df)) | trimws(names(df)) == ""
    names(df)[empty] <- paste0("Column", which(empty))
    df
  }

  # Mean number of separator characters per value: how "list-like" a column is
  separator_density <- function(x, seps = c(",", "/", ";", "|")) {
    x <- head(as.character(x[!is.na(x) & x != ""]), 200)
    if (length(x) == 0) return(setNames(rep(0, length(seps)), seps))
    sapply(seps, function(sp) mean(lengths(regmatches(x, gregexpr(sp, x, fixed = TRUE)))))
  }
  is_text_col    <- function(x) is.character(x) || is.factor(x)
  is_numeric_col <- function(x) mean(!is.na(suppressWarnings(as.numeric(as.character(x))))) > 0.5

  # First column (in file order) whose lower-cased name contains a keyword;
  # keywords are tried in priority order
  first_keyword_match <- function(cols, keywords) {
    for (k in keywords) {
      hit <- cols[grepl(k, tolower(cols), fixed = TRUE)]
      if (length(hit) > 0) return(hit[1])
    }
    NULL
  }

  # Pathway / term name: a text column that is not a member list
  detect_category_col <- function(df) {
    cols <- colnames(df)
    cand <- cols[sapply(df, is_text_col)]
    cand <- cand[sapply(cand, function(cl) max(separator_density(df[[cl]])) < 2)]
    if (length(cand) == 0) return(cols[1])
    hit <- first_keyword_match(cand, c("description", "term_name", "pathway", "term",
                                       "gene_set", "geneset", "name"))
    if (!is.null(hit)) return(hit)
    cand[which.max(sapply(cand, function(cl) length(unique(df[[cl]]))))]
  }

  # Member list: a text column holding delimited IDs; NULL if the file has none
  detect_molecules_col <- function(df) {
    cols <- colnames(df)
    text <- cols[-1][sapply(df[cols[-1]], is_text_col)]
    if (length(text) == 0) return(NULL)
    hit <- first_keyword_match(text, c("core_enrichment", "geneid", "gene_id", "molecules",
                                       "members", "intersection", "leadingedge", "leading_edge",
                                       "overlapping", "genes", "gene", "mol"))
    # A keyword hit must actually look like a list (leading_edge in clusterProfiler is "tags=..%")
    if (!is.null(hit) && max(separator_density(df[[hit]])) >= 1) return(hit)
    dens <- sapply(text, function(cl) max(separator_density(df[[cl]])))
    if (max(dens) >= 1) text[which.max(dens)] else NULL
  }

  # Significance: adjusted p-values preferred, then raw p-values, then -log10 columns
  detect_pvalue_col <- function(df) {
    cols <- colnames(df)
    num  <- cols[sapply(df, is_numeric_col)]
    hit  <- first_keyword_match(num, c(
      "padj", "p.adjust", "p_adjust", "adjusted p", "adjusted_p", "adj.p", "adj_p", "adj. p",
      "p_adj", "fdr", "benjamini", "qvalue", "q.value", "q-value", "q_value", "bonferroni",
      "pvalue", "p.value", "p-value", "p_value", "pval", "p value", "-log", "log(p", "log10"))
    if (!is.null(hit)) hit else if (length(num) > 0) num[1] else cols[min(2, length(cols))]
  }

  # ── 1. Load uploaded data ─────────────────────────────────────────────────
  norm_prot  <- reactive({ req(input$file_norm_prot);  read.csv(input$file_norm_prot$datapath,  check.names = FALSE) })
  ipa_funct  <- reactive({ req(input$file_ipa_funct);  read_enrichment_csv(input$file_ipa_funct$datapath) })
  annotation_raw <- reactive({ req(input$file_annotation); read.csv(input$file_annotation$datapath, check.names = FALSE) })
  split_keywords <- function(txt) {
    if (is.null(txt)) return(character(0))
    kw <- tolower(trimws(strsplit(trimws(txt), ",", fixed = TRUE)[[1]]))
    kw[nchar(kw) > 0]
  }
  keyword_filter_active <- reactive({
    length(split_keywords(input$pathway_keyword_filter)) > 0 ||
      length(split_keywords(input$molecule_keyword_filter)) > 0
  })

  # Global keyword-filtered pathway table (pathway name AND/OR molecule)
ipa_funct_filtered <- reactive({
    req(ipa_funct(), input$col_category)
    df <- ipa_funct()

# Filter 1: pathway name — literal substring match (no regex), case-insensitive
    kw_path <- split_keywords(input$pathway_keyword_filter)
    if (length(kw_path) > 0) {
      names_lc <- tolower(as.character(df[[input$col_category]]))
      keep <- Reduce(`|`, lapply(kw_path, function(k) grepl(k, names_lc, fixed = TRUE)))
      df <- df[keep, , drop = FALSE]
    }

    # Filter 2: molecule / gene symbol — exact whole-symbol match, case-insensitive
    # (so MAPK1 does not also match MAPK10 or MAPK11); "*" is a wildcard (SYN* -> SYN1, SYN3, ...)
    kw_mol <- split_keywords(input$molecule_keyword_filter)
    if (length(kw_mol) > 0) {
      validate(need(has_members(), members_msg))
      req(input$mol_separator)
      exact    <- kw_mol[!grepl("*", kw_mol, fixed = TRUE)]
      wild     <- kw_mol[grepl("*", kw_mol, fixed = TRUE)]
      # Escape regex characters, then turn * into .* and anchor to the whole symbol
      wild_pat <- if (length(wild) > 0)
        paste0("^(", paste(gsub("\\*", ".*", gsub("([.^$|()\\[\\]{}+?\\\\])", "\\\\\\1", wild, perl = TRUE)),
                           collapse = "|"), ")$")
      members <- strsplit(as.character(df[[input$col_molecules]]),
                          input$mol_separator, fixed = TRUE)
      keep <- vapply(members, function(m) {
        m <- tolower(trimws(m))
        any(m %in% exact) || (!is.null(wild_pat) && any(grepl(wild_pat, m)))
      }, logical(1))
      df <- df[keep, , drop = FALSE]
    }
    df
  })
  
  # Summary badge for keyword filter panel
  output$ui_keyword_filter_summary <- renderUI({
    req(ipa_funct())
    n_total    <- nrow(ipa_funct())
    kw_path    <- input$pathway_keyword_filter
    kw_mol     <- input$molecule_keyword_filter
    path_active <- !is.null(kw_path) && nchar(trimws(kw_path)) > 0
    mol_active  <- !is.null(kw_mol)  && nchar(trimws(kw_mol))  > 0
    if (!path_active && !mol_active) {
      helpText(paste0("\u2139\ufe0f  Using all ", n_total,
                      " pathways (no keyword filter active)."))
    } else {
      n_kept <- nrow(ipa_funct_filtered())
      lines  <- c()
      if (path_active) lines <- c(lines, paste0("pathway name: ", trimws(kw_path)))
      if (mol_active)  lines <- c(lines, paste0("molecule: ",     trimws(kw_mol)))
      helpText(paste0(
        "\u2705  ", n_kept, " / ", n_total,
        " pathways retained  \u2014  ",
        paste(lines, collapse = "  AND  ")
      ))
    }
  })
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
    default <- detect_category_col(ipa_funct())
    selectInput("col_category", label = paste0("Category column (detected: '", default, "')"),
                choices = colnames(ipa_funct()), selected = default)
  })
  
  output$ui_col_molecules <- renderUI({
    req(ipa_funct())
    default <- detect_molecules_col(ipa_funct())
    tagList(
      selectInput("col_molecules",
                  label    = paste0("Molecules column (detected: '",
                                    if (is.null(default)) "none" else default, "')"),
                  choices  = c("(none \u2014 no member list)" = "", colnames(ipa_funct())),
                  selected = if (is.null(default)) "" else default),
      if (is.null(default))
        helpText("\u26a0 No column with member genes / proteins / metabolites was found.",
                 "The Bar Plot works; Chord Diagram, Molecule Frequency, Heatmap and Boxplot",
                 "need a column listing each pathway's members.")
    )
  })

  # TRUE when a molecules (member list) column is selected
  has_members <- reactive(!is.null(input$col_molecules) && nzchar(input$col_molecules))
  members_msg <- paste("This view needs a column listing each pathway's member genes / proteins /",
                       "metabolites. Select it under 'Molecules column' in the sidebar.")

  # When the molecules column changes, pre-select the separator it actually uses
  # (the user can still change it afterwards)
  observeEvent(input$col_molecules, {
    req(ipa_funct(), has_members())
    dens <- separator_density(ipa_funct()[[input$col_molecules]])
    if (max(dens) > 0)
      updateRadioButtons(session, "mol_separator", selected = names(dens)[which.max(dens)])
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
    default <- detect_pvalue_col(ipa_funct())
    selectInput("col_pvalue", label = paste0("P-value or -log10(p-value) column (X-axis) (detected: '", default, "')"),
                choices = colnames(ipa_funct()), selected = default)
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
  
  # ── 6. Dynamic color inputs: one per group (uses filtered annotation) ────
  output$ui_group_colors <- renderUI({
    req(annotation(), input$col_group)
    groups <- unique(as.character(annotation()[[input$col_group]]))
    default_colors <- c("purple","darkorange","steelblue","forestgreen",
                        "firebrick","gold","hotpink","cyan4","sienna","gray40")
    tagList(lapply(seq_along(groups), function(i) {
      textInput(inputId = paste0("color_group_", i),
                label   = paste0("Color \u2014 ", groups[i], " (hex or name)"),
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
  # Server-side selectize: only matches for the typed text are sent to the browser
  observe({
    choices <- tryCatch({
      req(ipa_funct_filtered(), input$col_category)
      as.character(unique(ipa_funct_filtered()[[input$col_category]]))
    }, error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      character(0)
    })
    current <- isolate(input$selected_category)
    selected <- if (!is.null(current) && current %in% choices) current else head(choices, 1)
    updateSelectizeInput(session, "selected_category",
                         label    = paste0("Select: ", input$col_category, " (", length(choices), ")"),
                         choices  = choices, selected = selected, server = TRUE)
  })
  
  # ── 9. Gene list from selected category ──────────────────────────────────
  genes_interest <- reactive({
    validate(need(has_members(), members_msg))
    req(ipa_funct(), input$selected_category, input$col_category,
        input$col_molecules, input$mol_separator)
    ipa_funct() %>%
      filter(.data[[input$col_category]] == input$selected_category) %>%
      pull(.data[[input$col_molecules]]) %>%
      as.character() %>%
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
      "Molecules in list: ", length(genes_interest()),                      "\n",
      "Molecules matched: ", matched,                                       "\n",
      "Sample columns   : ", paste(head(sample_cols_resolved(), 5), collapse = ", "),
      if (length(sample_cols_resolved()) > 5) paste0(", ... (", length(sample_cols_resolved()), " in total)") else ""
    )
  })
  
  # ══════════════════════════════════════════════════════════════════════════
  # BAR / BUBBLE PLOT
  # ══════════════════════════════════════════════════════════════════════════
  output$ui_bubble_size_col <- renderUI({
    req(ipa_funct())
    cols     <- colnames(ipa_funct())
    # Bubble size must be numeric (skips text like Enrichr "Overlap" = "4/700")
    num_cols <- cols[sapply(ipa_funct(), is_numeric_col)]
    keywords <- c("count","size","gene","overlap","number","n_gene","ngene","ratio")
    detected <- num_cols[sapply(tolower(num_cols), function(x)
      any(sapply(keywords, function(k) grepl(k, x, fixed = TRUE))))]
    default  <- if (length(detected) > 0) detected[1] else if (length(num_cols) > 0) num_cols[1] else cols[1]
    selectInput("bubble_size_col",
                label   = paste0("Bubble SIZE column (detected: '", default, "')"),
                choices = cols, selected = default)
  })
  
  output$ui_bubble_color_col <- renderUI({
    req(ipa_funct())
    cols    <- colnames(ipa_funct())
    # Only numeric columns can color bubbles (text columns would give all-gray bubbles)
    num_cols <- cols[sapply(ipa_funct(), is_numeric_col)]
    default <- detect_pvalue_col(ipa_funct())
    if (!default %in% num_cols && length(num_cols) > 0) default <- num_cols[1]
    selectInput("bubble_color_col",
                label   = paste0("Bubble COLOR column (detected: '", default, "')"),
                choices = num_cols, selected = default)
  })
  
  barplot_data <- reactive({
    req(ipa_funct(), input$col_category, input$col_pvalue, pvalue_transform())
    df       <- ipa_funct_filtered()
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
      "Total pathways   : ", nrow(ipa_funct_filtered()),
      if (nrow(ipa_funct_filtered()) < nrow(ipa_funct()))
        paste0(" (keyword-filtered from ", nrow(ipa_funct()), ")") else "",                  "\n",
      "Showing top      : ", nrow(barplot_data())
    )
    if (input$barplot_type == "bar") {
      fill <- bar_fill_info(barplot_data())
      if (!is.null(fill))
        base <- paste0(base, "\nBar color col   : ", fill$label, " (", fill$scale, ")",
                       if (fill$col == input$col_pvalue) " -- same as the x-axis" else "",
                       if (any(is.na(fill$vals))) paste0("; ", sum(is.na(fill$vals)), " pathway(s) without a value shown in gray") else "")
    }
    if (input$barplot_type == "bubble" &&
        !is.null(input$bubble_size_col) && !is.null(input$bubble_color_col)) {
      base <- paste0(base,
                     "\nBubble size col  : ", input$bubble_size_col,
                     "\nBubble color col: ", input$bubble_color_col)
    }
    base
  })
  
  # Optional numeric column used to color the bars ("thermal" bar plot)
  output$ui_bar_fill_col <- renderUI({
    req(ipa_funct())
    num_cols <- colnames(ipa_funct())[sapply(ipa_funct(), is_numeric_col)]
    selectInput("bar_fill_col", "Color bars by",
                choices  = c("None (single color)" = "", num_cols),
                selected = "")
  })

  # TRUE for columns whose name says they hold p-values / q-values (not ratios etc.)
  is_pvalue_name <- function(nm) {
    keys <- c("padj", "p.adjust", "p_adjust", "adjusted p", "adjusted_p", "adj.p", "adj_p",
              "p_adj", "fdr", "benjamini", "qvalue", "q.value", "q-value", "q_value", "bonferroni",
              "pvalue", "p.value", "p-value", "p_value", "pval", "p value")
    nm <- tolower(nm)
    any(sapply(keys, function(k) grepl(k, nm, fixed = TRUE))) && !grepl("log", nm, fixed = TRUE)
  }

  # Values, legend label and scale type for the bar color column
  bar_fill_info <- function(df) {
    col <- input$bar_fill_col
    if (is.null(col) || col == "" || !col %in% colnames(df)) return(NULL)
    vals  <- suppressWarnings(as.numeric(df[[col]]))
    label <- col
    if (is_pvalue_name(col) && all(vals >= 0 & vals <= 1, na.rm = TRUE)) {
      vals  <- -log10(pmax(vals, .Machine$double.xmin))
      label <- paste0("-log10(", col, ")")
    }
    # Signed metrics (z-score, NES, ...) are diverging even if the bars shown share one sign,
    # so white always means 0 and the strongest values are always the darkest
    signed <- any(suppressWarnings(as.numeric(ipa_funct()[[col]])) < 0, na.rm = TRUE)
    scale  <- switch(if (is.null(input$bar_fill_scale)) "auto" else input$bar_fill_scale,
                     gradient = "gradient", diverging = "diverging",
                     if (signed) "diverging" else "gradient")
    list(col = col, vals = vals, label = label, scale = scale)
  }

  # Shortens long pathway names for axis labels. Wrap breaks at spaces, "_" or "-"
  # (GSEA names such as GOBP_..._..._ have no spaces); at most `max_lines` lines.
  # Save a plot as a 300 dpi PNG of the user-set size; `scaling` enlarges/shrinks all
  # text and lines the same way for ggplot2, pheatmap (grid) and circlize (base) plots
  dl_png <- function(file, key, draw) {
    d <- dl_defaults[[key]]
    num <- function(id, def) { v <- input[[id]]; if (is.null(v) || is.na(v) || v <= 0) def else v }
    sc  <- num("dl_text_scale", 1)
    ragg::agg_png(file, width = num(paste0("dlw_", key), d[1]), height = num(paste0("dlh_", key), d[2]),
                  units = "in", res = 300, scaling = sc)
    on.exit(dev.off(), add = TRUE)
    draw()
    invisible(TRUE)
  }

  # Pathway names as displayed: optional user prefixes (e.g. "GOBP_") removed from the start
  display_name <- function(x) {
    x <- as.character(x)
    pre <- trimws(strsplit(if (is.null(input$strip_prefix)) "" else input$strip_prefix, ",")[[1]])
    pre <- pre[nzchar(pre)]
    for (p in pre) {
      hit <- startsWith(x, p) & nchar(x) > nchar(p)
      x[hit] <- substring(x[hit], nchar(p) + 1)
    }
    x
  }

  shorten_labels <- function(x, width = 50, mode = "wrap", max_lines = 3) {
    x <- as.character(x)
    if (is.null(mode) || mode == "full" || is.null(width) || is.na(width)) return(x)
    width <- max(10, width)
    ell   <- "\u2026"
    trunc1 <- function(s) if (nchar(s) > width) paste0(substr(s, 1, width - 1), ell) else s
    if (mode == "truncate") return(vapply(x, trunc1, character(1), USE.NAMES = FALSE))
    vapply(x, function(s) {
      if (nchar(s) <= width) return(s)
      tokens <- regmatches(s, gregexpr("[^ _-]+[ _-]*", s))[[1]]
      lines <- character(0); cur <- ""
      for (tk in tokens) {
        if (nchar(cur) > 0 && nchar(cur) + nchar(trimws(tk, "right")) > width) {
          lines <- c(lines, cur); cur <- tk
        } else cur <- paste0(cur, tk)
      }
      lines <- trimws(c(lines, cur), "right")
      lines <- vapply(lines, trunc1, character(1), USE.NAMES = FALSE)   # a single very long word
      if (length(lines) > max_lines)
        lines <- c(lines[seq_len(max_lines - 1)], trunc1(paste0(lines[max_lines], ell)))
      paste(lines, collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
  }
  pathway_label_scale <- function()
    scale_y_discrete(labels = function(x) shorten_labels(display_name(x), input$label_width, input$label_mode))

  make_barplot <- function() {
    req(barplot_data(), input$barplot_type)
    df      <- barplot_data()
    x_label <- if (pvalue_transform() == "raw")
      paste0("-log10(", input$col_pvalue, ")") else input$col_pvalue
    
    if (input$barplot_type == "bar") {
      bar_theme <- list(
        labs(title = paste0("Top ", nrow(df), " Enriched Pathways / Functions"),
             x = x_label, y = NULL),
        pathway_label_scale(),
        theme_bw(base_size = 12),
        theme(axis.text.y         = element_text(size = 10, lineheight = 0.8),
              axis.text.x         = element_text(size = 10),
              plot.title          = element_text(face = "bold", size = 13),
              panel.grid.major.y  = element_blank())
      )
      fill <- bar_fill_info(df)
      if (is.null(fill)) {
        bar_color <- if (is.null(input$bar_color) || input$bar_color == "") "steelblue" else input$bar_color
        return(ggplot(df, aes(x = .neg_log10, y = .category)) +
                 geom_col(fill = bar_color, width = 0.7) + bar_theme)
      }
      
      # Color bars by a numeric column
      df$.bar_fill <- fill$vals
      p <- ggplot(df, aes(x = .neg_log10, y = .category, fill = .bar_fill)) +
        geom_col(width = 0.7, colour = "grey40", linewidth = 0.2) + bar_theme
      if (fill$scale == "diverging") {
        neg <- if (is.null(input$bar_neg_color) || input$bar_neg_color == "") "steelblue"  else input$bar_neg_color
        pos <- if (is.null(input$bar_pos_color) || input$bar_pos_color == "") "darkorange" else input$bar_pos_color
        lim <- max(abs(fill$vals), na.rm = TRUE)
        if (!is.finite(lim) || lim == 0) lim <- 1
        p + scale_fill_gradient2(low = neg, mid = "white", high = pos, midpoint = 0,
                                 limits = c(-lim, lim), na.value = "grey80", name = fill$label)
      } else {
        palette <- if (is.null(input$bar_palette) || input$bar_palette == "") "Blues" else input$bar_palette
        if (palette %in% c("viridis", "magma", "plasma"))
          p + scale_fill_viridis_c(option = palette, na.value = "grey80", name = fill$label)
        else
          p + scale_fill_distiller(palette = palette, direction = 1, na.value = "grey80", name = fill$label)
      }
      
    } else {
      req(input$bubble_size_col, input$bubble_color_col)
      size_raw  <- suppressWarnings(as.numeric(ipa_funct()[[input$bubble_size_col]]))
      color_raw <- suppressWarnings(as.numeric(ipa_funct()[[input$bubble_color_col]]))
      df$.bub_size  <- size_raw [match(as.character(df[[input$col_category]]),
                                       as.character(ipa_funct()[[input$col_category]]))]
      df$.bub_color <- color_raw[match(as.character(df[[input$col_category]]),
                                       as.character(ipa_funct()[[input$col_category]]))]
      # -log10 only for p-value-like columns (a 0-1 column such as IPA "Ratio" stays as is)
      color_needs_transform <- is_pvalue_name(input$bubble_color_col) &&
        all(df$.bub_color >= 0 & df$.bub_color <= 1, na.rm = TRUE)
      if (color_needs_transform) df$.bub_color <- -log10(df$.bub_color)
      color_label <- if (color_needs_transform)
        paste0("-log10(", input$bubble_color_col, ")") else input$bubble_color_col
      palette     <- if (is.null(input$bubble_palette) || input$bubble_palette == "") "Blues" else input$bubble_palette
      use_viridis <- palette %in% c("viridis", "magma", "plasma")
      p <- ggplot(df, aes(x = .neg_log10, y = .category,
                          size = .bub_size, fill = .bub_color)) +
        geom_point(shape = 21, alpha = 0.85, colour = "grey30") +
        scale_size(name  = input$bubble_size_col, range = input$bubble_size_range) +
        pathway_label_scale() +
        labs(title = paste0("Top ", nrow(df), " Enriched Pathways / Functions"),
             x = x_label, y = NULL, fill = color_label) +
        theme_bw(base_size = 12) +
        theme(axis.text.y        = element_text(size = 10, lineheight = 0.8),
              axis.text.x        = element_text(size = 10),
              plot.title         = element_text(face = "bold", size = 13),
              panel.grid.major.y = element_blank(),
              legend.position    = "right")
      signed <- any(df$.bub_color < 0, na.rm = TRUE)
      scale  <- switch(if (is.null(input$bubble_fill_scale)) "auto" else input$bubble_fill_scale,
                       gradient = "gradient", diverging = "diverging",
                       if (signed) "diverging" else "gradient")
      if (scale == "diverging") {
        neg <- if (is.null(input$bubble_neg_color) || input$bubble_neg_color == "") "steelblue"  else input$bubble_neg_color
        pos <- if (is.null(input$bubble_pos_color) || input$bubble_pos_color == "") "darkorange" else input$bubble_pos_color
        lim <- max(abs(df$.bub_color), na.rm = TRUE)
        if (!is.finite(lim) || lim == 0) lim <- 1
        p <- p + scale_fill_gradient2(low = neg, mid = "white", high = pos, midpoint = 0,
                                      limits = c(-lim, lim), na.value = "grey80", name = color_label)
      } else if (use_viridis) {
        p <- p + scale_fill_viridis_c(option = palette, na.value = "grey80", name = color_label)
      } else {
        p <- p + scale_fill_distiller(palette = palette, direction = 1, na.value = "grey80", name = color_label)
      }
      p
    }
  }
  
  output$barplot <- renderPlot({ req(barplot_data()); make_barplot() },
                               height = function() input$bar_plot_height)
  
  dl_barplot_handler <- function(file) {
    dl_png(file, "bar", function() print(make_barplot()))
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
    validate(need(has_members(), members_msg))
    req(ipa_funct(), input$col_category, input$col_molecules,
        input$col_pvalue, pvalue_transform(), input$mol_separator)
    df <- ipa_funct_filtered()
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
      "Link curve depth : ", input$chord_inner_radius,
      if (n_distinct(cd$molecules) > 100)
        "\nTip              : many molecules — raise 'Minimum appearances' for readable labels"
      else ""
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
    # Shades from a darker to a lighter tone of each base color (the old ramp to "lightyellow"
    # made the last sectors and their links nearly white)
    shade_ramp  <- function(col, n) {
      dark  <- colorRampPalette(c(col, "black"))(5)[2]
      light <- colorRampPalette(c(col, "white"))(5)[3]
      colorRampPalette(c(dark, col, light))(n)
    }
    path_cols   <- shade_ramp(pathway_color,  n_path)
    mol_cols    <- shade_ramp(molecule_color, n_mol)
    grid_colors <- c(setNames(path_cols, rownames(mat)), setNames(mol_cols, colnames(mat)))
    lab_width   <- if (is.null(input$chord_label_width) || is.na(input$chord_label_width)) 35 else input$chord_label_width
    circos.clear()
    par(mar = c(2, 2, 3, 2), xpd = NA)   # xpd = NA: labels may extend into the margins instead of being clipped
    # A fixed 2-degree gap per sector fails when there are more than ~180 sectors
    # (circlize: "gap.degree is too large"), so the gaps take at most 1/4 of the circle.
    gap <- min(2, 90 / (n_path + n_mol))
    # Extra canvas around the circle so outward-facing labels are not clipped (user-adjustable:
    # larger values leave more room for long or large labels and draw a smaller circle)
    canvas <- if (is.null(input$chord_canvas) || is.na(input$chord_canvas)) 1.5 else input$chord_canvas
    circos.par("gap.after" = gap, "start.degree" = 90,
               "canvas.xlim" = c(-canvas, canvas), "canvas.ylim" = c(-canvas, canvas))
    chordDiagram(mat, h = input$chord_inner_radius, grid.col = grid_colors,
                 transparency = 0.4, annotationTrack = "grid",
                 preAllocateTracks = list(track.height = 0.08), directional = FALSE)
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      sector_name <- get.cell.meta.data("sector.index")
      xlim        <- get.cell.meta.data("xlim")
      ylim        <- get.cell.meta.data("ylim")
      is_pathway  <- sector_name %in% rownames(mat)
      lab_mode    <- if (identical(input$chord_label_mode, "wrap")) "wrap" else "truncate"
      label       <- if (is_pathway) shorten_labels(display_name(sector_name), lab_width, lab_mode) else sector_name
      circos.text(mean(xlim), ylim[1] + 0.1, label,
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                  cex = if (is_pathway) input$chord_label_size_path else input$chord_label_size_mol)
    }, bg.border = NA)
    if (!isFALSE(input$chord_show_title))
      title(paste0("Chord Diagram \u2014 Top ", n_distinct(cd$pathway), " Pathways & Molecules"), cex.main = 1.1)
    legend("bottomleft", legend = c("Pathways / Functions", "Molecules"),
           fill = c(pathway_color, molecule_color), border = NA, bty = "n", cex = 0.85)
    circos.clear()
  }
  
  output$chord_plot <- renderPlot({ req(chord_data()); make_chord() },
                                  height = function() input$chord_plot_height)
  dl_chord_handler  <- function(file) dl_png(file, "chord", make_chord)
  output$download_chord      <- downloadHandler(filename = function() "Chord_Diagram_Pathways_Molecules.png", content = dl_chord_handler)
  output$download_chord_main <- downloadHandler(filename = function() "Chord_Diagram_Pathways_Molecules.png", content = dl_chord_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # MOLECULE FREQUENCY BAR PLOT (All pathways OR Top N pathways)
  # ══════════════════════════════════════════════════════════════════════════
  gene_freq_data <- reactive({
    validate(need(has_members(), members_msg))
    req(ipa_funct(), input$col_category, input$col_molecules,
        input$mol_separator, input$freq_pathway_scope)
    
    # Keyword-filtered pathways, or the full table if the user unticks the option
    df <- if (isTRUE(input$freq_use_keyword_filter)) ipa_funct_filtered() else ipa_funct()

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
    
    use_kw   <- isTRUE(input$freq_use_keyword_filter)
    base_df  <- if (use_kw) ipa_funct_filtered() else ipa_funct()
    kw_line  <- if (!keyword_filter_active()) "Keyword filter  : none active\n"
                else if (use_kw) paste0("Keyword filter  : applied (", nrow(ipa_funct_filtered()),
                                        " of ", nrow(ipa_funct()), " pathways)\n")
                else paste0("Keyword filter  : ignored (all ", nrow(ipa_funct()), " pathways used)\n")
    scope_line <- if (input$freq_pathway_scope == "all") {
      paste0("Pathways used   : ALL (", n_distinct(base_df[[input$col_category]]), ")")
    } else {
      paste0("Pathways used   : TOP ", input$freq_top_n_pathways,
             " (ranked by ", input$col_pvalue, ")")
    }
    
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter : ", input$filter_col, " = ",
             paste(input$filter_vals, collapse = ", "),
             " (not applied: frequency uses the pathway table only)\n")
    else ""

    paste0(
      filter_line,
      kw_line,
      scope_line, "\n",
      "Molecules total : ", nrow(gf), "\n",
      "Molecules with >= ", input$freq_min_pathways, " pathways: ", nrow(gf_filt), "\n",
      "Showing top     : ", min(input$freq_top_n_genes, nrow(gf_filt)),
      freq_lfc_line()
    )
  })
  
  # log2 fold changes from the DE file, used to color the frequency bars (NULL if not available)
  freq_lfc <- reactive({
    if (!isTRUE(input$freq_color_lfc) || is.null(input$file_de)) return(NULL)
    lfc <- tryCatch(de_data(), error = function(e) NULL)
    if (is.null(lfc) || nrow(lfc) == 0) return(NULL)
    lfc <- lfc[!duplicated(lfc$gene), c("gene", "lfc")]
    setNames(lfc$lfc, lfc$gene)
  })

  freq_lfc_line <- function() {
    if (!isTRUE(input$freq_color_lfc)) return("")
    lfc <- freq_lfc()
    if (is.null(lfc)) return("\nBar color      : single color (no differential analysis file loaded)")
    gf <- gene_freq_data() %>% filter(pathway_count >= input$freq_min_pathways) %>%
      slice_head(n = input$freq_top_n_genes)
    n_match <- sum(gf$molecules %in% names(lfc))
    paste0("\nBar color      : log2 fold change (", input$de_col_lfc, "); ",
           n_match, " of ", nrow(gf), " shown molecules found in the DE file",
           if (n_match < nrow(gf)) " (others in gray)" else "")
  }

  make_gene_freq_plot <- function() {
    req(gene_freq_data(), input$freq_min_pathways, input$freq_top_n_genes)
    
    gf <- gene_freq_data() %>%
      filter(pathway_count >= input$freq_min_pathways) %>%
      slice_head(n = input$freq_top_n_genes)
    
    validate(need(nrow(gf) > 0,
                  "No molecules meet the frequency cutoff. Lower the cutoff or expand pathway scope."))
    
    lfc <- freq_lfc()
    p <- ggplot(gf, aes(x = reorder(molecules, -pathway_count), y = pathway_count))
    if (is.null(lfc)) {
      p <- p + geom_col(fill = "steelblue", color = "white")
    } else {
      # Diverging scale centered at 0 with the volcano plot's Down / Up colors
      up_col   <- if (is.null(input$volcano_up_color)   || input$volcano_up_color   == "") "firebrick" else input$volcano_up_color
      down_col <- if (is.null(input$volcano_down_color) || input$volcano_down_color == "") "steelblue" else input$volcano_down_color
      gf$lfc <- unname(lfc[gf$molecules])
      lim <- max(abs(gf$lfc), na.rm = TRUE)
      if (!is.finite(lim) || lim == 0) lim <- 1
      p <- p + geom_col(data = gf, aes(fill = lfc), color = "white") +
        scale_fill_gradient2(low = down_col, mid = "grey95", high = up_col, midpoint = 0,
                             limits = c(-lim, lim), na.value = "grey60",
                             name = "log2 fold\nchange")
    }
    p +
      labs(
        title = paste0("Molecules appearing in \u2265 ", input$freq_min_pathways, " pathways"),
        x = "Molecule",
        y = "Number of pathways"
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  output$gene_freq_plot <- renderPlot({
    make_gene_freq_plot()
  }, height = function() input$freq_plot_height)
  
  dl_gene_freq_filename <- function() {
    paste0("Molecule_Frequency_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  }
  dl_gene_freq_handler <- function(file) {
    dl_png(file, "freq", function() print(make_gene_freq_plot()))
  }
  
  output$download_gene_freq      <- downloadHandler(filename = dl_gene_freq_filename, content = dl_gene_freq_handler)
  output$download_gene_freq_main <- downloadHandler(filename = dl_gene_freq_filename, content = dl_gene_freq_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # HEATMAP  — now uses filtered annotation + filtered sample columns
  # ══════════════════════════════════════════════════════════════════════════
  make_heatmap <- function(gene_filter, title_text, log2_transform = FALSE) {
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
    if (log2_transform) {
      validate(need(all(mat >= 0, na.rm = TRUE),
                    "Data contain negative values (already log-transformed?). Untick 'Apply log2(x + 1)'."))
      mat <- log2(mat + 1)
    }
    
    row_vars <- apply(mat, 1, var, na.rm = TRUE)
    mat      <- mat[!is.na(row_vars) & row_vars > 0, , drop = FALSE]
    
    validate(need(nrow(mat) > 1,
                  "All matched genes have zero variance across samples \u2014 cannot draw heatmap."))
    
    mat_scaled <- t(scale(t(mat)))
     
    breaks_sym <- seq(-input$heatmap_scale_limit, input$heatmap_scale_limit, length.out = 101)
    
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
      cluster_rows      = isTRUE(input$heatmap_cluster_rows),   # ← NEW
      cluster_cols      = isTRUE(input$heatmap_cluster_cols),   # ← NEW
      show_rownames     = nrow(mat) <= 100, show_colnames = ncol(mat) <= 60,
      scale             = "none",
      color             = colorRampPalette(c("darkblue", "white", "darkred"))(100),
      breaks            = breaks_sym,
      border_color      = NA,
      fontsize_row      = if (is.null(input$heatmap_fontsize_row) || is.na(input$heatmap_fontsize_row)) 8 else input$heatmap_fontsize_row,
      fontsize_col      = 10,
      main              = paste0(title_text, filter_subtitle)
    )
  }
  
  output$heatmap_plot <- renderPlot({
    req(norm_prot(), annotation(), genes_interest(), input$selected_category)
    make_heatmap(gene_filter = genes_interest(),
                 title_text  = paste0(shorten_labels(display_name(input$selected_category), 30, "wrap", max_lines = 4), "\n(Z-score, by Group)"))
  }, height = function() input$plot_height)
  
  dl_single_filename <- function()
    paste0("Heatmap_", gsub("[^A-Za-z0-9]", "_", input$selected_category), ".png")
  dl_single_handler  <- function(file) {
    dl_png(file, "heat", function()
      make_heatmap(gene_filter = genes_interest(),
                   title_text  = paste0(shorten_labels(display_name(input$selected_category), 30, "wrap", max_lines = 4), "\n(Z-score, by Group)")))
  }
  output$download_single      <- downloadHandler(filename = dl_single_filename, content = dl_single_handler)
  output$download_single_main <- downloadHandler(filename = dl_single_filename, content = dl_single_handler)
  
  # Batch export is limited to keyword-filtered pathways, capped to keep downloads fast
  max_batch_heatmaps <- 100
  batch_categories <- reactive({
    req(ipa_funct_filtered(), input$col_category)
    unique(ipa_funct_filtered()[[input$col_category]])
  })

  batch_download_ui <- function(id, prefix = "") {
    if (is.null(input$file_ipa_funct)) return(NULL)
    if (!keyword_filter_active())
      return(helpText("\u2139\ufe0f  Batch download of heatmaps is available after applying a",
                      "keyword filter (Keyword Filter \u2014 Pathways panel)."))
    n <- length(batch_categories())
    if (n == 0)
      return(helpText("\u26a0 No pathways match the keyword filter."))
    if (n > max_batch_heatmaps)
      return(helpText(paste0("\u26a0 ", n, " pathways match the keyword filter. Narrow it to ",
                             max_batch_heatmaps, " or fewer to enable batch download.")))
    downloadButton(id, paste0(prefix, "Download Filtered Heatmaps (", n, ") (.zip)"))
  }
  output$ui_download_all      <- renderUI(batch_download_ui("download_all"))
  output$ui_download_all_main <- renderUI(batch_download_ui("download_all_main", "\u2b07 "))

  download_all_handler <- function(file) {
    validate(need(has_members(), members_msg))
    req(input$col_category, input$col_molecules, input$col_gene, input$mol_separator)
    cats <- batch_categories()
    validate(need(keyword_filter_active(), "Apply a keyword filter before batch-downloading heatmaps."),
             need(length(cats) <= max_batch_heatmaps,
                  paste0("Too many pathways (", length(cats), "); narrow the keyword filter.")))
# Fresh sub-directory so files from earlier downloads are never mixed in
    tmp_dir <- tempfile("heatmaps_")
    dir.create(tmp_dir)
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    saved   <- c()
    for (cat_val in cats) {
genes <- ipa_funct() %>%
        filter(.data[[input$col_category]] == cat_val) %>%
        pull(.data[[input$col_molecules]]) %>%
        as.character() %>%
        strsplit(input$mol_separator, fixed = TRUE) %>%
        unlist() %>% trimws() %>% unique()
      if (length(genes) == 0) next
      data <- norm_prot()[as.character(norm_prot()[[input$col_gene]]) %in% genes, ]
      if (nrow(data) <= 1) next
      clean    <- gsub("[^A-Za-z0-9]", "_", cat_val)
      png_path <- file.path(tmp_dir, paste0("Heatmap_", clean, ".png"))
      # A category that cannot be drawn (e.g. all genes zero-variance) is skipped
      # instead of aborting the whole archive; the device is always closed.
      ok <- tryCatch(dl_png(png_path, "heat", function()
        make_heatmap(gene_filter = genes, title_text = paste0(shorten_labels(display_name(cat_val), 30, "wrap", max_lines = 4), "\n(Z-score, by Group)"))),
        error = function(e) FALSE)
      if (ok) saved <- c(saved, png_path) else unlink(png_path)
    }
    validate(need(length(saved) > 0,
                  "No heatmaps could be generated for any category."))
    zip(zipfile = file, files = saved, flags = "-j")
  }
  dl_all_filename <- function() paste0("Filtered_Heatmaps_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
  output$download_all      <- downloadHandler(filename = dl_all_filename, content = download_all_handler)
  output$download_all_main <- downloadHandler(filename = dl_all_filename, content = download_all_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # BOXPLOT / VIOLIN PLOT  — uses filtered annotation + filtered samples
  # ══════════════════════════════════════════════════════════════════════════
  all_molecules <- reactive({
    validate(need(has_members(), members_msg))
    req(ipa_funct(), input$col_molecules, input$mol_separator)
    ipa_funct()[[input$col_molecules]] %>%
      as.character() %>%
      strsplit(input$mol_separator, fixed = TRUE) %>%
      unlist() %>% trimws() %>% unique() %>% sort()
  })
  
  boxplot_molecules <- reactive({
    req(all_molecules(), norm_prot(), input$col_gene)
    gene_col <- as.character(norm_prot()[[input$col_gene]])
    all_molecules()[all_molecules() %in% gene_col]
  })

  # Server-side selectize: only matches for the typed text are sent to the browser
  observe({
    mols <- tryCatch(boxplot_molecules(), error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)   # keep normal req() behavior
      character(0)                                     # bad column choice -> empty list, no crash
    })
    current <- isolate(input$boxplot_protein)
    selected <- if (!is.null(current) && current %in% mols) current else head(mols, 1)
    updateSelectizeInput(session, "boxplot_protein",
                         label    = paste0("Search protein / gene / metabolite (", length(mols), " matched)"),
                         choices  = mols, selected = selected, server = TRUE)
  })

  output$ui_boxplot_protein_msg <- renderUI({
    mols <- tryCatch(boxplot_molecules(), error = function(e) {
      if (inherits(e, "shiny.silent.error")) stop(e)
      character(0)
    })
    if (length(mols) == 0)
      helpText("\u26a0 No molecules from the pathway file found in the normalized data.",
               "Check the Molecules column, the separator, and the Gene/Protein ID column.")
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
      "Selected molecule: ", input$boxplot_protein,                    "\n",
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
        # trim = TRUE: the violin ends at the observed min / max (no tails below 0 for counts)
        geom_violin(alpha = 0.5, trim = TRUE, width = 0.8) +
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
           x = input$col_group, y = "Normalized abundance / counts") +
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
    dl_png(file, "box", function() print(make_boxplot()))
  
  output$download_boxplot_sidebar <- downloadHandler(filename = dl_boxplot_filename, content = dl_boxplot_handler)
  output$download_boxplot_main    <- downloadHandler(filename = dl_boxplot_filename, content = dl_boxplot_handler)
  
  # ══════════════════════════════════════════════════════════════════════════
  # DIFFERENTIAL ANALYSIS — volcano, PCA, DE heatmap
  # ══════════════════════════════════════════════════════════════════════════
  de_raw <- reactive({ req(input$file_de); read.csv(input$file_de$datapath, check.names = FALSE) })
  
  # First column whose lower-cased name contains a keyword; keywords are tried in priority order
  detect_col <- function(cols, keywords, fallback) {
    for (k in keywords) {
      hit <- cols[grepl(k, tolower(cols), fixed = TRUE)]
      if (length(hit) > 0) return(hit[1])
    }
    fallback
  }
  
  output$ui_de_cols <- renderUI({
    req(de_raw())
    cols <- colnames(de_raw())
    # Prefer the same ID column name as the normalized data, so identifiers match
    gene_default <- if (!is.null(input$col_gene) && input$col_gene %in% cols) input$col_gene else
      detect_col(cols, c("symbol", "gene_name", "genename", "gene", "protein", "uniprot", "accession",
                         "metabolite", "compound", "hmdb", "kegg", "feature", "name", "id"), cols[1])
    lfc_default  <- detect_col(cols, c("log2foldchange", "log2fc", "logfc", "log2_fc", "lfc", "fold"), cols[2])
    padj_default <- detect_col(cols, c("padj", "p.adjust", "adj.p", "adj_p", "fdr", "qvalue", "q.value", "pvalue", "p.value"), cols[3])
    tagList(
      selectInput("de_col_gene", paste0("Feature ID column \u2014 gene / protein / metabolite (detected: '", gene_default, "')"),
                  choices = cols, selected = gene_default),
      selectInput("de_col_lfc",  paste0("log2 fold-change column (detected: '", lfc_default, "')"),
                  choices = cols, selected = lfc_default),
      selectInput("de_col_padj", paste0("Adjusted p-value column (detected: '", padj_default, "')"),
                  choices = cols, selected = padj_default)
    )
  })
  
  # DE table with standardised columns and Up / Down / NS status
  de_data <- reactive({
    req(de_raw(), input$de_col_gene, input$de_col_lfc, input$de_col_padj,
        input$de_lfc_cut, input$de_padj_cut)
    de <- de_raw()
    df <- data.frame(
      gene = trimws(as.character(de[[input$de_col_gene]])),
      lfc  = suppressWarnings(as.numeric(de[[input$de_col_lfc]])),
      padj = suppressWarnings(as.numeric(de[[input$de_col_padj]])),
      stringsAsFactors = FALSE
    )
    # NA padj (e.g. DESeq2 independent filtering) cannot be tested -> dropped
    df <- df[!is.na(df$lfc) & !is.na(df$padj), , drop = FALSE]
    validate(need(nrow(df) > 0, "No rows with numeric fold change and adjusted p-value. Check the DE column selection."))
    # padj of exactly 0 would give Inf on the -log10 scale
    df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))
    df$status <- ifelse(df$padj < input$de_padj_cut & df$lfc >=  input$de_lfc_cut, "Up",
                 ifelse(df$padj < input$de_padj_cut & df$lfc <= -input$de_lfc_cut, "Down", "NS"))
    df
  })
  
  output$ui_de_summary <- renderUI({
    req(de_data())
    df <- de_data()
    helpText(paste0("✅  ", sum(df$status == "Up"), " up  /  ", sum(df$status == "Down"),
                    " down  of ", nrow(df), " tested features"))
  })
  
  # Significant features (genes / proteins / metabolites), ordered by padj
  de_sig_genes <- function(direction = "both") {
    df  <- de_data()
    keep <- switch(direction, up = "Up", down = "Down", c("Up", "Down"))
    df  <- df[df$status %in% keep & !is.na(df$gene) & df$gene != "", , drop = FALSE]
    unique(df$gene[order(df$padj)])
  }
  
  # ── Volcano ──────────────────────────────────────────────────────────────
  output$volcano_status <- renderText({
    req(de_data())
    df <- de_data()
    paste0(
      "Feature ID column: ", input$de_col_gene, "\n",
      "Fold-change col  : ", input$de_col_lfc,  "\n",
      "Adj. p-value col : ", input$de_col_padj, "\n",
      "Cutoffs          : |log2FC| >= ", input$de_lfc_cut, " and padj < ", input$de_padj_cut, "\n",
      "Features tested  : ", nrow(df), " (rows with missing padj removed: ", nrow(de_raw()) - nrow(df), ")\n",
      "Up / Down        : ", sum(df$status == "Up"), " / ", sum(df$status == "Down")
    )
  })
  
  make_volcano <- function() {
    req(de_data())
    df <- de_data()
    up_col   <- if (is.null(input$volcano_up_color)   || input$volcano_up_color   == "") "firebrick" else input$volcano_up_color
    down_col <- if (is.null(input$volcano_down_color) || input$volcano_down_color == "") "steelblue" else input$volcano_down_color
    
    # Labels: top N significant by padj + any user-typed genes (case-insensitive)
    sig      <- df[df$status != "NS", , drop = FALSE]
    top_lab  <- head(sig$gene[order(sig$padj)], max(0, input$volcano_top_labels))
    typed_txt <- if (is.null(input$volcano_label_genes)) "" else input$volcano_label_genes
    typed    <- tolower(trimws(strsplit(typed_txt, ",", fixed = TRUE)[[1]]))
    typed    <- typed[nchar(typed) > 0]
    lab_df   <- df[df$gene %in% top_lab | tolower(df$gene) %in% typed, , drop = FALSE]
    lab_df   <- lab_df[!is.na(lab_df$gene) & lab_df$gene != "" & !duplicated(lab_df$gene), , drop = FALSE]
    
    n_up <- sum(df$status == "Up"); n_down <- sum(df$status == "Down")
    df$status <- factor(df$status, levels = c("Down", "NS", "Up"))
    
    ggplot(df, aes(x = lfc, y = neg_log10_padj, colour = status)) +
      geom_point(size = 1.4, alpha = 0.7) +
      geom_vline(xintercept = c(-input$de_lfc_cut, input$de_lfc_cut), linetype = "dashed", colour = "grey40") +
      geom_hline(yintercept = -log10(input$de_padj_cut), linetype = "dashed", colour = "grey40") +
      geom_text_repel(data = lab_df, aes(label = gene), colour = "black", size = 3.2,
                      max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4) +
      scale_colour_manual(values = c(Down = down_col, NS = "grey75", Up = up_col),
                          labels = c(Down = paste0("Down (", n_down, ")"), NS = "Not significant",
                                     Up = paste0("Up (", n_up, ")")),
                          drop = FALSE, name = NULL) +
      labs(title = "Volcano Plot",
           x = paste0("log2 fold change (", input$de_col_lfc, ")"),
           y = paste0("-log10(", input$de_col_padj, ")")) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 14), legend.position = "top")
  }
  
  output$volcano_plot <- renderPlot({ make_volcano() }, height = function() input$volcano_height)
  dl_volcano_handler <- function(file) dl_png(file, "volcano", function() print(make_volcano()))
  dl_volcano_filename <- function() paste0("Volcano_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  output$download_volcano      <- downloadHandler(filename = dl_volcano_filename, content = dl_volcano_handler)
  output$download_volcano_main <- downloadHandler(filename = dl_volcano_filename, content = dl_volcano_handler)
  
  # ── PCA ──────────────────────────────────────────────────────────────────
  pca_result <- reactive({
    req(norm_prot(), annotation(), input$col_gene, input$col_sample_id, input$col_group,
        sample_cols_resolved(), input$pca_gene_set)
    data        <- norm_prot()
    sample_cols <- sample_cols_resolved()
    validate(need(length(sample_cols) >= 3, "PCA needs at least 3 samples after metadata filtering."))
    
    mat <- as.matrix(data[, sample_cols, drop = FALSE])
    rownames(mat) <- as.character(data[[input$col_gene]])
    mat <- mat[complete.cases(mat), , drop = FALSE]
    if (isTRUE(input$de_log_transform)) {
      validate(need(all(mat >= 0),
                    "Data contain negative values (already log-transformed?). Untick 'Apply log2(x + 1)'."))
      mat <- log2(mat + 1)
    }
    vars <- apply(mat, 1, var)
    mat  <- mat[vars > 0, , drop = FALSE]
    vars <- vars[vars > 0]
    
    if (input$pca_gene_set == "top") {
      req(input$pca_top_n)
      mat <- mat[order(vars, decreasing = TRUE)[seq_len(min(input$pca_top_n, nrow(mat)))], , drop = FALSE]
    } else if (input$pca_gene_set == "de") {
      validate(need(!is.null(input$file_de), "Upload the differential analysis file to use DE-significant features."))
      mat <- mat[rownames(mat) %in% de_sig_genes("both"), , drop = FALSE]
    }
    validate(need(nrow(mat) >= 2, "Fewer than 2 features available for PCA with the current settings."))
    
    pca <- prcomp(t(mat), center = TRUE, scale. = isTRUE(input$pca_scale))
    ann <- annotation()
    grp <- as.character(ann[[input$col_group]])[match(rownames(pca$x), as.character(ann[[input$col_sample_id]]))]
    list(pca = pca, n_genes = nrow(mat), group = grp,
         var_pct = 100 * pca$sdev^2 / sum(pca$sdev^2))
  })
  
  output$pca_status <- renderText({
    req(pca_result())
    r <- pca_result()
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter  : ", input$filter_col, " = ", paste(input$filter_vals, collapse = ", "), "\n")
    else ""
    set_label <- switch(input$pca_gene_set,
                        all = "all features", top = paste0("top ", input$pca_top_n, " most variable"),
                        de  = "DE-significant features")
    paste0(
      filter_line,
      "Samples          : ", nrow(r$pca$x), "\n",
      "Features used    : ", r$n_genes, " (", set_label, ")\n",
      "Transformation   : ", if (isTRUE(input$de_log_transform)) "log2(x + 1)" else "none",
      if (isTRUE(input$pca_scale)) ", scaled to unit variance" else ", centered only", "\n",
      "Variance (PC1-5) : ", paste0(round(head(r$var_pct, 5), 1), "%", collapse = ", ")
    )
  })
  
  make_pca <- function() {
    req(pca_result(), boxplot_group_colors(), input$pca_x, input$pca_y)
    r <- pca_result()
    n_pc <- ncol(r$pca$x)
    validate(need(input$pca_x <= n_pc && input$pca_y <= n_pc && input$pca_x != input$pca_y,
                  paste0("Choose two different PCs between 1 and ", n_pc, ".")))
    df <- data.frame(x = r$pca$x[, input$pca_x], y = r$pca$x[, input$pca_y],
                     sample = rownames(r$pca$x), group = r$group, stringsAsFactors = FALSE)
    grp_colors <- boxplot_group_colors()
    filter_subtitle <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0(input$filter_col, " = ", paste(input$filter_vals, collapse = ", ")) else NULL
    
    p <- ggplot(df, aes(x = x, y = y, colour = group)) +
      geom_point(size = 3.5, alpha = 0.9) +
      scale_colour_manual(values = grp_colors, name = input$col_group) +
      labs(title    = "PCA of Normalized Data",
           subtitle = filter_subtitle,
           x = paste0("PC", input$pca_x, " (", round(r$var_pct[input$pca_x], 1), "%)"),
           y = paste0("PC", input$pca_y, " (", round(r$var_pct[input$pca_y], 1), "%)")) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold", size = 14))
    if (isTRUE(input$pca_show_labels))
      p <- p + geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE, max.overlaps = Inf)
    p
  }
  
  output$pca_plot <- renderPlot({ make_pca() }, height = function() input$pca_height)
  dl_pca_handler  <- function(file) dl_png(file, "pca", function() print(make_pca()))
  dl_pca_filename <- function() paste0("PCA_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  output$download_pca      <- downloadHandler(filename = dl_pca_filename, content = dl_pca_handler)
  output$download_pca_main <- downloadHandler(filename = dl_pca_filename, content = dl_pca_handler)
  
  # ── DE heatmap (reuses make_heatmap) ─────────────────────────────────────
  de_heatmap_genes <- reactive({
    req(de_data(), norm_prot(), input$col_gene, input$de_heatmap_direction, input$de_heatmap_top_n)
    genes <- de_sig_genes(input$de_heatmap_direction)
    genes <- genes[genes %in% as.character(norm_prot()[[input$col_gene]])]
    head(genes, input$de_heatmap_top_n)
  })
  
  output$de_heatmap_status <- renderText({
    req(de_data(), de_heatmap_genes(), sample_cols_resolved())
    n_sig <- length(de_sig_genes(input$de_heatmap_direction))
    filter_line <- if (!is.null(input$filter_col) && input$filter_col != "__none__")
      paste0("Metadata filter  : ", input$filter_col, " = ", paste(input$filter_vals, collapse = ", "), "\n")
    else ""
    paste0(
      filter_line,
      "Cutoffs          : |log2FC| >= ", input$de_lfc_cut, " and padj < ", input$de_padj_cut, "\n",
      "Direction        : ", input$de_heatmap_direction, "\n",
      "Significant      : ", n_sig, "\n",
      "Shown (in data)  : ", length(de_heatmap_genes()), " (max ", input$de_heatmap_top_n, ", by padj)\n",
      "Samples used     : ", length(sample_cols_resolved()), "\n",
      "Transformation   : ", if (isTRUE(input$de_log_transform)) "log2(x + 1), then row Z-score" else "row Z-score"
    )
  })
  
  de_heatmap_title <- function()
    paste0("DE features: |log2FC| ≥ ", input$de_lfc_cut, ", padj < ", input$de_padj_cut,
           "\n(Z-score, by Group)")
  
  draw_de_heatmap <- function() {
    genes <- de_heatmap_genes()
    validate(need(length(genes) >= 2,
                  "Fewer than 2 significant features found in the normalized data. Relax the cutoffs or check that the DE feature ID column matches the normalized data ID column."))
    make_heatmap(gene_filter = genes, title_text = de_heatmap_title(),
                 log2_transform = isTRUE(input$de_log_transform))
  }
  
  output$de_heatmap_plot <- renderPlot({ draw_de_heatmap() }, height = function() input$de_heatmap_height)
  dl_de_heatmap_handler <- function(file) {
    dl_png(file, "deheat", draw_de_heatmap)
  }
  dl_de_heatmap_filename <- function() paste0("DE_Heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
  output$download_de_heatmap      <- downloadHandler(filename = dl_de_heatmap_filename, content = dl_de_heatmap_handler)
  output$download_de_heatmap_main <- downloadHandler(filename = dl_de_heatmap_filename, content = dl_de_heatmap_handler)
}

# ══════════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)