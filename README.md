# EnrichViz Shiny App — v1.1.0

Interactive visualization of enriched pathways or biological functions from
normalized proteomics or transcriptomics data.

Copyright 2026 RGM  
Licensed under the MIT License — see [License](#license) section below.

---

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Input Files](#input-files)
5. [App Structure](#app-structure)
6. [Tabs and Features](#tabs-and-features)
   - [Tab 1 -- Bar Plot](#tab-1--bar-plot)
   - [Tab 2 -- Chord Diagram](#tab-2--chord-diagram)
   - [Tab 3 -- Heatmap](#tab-3--heatmap)
   - [Tab 4 -- Boxplot](#tab-4--boxplot)
7. [Sidebar Controls Reference](#sidebar-controls-reference)
8. [Download Outputs](#download-outputs)
9. [Deployment](#deployment)
10. [Tips and Troubleshooting](#tips-and-troubleshooting)
11. [License](#license)

---

## Overview

EnrichViz (v1.1.0) is an R Shiny application for interactive exploration and
visualization of pathway or functional enrichment results derived from
normalized proteomics or transcriptomics experiments.

The app accepts two plain-text input files (a data matrix and an annotation
table) and produces four coordinated views:

- A **Bar Plot** (optionally rendered as a bubble chart) ranking pathways or
  functions by a user-selected metric.
- A **Chord Diagram** showing co-membership of molecules across pathways or
  functions.
- A **Heatmap** of molecule-level abundances grouped by pathway or function.
- A **Boxplot** (optionally rendered as a violin plot) displaying the
  distribution of abundance values for a selected molecule across groups.

All plots can be downloaded as publication-ready PDF or PNG files directly
from the sidebar.

The browser window title is set to:
*"Heatmap, Bar Plot and Chord Diagram of Normalized Data by Pathway or
Function"*

---

## Requirements

| Dependency | Minimum version |
|---|---|
| R | >= 4.2.0 |
| shiny | >= 1.7.4 |
| ggplot2 | >= 3.4.0 |
| dplyr | >= 1.1.0 |
| tidyr | >= 1.3.0 |
| readr | >= 2.1.0 |
| stringr | >= 1.5.0 |
| circlize | >= 0.4.15 |
| ComplexHeatmap | >= 2.14.0 |
| viridis | >= 0.6.3 |
| RColorBrewer | >= 1.1-3 |
| scales | >= 1.2.1 |
| ggforce | >= 0.4.1 |
| shinycssloaders | >= 1.0.0 |

Install all R dependencies at once:

```r
install.packages(c(
  "shiny", "ggplot2", "dplyr", "tidyr", "readr",
  "stringr", "circlize", "viridis", "RColorBrewer",
  "scales", "ggforce", "shinycssloaders"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

---

## Installation

1. Clone or download this repository.
2. Place `app.R` (v1.1.0) in a dedicated project folder.
3. Install all dependencies listed above.
4. Launch the app from R or RStudio:

```r
shiny::runApp("path/to/project/folder")
```

Or open `app.R` in RStudio and click **Run App**.

---

## Input Files

EnrichViz expects two tab-delimited plain-text files (`.txt` or `.tsv`).

### 1. Data Matrix File

A numeric abundance matrix where:

- **Row 1** is a header row.
- **Column 1** contains molecule identifiers (protein or gene IDs / symbols).
- **Remaining columns** are sample abundance values (log-transformed or
  normalized counts).

Example layout:

```
Protein   Sample1  Sample2  Sample3  Sample4
ProtA     12.3     11.8     13.1     12.6
ProtB      9.4      9.9      8.7      9.2
...
```

### 2. Annotation File

A mapping table where:

- **Column 1** contains pathway or function names.
- **Column 2** contains molecule identifiers matching column 1 of the data
  matrix.
- Additional columns are ignored.

Example layout:

```
Pathway            Molecule
Glycolysis         ProtA
Glycolysis         ProtC
Oxidative_Phospho  ProtB
Oxidative_Phospho  ProtA
...
```

> **Note:** Molecule identifiers must match exactly (case-sensitive) between
> the two files.

---

## App Structure

The app uses a standard Shiny `fluidPage` layout with a fixed left sidebar
(`sidebarPanel`) and a main content area (`mainPanel`).

The main panel contains a `tabsetPanel` with four tabs:

| Tab index | Tab label |
|---|---|
| 1 | Bar Plot |
| 2 | Chord Diagram |
| 3 | Heatmap |
| 4 | Boxplot |

Each tab renders its plot inside a `withSpinner` loading indicator.

The sidebar is divided into clearly labelled collapsible sections (implemented
as `conditionalPanel` or plain `wellPanel` blocks):

- **File Inputs** -- upload the two required data files.
- **Global Filters** -- select pathways/functions and sample groups.
- **Bar Plot Settings** -- controls specific to the Bar Plot tab.
- **Chord Diagram Settings** -- controls specific to the Chord Diagram tab.
- **Heatmap Settings** -- controls specific to the Heatmap tab.
- **Boxplot / Violin Settings** -- controls specific to the Boxplot tab.
- **Download** -- buttons to save the currently visible plot.

---

## Tabs and Features

### Tab 1 -- Bar Plot

Displays pathways or functions ranked by a chosen metric (e.g., molecule
count, mean abundance, enrichment score).

| Feature | Details |
|---|---|
| Plot type toggle | Switch between a classic horizontal bar chart and a bubble chart (`checkboxInput`) |
| Metric selector | Choose the y-axis / size metric from a `selectInput` |
| Colour palette | Pick a discrete colour palette (`selectInput`) |
| Number of terms | Slider to cap the number of displayed pathways (`sliderInput`) |
| Sort order | Ascending or descending sort (`radioButtons`) |
| Status line | Displays: selected metric, number of terms shown, and current sort order |

The status line (`bar_status`) renders text in the format:

```
Metric           : <metric>
Terms shown      : <n>
Sort             : <order>
```

### Tab 2 -- Chord Diagram

Renders an interactive chord diagram (via `circlize`) showing which molecules
are shared across pathways or functions.

| Feature | Details |
|---|---|
| Minimum frequency filter | Show only molecules appearing in at least N pathways (`numericInput`, default 2) |
| Pathway colour palette | Colour sectors by pathway (`selectInput`) |
| Label size -- pathways | Font size for pathway sector labels (`numericInput`, default 0.55) |
| Label size -- molecules | Font size for molecule labels (`numericInput`, default 0.70) |
| Gap between sectors | Angular gap between chord sectors in degrees (`numericInput`) |
| Highlight top N molecules | Emphasise the N most-connected molecules (`numericInput`) |

### Tab 3 -- Heatmap

Renders a clustered heatmap (via `ComplexHeatmap`) of molecule abundances
within selected pathways or functions.

| Feature | Details |
|---|---|
| Colour scheme | Sequential or diverging palette (`selectInput`) |
| Row clustering | Toggle hierarchical clustering of rows (`checkboxInput`) |
| Column clustering | Toggle hierarchical clustering of columns (`checkboxInput`) |
| Show row names | Display molecule identifiers on the heatmap rows (`checkboxInput`) |
| Show column names | Display sample names on the heatmap columns (`checkboxInput`) |
| Row name font size | Numeric size for row labels (`numericInput`) |
| Column name font size | Numeric size for column labels (`numericInput`) |
| Scale rows | Z-score scale each row before plotting (`checkboxInput`) |

### Tab 4 -- Boxplot

Displays the distribution of abundance values for a single selected molecule
across sample groups, rendered as a box plot or violin plot.

| Feature | Details |
|---|---|
| Protein / molecule selector | `selectizeInput` with type-ahead search capability; supports up to 2000 options. Begin typing to filter the list. |
| Plot type toggle | Switch between boxplot and violin plot (`checkboxInput` or `radioButtons`) |
| Show data points | Overlay individual data points (jitter) on the plot (`checkboxInput`) |
| Colour by group | Colour boxes/violins by sample group (`checkboxInput`) |
| Status line | Displays summary statistics for the selected molecule |

The status line (`boxplot_status`) renders text in the format:

```
Molecule         : <id>
Groups           : <n>
Observations     : <n>
Values range     : <min> -- <max>
```

---

## Sidebar Controls Reference

Complete reference of every input widget present in the sidebar, grouped by
section.

### File Input Controls

| Input ID | Widget type | Label | Notes |
|---|---|---|---|
| `data_file` | `fileInput` | Data matrix file | Accepts `.txt` / `.tsv` |
| `annot_file` | `fileInput` | Annotation file | Accepts `.txt` / `.tsv` |

### Global Filter Controls

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `selected_pathways` | `selectizeInput` | Select pathways / functions | All selected |
| `selected_groups` | `checkboxGroupInput` | Select sample groups | All selected |

### Bar Plot Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `bar_metric` | `selectInput` | Metric | `"count"` |
| `bar_bubble` | `checkboxInput` | Show as bubble chart | `FALSE` |
| `bar_palette` | `selectInput` | Colour palette | `"Set2"` |
| `bar_n_terms` | `sliderInput` | Number of terms | `20` |
| `bar_sort_desc` | `checkboxInput` | Sort descending | `TRUE` |

### Chord Diagram Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `chord_min_freq` | `numericInput` | Minimum molecule frequency | `2` |
| `chord_palette` | `selectInput` | Sector colour palette | `"Set1"` |
| `chord_label_size_path` | `numericInput` | Label size -- pathways | `0.55` |
| `chord_label_size_mol` | `numericInput` | Label size -- molecules | `0.70` |
| `chord_gap` | `numericInput` | Gap between sectors (degrees) | `2` |
| `chord_top_n` | `numericInput` | Highlight top N molecules | `10` |

### Heatmap Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `hm_palette` | `selectInput` | Colour scheme | `"viridis"` |
| `hm_cluster_rows` | `checkboxInput` | Cluster rows | `TRUE` |
| `hm_cluster_cols` | `checkboxInput` | Cluster columns | `TRUE` |
| `hm_show_rownames` | `checkboxInput` | Show row names | `TRUE` |
| `hm_show_colnames` | `checkboxInput` | Show column names | `TRUE` |
| `hm_row_fontsize` | `numericInput` | Row name font size | `8` |
| `hm_col_fontsize` | `numericInput` | Column name font size | `8` |
| `hm_scale_rows` | `checkboxInput` | Scale rows (Z-score) | `FALSE` |

### Boxplot / Violin Settings

| Input ID | Widget type | Label | Default / Notes |
|---|---|---|---|
| `box_protein` | `selectizeInput` | Select protein / molecule | First in list; type-ahead search, up to 2000 options |
| `box_violin` | `checkboxInput` | Show as violin plot | `FALSE` |
| `box_points` | `checkboxInput` | Overlay data points | `TRUE` |
| `box_colour_group` | `checkboxInput` | Colour by group | `TRUE` |

### Download Controls

| Input ID | Widget type | Label | Notes |
|---|---|---|---|
| `dl_format` | `radioButtons` | File format | `"PDF"` / `"PNG"` |
| `dl_width` | `numericInput` | Width (inches) | `8` |
| `dl_height` | `numericInput` | Height (inches) | `6` |
| `dl_res` | `numericInput` | Resolution (DPI, PNG only) | `300` |
| `download_plot` | `downloadButton` | Download current plot | -- |

---

## Download Outputs

Click **Download current plot** in the sidebar to save the plot currently
visible in the active tab.

- **PDF** output uses the default system font and vector graphics; recommended
  for publication figures.
- **PNG** output uses the DPI value set in the sidebar (default 300).

File names are auto-generated in the format:

```
EnrichViz_<tab>_<YYYY-MM-DD>.pdf
EnrichViz_<tab>_<YYYY-MM-DD>.png
```

where `<tab>` is one of `BarPlot`, `ChordDiagram`, `Heatmap`, or `Boxplot`.

---

## Deployment

### Local (development)

```r
shiny::runApp("path/to/project")
```

### shinyapps.io

```r
library(rsconnect)
rsconnect::deployApp("path/to/project")
```

### Shiny Server / Posit Connect

Copy the project folder to the server apps directory and follow the standard
server deployment documentation.

---

## Tips and Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Chord diagram is empty | No molecules exceed the minimum frequency threshold | Lower `chord_min_freq` (default is 2) |
| Heatmap is blank | Selected pathways contain no overlapping molecules with the data matrix | Check identifier matching between input files |
| Boxplot shows a flat line | Only one observation per group | Verify group assignments in the annotation file |
| App crashes on upload | Wrong file delimiter or encoding | Ensure files are tab-delimited UTF-8 plain text |
| Molecule not found in selector | Identifier mismatch between files | Check for trailing spaces or case differences |
| Download produces empty PDF | Plot has not rendered yet | Switch to the desired tab and wait for the plot to appear before downloading |

---

## License

MIT License

Copyright 2026 RGM

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
