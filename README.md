# EnrichViz Shiny App — v1.2.0

Interactive visualization of enriched pathways or biological functions from
normalized proteomics or transcriptomics data.

Copyright 2026 RGM  
Licensed under the MIT License — see [License](#license) section below.

---

## Table of Contents

1. [Overview](#overview)
2. [What's New in v1.2.0](#whats-new-in-v120)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Input Files](#input-files)
6. [App Structure](#app-structure)
7. [Tabs and Features](#tabs-and-features)
   - [Tab 1 -- Bar Plot](#tab-1--bar-plot)
   - [Tab 2 -- Chord Diagram](#tab-2--chord-diagram)
   - [Tab 3 -- Heatmap](#tab-3--heatmap)
   - [Tab 4 -- Boxplot](#tab-4--boxplot)
8. [Sidebar Controls Reference](#sidebar-controls-reference)
9. [Download Outputs](#download-outputs)
10. [Deployment](#deployment)
11. [Tips and Troubleshooting](#tips-and-troubleshooting)
12. [License](#license)

---

## Overview

EnrichViz (v1.2.0) is an R Shiny application for interactive exploration and
visualization of pathway or functional enrichment results derived from
normalized proteomics or transcriptomics experiments.

The app accepts three CSV input files (a normalized data matrix, a sample
metadata/annotation table, and an enrichment results table) and produces four
coordinated views:

- A **Bar Plot** (optionally rendered as a bubble chart) ranking pathways or
  functions by a user-selected metric.
- A **Chord Diagram** showing co-membership of molecules across pathways or
  functions.
- A **Heatmap** of molecule-level abundances grouped by pathway or function.
- A **Boxplot** (optionally rendered as a violin plot) displaying the
  distribution of abundance values for a selected molecule across groups.

All plots can be downloaded as publication-ready PNG files directly
from the sidebar or from within each tab.

The browser window title is set to:
*"Heatmap, Bar Plot and Chord Diagram of Normalized Data by Pathway or
Function"*

---

## What's New in v1.2.0

### Metadata Slice Filter

The most significant addition in this release is a **Filter Samples by
Metadata** panel, highlighted in blue in the sidebar. It allows users to
subset the sample space before any visualization is rendered — without
modifying the uploaded files.

#### How it works

1. Upload the metadata/annotation CSV file.
2. In the **Filter Samples by Metadata** panel, choose any column from the
   metadata file using the **Filter column** drop-down (e.g. `sex`, `tissue`,
   `batch`).
3. A checkbox group is automatically populated with all unique values found in
   that column. All values are selected by default (no filtering).
4. Deselect one or more values to restrict the analysis to the remaining
   samples (e.g. keep only `male`, or keep only `control` and `disease`).
5. Use **Select all** / **Clear all** buttons to quickly toggle the entire
   set.
6. A live summary badge immediately below the checkboxes confirms how many
   samples are retained out of the total (e.g. `✅ 24 / 48 samples retained
   — sex: male`).
7. Choose `(none)` in the filter column drop-down to deactivate all filtering
   and restore the full sample set.

#### What is affected

The filter propagates automatically to **all four plots**:

| Plot | Effect |
|---|---|
| Bar Plot | Status line reports the active filter |
| Chord Diagram | Status line reports the active filter |
| Heatmap | Only filtered samples appear as columns; filter is shown in the plot title |
| Boxplot / Violin | Only filtered samples contribute data points; filter is shown in the plot subtitle |

#### Typical use case

A metadata file contains columns `sex` (`male` / `female`) and `diagnosis`
(`disease` / `healthy`). To compare disease vs. healthy **within males only**:

1. Set **Filter column** → `sex`.
2. Uncheck `female`, keep `male`.
3. Set **Group column** → `diagnosis`.
4. All plots now show only male samples split by diagnosis.

To then repeat the comparison for females, change the checked value to
`female` — no file re-upload is needed.

#### Implementation notes

- The raw uploaded annotation is stored in an internal `annotation_raw()`
  reactive and is never modified.
- The public `annotation()` reactive returns the filtered subset and is the
  sole source of sample metadata consumed by all plots and helper reactives
  (including `sample_cols_resolved()` and `boxplot_group_colors()`).
- String representations of missing values (`"NA"`, `"NaN"`, `""`) are
  excluded from the checkbox list automatically.
- The filter panel uses `checkboxGroupInput` for categorical columns.
  **Select all** and **Clear all** are `actionButton` controls that call
  `updateCheckboxGroupInput()` without triggering a file reload.

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
| pheatmap | >= 1.0.12 |
| viridis | >= 0.6.3 |
| RColorBrewer | >= 1.1-3 |
| scales | >= 1.2.1 |

Install all R dependencies at once:

```r
install.packages(c(
  "shiny", "ggplot2", "dplyr", "tidyr", "readr",
  "stringr", "circlize", "pheatmap",
  "viridis", "RColorBrewer", "scales"
))
```

> **Note:** v1.2.0 replaces `ComplexHeatmap` (Bioconductor) with `pheatmap`
> (CRAN), removing the Bioconductor dependency and simplifying installation.

---

## Installation

1. Clone or download this repository.
2. Place `app.R` (v1.2.0) in a dedicated project folder.
3. Install all dependencies listed above.
4. Launch the app from R or RStudio:

```r
shiny::runApp("path/to/project/folder")
```

Or open `app.R` in RStudio and click **Run App**.

---

## Input Files

EnrichViz expects **three CSV files** (`.csv`). All files must have a header
row. Maximum upload size is 500 MB per file.

### 1. Normalized Data File

A numeric abundance matrix where:

- **Column 1** contains molecule identifiers (protein or gene IDs / symbols).
- **Remaining columns** are sample abundance values (log-transformed or
  normalized counts). Column names must match the sample IDs in the metadata
  file.

Example layout:

```
gene_symbol  Sample1  Sample2  Sample3  Sample4
ProtA        12.3     11.8     13.1     12.6
ProtB         9.4      9.9      8.7      9.2
```

### 2. Metadata / Sample Annotation File

A sample-level table where:

- One column contains **sample IDs** matching the column names in the
  normalized data file.
- One column contains the **group or condition** label used for colouring
  plots.
- **Additional columns** (e.g. `sex`, `age`, `batch`, `tissue`) are fully
  supported and can be used as **metadata slice filters** (new in v1.2.0).

Example layout:

```
sample_id  diagnosis  sex     age
Sample1    disease    male    62
Sample2    healthy    female  58
Sample3    disease    female  71
Sample4    healthy    male    55
```

> The app auto-detects the sample ID column (column 1) and the group column
> (column 2) but both can be overridden via sidebar drop-downs.

### 3. Enriched Pathways / Functions File

A results table from a pathway or functional enrichment analysis where:

- **Column 1** contains pathway or function names.
- One column contains a **p-value or -log10(p-value)**; the app auto-detects
  the scale.
- One column contains **member molecules / genes** as a delimited string
  (comma, slash, semicolon, or pipe — selectable in the sidebar).

Example layout:

```
Pathway            pvalue   Molecules
Glycolysis         0.001    ProtA,ProtB,ProtC
Oxidative_Phospho  0.003    ProtB,ProtD
```

> **Note:** Molecule identifiers must match exactly (case-sensitive) the
> identifiers in the normalized data file.

---

## App Structure

The app uses a standard Shiny `fluidPage` layout with a fixed left sidebar
(`sidebarPanel`, width = 3) and a fixed main content area (`mainPanel`,
width = 9). Both panels scroll independently.

A fixed blurb bar below the title summarises the three required input files.

The main panel contains a `tabsetPanel` with four tabs:

| Tab index | Tab label |
|---|---|
| 1 | Bar Plot |
| 2 | Chord Diagram |
| 3 | Heatmap |
| 4 | Boxplot |

The sidebar is divided into clearly labelled sections:

- **Upload Data Files** — three `fileInput` widgets.
- **Filter Samples by Metadata** *(new in v1.2.0)* — column picker,
  value checkboxes, Select all / Clear all buttons, live sample count badge.
- **Select Pathway / Functions — Column** — category and molecule column
  selectors, separator radio buttons.
- **Normalized counts or protein abundance — Column** — gene/protein ID
  column selector.
- **Sample Annotation — Column** — sample ID and group column selectors.
- **Bar / Bubble Plot Settings** — plot type, p-value column, top N,
  colour / palette controls.
- **Chord Diagram Settings** — top N, minimum frequency, inner radius,
  label sizes, sector colours.
- **Heatmap Settings** — category selector, per-group colour inputs, height,
  download buttons.
- **Boxplot / Violin Settings** — molecule selector, plot type, height,
  download button.

---

## Tabs and Features

### Tab 1 -- Bar Plot

Displays pathways or functions ranked by a chosen p-value metric.

| Feature | Details |
|---|---|
| Plot type toggle | Switch between a horizontal bar chart and a bubble chart (`radioButtons`) |
| P-value column | Select the significance column from the enrichment file (`selectInput`); raw p-values (0–1) are auto-transformed with −log10 |
| Auto-transform indicator | `helpText` confirms whether −log10 will be applied |
| Top N selector | Cap the number of displayed pathways (`numericInput`, default 20) |
| Bar colour | Hex or named colour for bar fill (`textInput`, default `steelblue`) |
| Bubble size column | Column mapped to bubble area (`selectInput`) |
| Bubble colour column | Column mapped to bubble fill intensity (`selectInput`) |
| Bubble size range | Min/max point size (`sliderInput`, default 3–15) |
| Bubble palette | Sequential or diverging colour palette (`selectInput`) |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 600) |
| Status line | Reports active metadata filter, plot type, p-value column, transformation applied, total and displayed pathway counts |

### Tab 2 -- Chord Diagram

Renders a chord diagram (via `circlize`) showing molecule co-membership across
the top N pathways or functions.

| Feature | Details |
|---|---|
| Top N pathways | Number of highest-ranked pathways to include (`numericInput`, default 10) |
| Minimum frequency | Show only molecules appearing in at least N pathways (`numericInput`, default 2) |
| Inner circle size | Radius of the central hole (`sliderInput`, default 0.4) |
| Pathway label size | Font size for pathway sector labels (`numericInput`, default 0.55) |
| Molecule label size | Font size for molecule labels (`numericInput`, default 0.70) |
| Pathway sector colour | Hex or named colour for pathway arcs (`textInput`, default `tomato`) |
| Molecule sector colour | Hex or named colour for molecule arcs (`textInput`, default `steelblue`) |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 750) |
| Status line | Reports active metadata filter, pathway count, molecule count, total connections, minimum frequency, and inner radius |

### Tab 3 -- Heatmap

Renders a clustered heatmap (via `pheatmap`) of Z-score normalised molecule
abundances for a selected pathway or function, using only the filtered samples.

| Feature | Details |
|---|---|
| Category selector | Choose one pathway or function to visualise (`selectInput`) |
| Per-group colours | One `textInput` per group (hex or name) for the column annotation bar |
| Row clustering | Hierarchical clustering of molecules (always on) |
| Column clustering | Disabled — columns are ordered by group for comparability |
| Row name display | Shown when ≤ 100 molecules are matched |
| Z-score scaling | Applied row-wise before rendering |
| Plot title | Includes the selected category name and active metadata filter |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 700) |
| Status line | Reports active metadata filter, samples used, category name, gene list size, matched proteins, and sample column names |

### Tab 4 -- Boxplot

Displays the distribution of normalized abundance values for a single selected
molecule across sample groups, using only the filtered samples.

| Feature | Details |
|---|---|
| Molecule selector | `selectizeInput` with type-ahead search; lists all molecules present in both the enrichment file and the normalized data (up to 2000 options) |
| Plot type toggle | Switch between boxplot and violin plot (`radioButtons`) |
| Data points overlay | Individual sample points overlaid via `geom_jitter` (always on) |
| Colour by group | Per-group colours shared with the Heatmap Settings section |
| Plot subtitle | Shows the active metadata filter when one is set |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 500) |
| Status line | Reports active metadata filter, plot type, selected molecule, sample count, groups, and value range |

---

## Sidebar Controls Reference

### File Input Controls

| Input ID | Widget type | Label | Accepts |
|---|---|---|---|
| `file_norm_prot` | `fileInput` | Normalized Data (.csv) | `.csv` |
| `file_annotation` | `fileInput` | Metadata or Sample Annotation (.csv) | `.csv` |
| `file_ipa_funct` | `fileInput` | Enriched Pathway or Functions (.csv) | `.csv` |

### Metadata Filter Controls *(new in v1.2.0)*

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `filter_col` | `selectInput` | Filter column | `(none)` |
| `filter_vals` | `checkboxGroupInput` | Keep values in `<col>` | All values selected |
| `filter_select_all` | `actionButton` | Select all | — |
| `filter_select_none` | `actionButton` | Clear all | — |

### Pathway / Function Column Controls

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `col_category` | `selectInput` | Category column | First column of enrichment file |
| `col_molecules` | `selectInput` | Molecules column | Auto-detected by keyword |
| `mol_separator` | `radioButtons` | Gene/molecule separator | `,` |

### Normalized Data Column Controls

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `col_gene` | `selectInput` | Gene/Protein ID column | Auto-detected by keyword |

### Sample Annotation Column Controls

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `col_sample_id` | `selectInput` | Sample ID column | First column of annotation file |
| `col_group` | `selectInput` | Group column | Second column of annotation file |

### Bar / Bubble Plot Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `barplot_type` | `radioButtons` | Plot type | `bar` |
| `col_pvalue` | `selectInput` | P-value column | Auto-detected by keyword |
| `bar_top_n` | `numericInput` | Show top N pathways | `20` |
| `bar_color` | `textInput` | Bar fill colour | `steelblue` |
| `bubble_size_col` | `selectInput` | Bubble SIZE column | Auto-detected by keyword |
| `bubble_color_col` | `selectInput` | Bubble COLOUR column | Auto-detected by keyword |
| `bubble_size_range` | `sliderInput` | Bubble size range | `3` – `15` |
| `bubble_palette` | `selectInput` | Colour palette | `Blues` |
| `bar_plot_height` | `numericInput` | Plot height (px) | `600` |

### Chord Diagram Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `chord_top_n` | `numericInput` | Top N pathways | `10` |
| `chord_min_freq` | `numericInput` | Minimum molecule frequency | `2` |
| `chord_inner_radius` | `sliderInput` | Inner circle size | `0.4` |
| `chord_label_size_path` | `numericInput` | Pathway label font size | `0.55` |
| `chord_label_size_mol` | `numericInput` | Molecule label font size | `0.70` |
| `chord_pathway_color` | `textInput` | Pathway sector colour | `tomato` |
| `chord_molecule_color` | `textInput` | Molecule sector colour | `steelblue` |
| `chord_plot_height` | `numericInput` | Chord plot height (px) | `750` |

### Heatmap Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `selected_category` | `selectInput` | Select category | First category in enrichment file |
| `color_group_<i>` | `textInput` | Colour — `<group>` | Cycles through preset palette |
| `plot_height` | `numericInput` | Heatmap height (px) | `700` |

### Boxplot / Violin Settings

| Input ID | Widget type | Label | Default / Notes |
|---|---|---|---|
| `plot_type` | `radioButtons` | Plot type | `box` |
| `boxplot_protein` | `selectizeInput` | Search protein / gene | First matched molecule; type-ahead, up to 2000 options |
| `boxplot_height` | `numericInput` | Plot height (px) | `500` |

### Download Controls

Each tab exposes two download buttons (sidebar and in-tab footer) for its
respective plot:

| Output ID | Label | Format |
|---|---|---|
| `download_barplot` / `download_barplot_main` | Download Plot (.png) | PNG, 4200 × 2400 px, 300 DPI |
| `download_chord` / `download_chord_main` | Download Chord Diagram (.png) | PNG, 10 × 10 in, 300 DPI |
| `download_single` / `download_single_main` | Download Selected Heatmap (.png) | PNG, 8 × 10 in, 300 DPI |
| `download_all` / `download_all_main` | Download ALL Heatmaps (.zip) | ZIP of PNG files, one per category |
| `download_boxplot_sidebar` / `download_boxplot_main` | Download Plot (.png) | PNG, 8 × 6 in, 300 DPI |

File names are auto-generated in the format:

```
Barplot_Enriched_Pathways_YYYYMMDD_HHMMSS.png
Bubbleplot_Enriched_Pathways_YYYYMMDD_HHMMSS.png
Chord_Diagram_Pathways_Molecules.png
Heatmap_<Category>.png
All_Heatmaps.zip
Boxplot_<Molecule>.png
Violin_<Molecule>.png
```

---

## Deployment

### Local (development)

```r
shiny::runApp("path/to/project")
```

### shinyapps.io

```r
library(rsconnect)

rsconnect::setAccountInfo(
  name   = "your-account-name",
  token  = "YOUR_TOKEN",
  secret = "YOUR_SECRET"
)

rsconnect::deployApp(
  appDir  = "path/to/project",
  appName = "EnrichViz",
  account = "your-account-name"
)
```

> ⚠️ **Cairo graphics** are not available on the shinyapps.io free tier.
> Remove `type = "cairo-png"` from all `png()` calls in the download handlers
> before deploying.

### Shiny Server / Posit Connect

Copy the project folder to the server apps directory and follow the standard
server deployment documentation.

---

## Tips and Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Chord diagram is empty | No molecules exceed the minimum frequency threshold | Lower `chord_min_freq` (default is 2) |
| Heatmap is blank | Selected pathway contains no molecules matching the data matrix | Check identifier matching between input files |
| Heatmap shows fewer columns than expected | Metadata filter is active | Check the filter panel badge; broaden the selection or choose `(none)` |
| Boxplot shows a flat line or single point | Only one sample per group after filtering | Broaden the metadata filter or check group assignments |
| Filter checkbox group is empty | Selected filter column contains only `NA` / blank values | Choose a different filter column |
| App crashes on upload | Wrong file format or encoding | Ensure files are comma-delimited UTF-8 CSV |
| Molecule not found in boxplot selector | Identifier mismatch between files | Check for trailing spaces or case differences |
| Download produces a blank PNG | Plot has not rendered yet | Switch to the desired tab and wait for the plot before downloading |
| `cairo-png` error on shinyapps.io | Cairo not installed on server | Remove `type = "cairo-png"` from all `png()` calls |
| `object 'X' not found` on deploy | Unquoted string in `deployApp()` call | Wrap `appName` and `account` values in quotes |

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
