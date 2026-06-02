# EnrichViz Shiny App — v1.0.9

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
   - [Tab 1 — Bar / Bubble Plot](#tab-1--bar--bubble-plot)
   - [Tab 2 — Chord Diagram](#tab-2--chord-diagram)
   - [Tab 3 — Heatmap](#tab-3--heatmap)
   - [Tab 4 — Boxplot / Violin Plot](#tab-4--boxplot--violin-plot)
7. [Sidebar Controls Reference](#sidebar-controls-reference)
8. [Download Outputs](#download-outputs)
9. [Deployment](#deployment)
10. [Tips and Troubleshooting](#tips-and-troubleshooting)
11. [License](#license)

---

## Overview

EnrichViz is a self-contained R Shiny application that takes three CSV files
as input and produces four interactive, publication-ready visualizations:

| Visualization | What it shows |
|---|---|
| **Bar / Bubble Plot** | Top N enriched pathways ranked by significance, as a horizontal bar chart or a bubble chart with configurable size and colour intensity |
| **Chord Diagram** | Pathway–molecule connectivity and shared molecules |
| **Heatmap** | Z-score expression patterns per pathway/function |
| **Boxplot / Violin Plot** | Per-group normalized abundance for a single protein/gene, as a boxplot or violin plot |

All plots are adjustable in real time via the sidebar and can be downloaded
as high-resolution image files.

---

## Requirements

### R version
R ≥ 4.2.0 recommended.

### R packages

| Package | Version tested | Purpose |
|---|---|---|
| `shiny` | ≥ 1.8 | App framework |
| `tidyverse` | ≥ 2.0 | Data wrangling and ggplot2 plotting |
| `pheatmap` | ≥ 1.0.12 | Clustered heatmap rendering |
| `circlize` | ≥ 0.4.15 | Chord diagram rendering |

Install all dependencies at once:

```r
install.packages(c("shiny", "tidyverse", "pheatmap", "circlize"))
```

---

## Installation

1. Clone or download this repository.
2. Place `app.R` in a folder of your choice.
3. Open R or RStudio and run:

```r
shiny::runApp("path/to/your/folder")
```

Or open `app.R` in RStudio and click **Run App**.

> **Upload limit:** The app accepts files up to **500 MB** per upload.

---

## Input Files

The app requires **three CSV files**. All must have a header row.

---

### 1. Normalized Data file

Rows = proteins or genes identifiers. Columns = sample measurements plus an identifier column.

| Column type | Description |
|---|---|
| **Gene/Protein ID** | One column with gene symbols or protein IDs (e.g. `Gene.Symbol`). Selected in sidebar under *Normalized counts or protein abundance — Column*. |
| **Sample columns** | One column per sample. Column names must match the Sample ID values in the Annotation file. |

**Minimal example:**

```
Gene.Symbol, Sample_A1, Sample_A2, Sample_B1, Sample_B2
TP53,        12.4,      11.9,      8.1,        7.8
EGFR,        9.0,       9.3,       14.2,       13.9
```

---

### 2. Sample Annotation file

Rows = samples. Maps each sample to a group.

| Column type | Description |
|---|---|
| **Sample ID** | Must match the column names in the Normalized Data file |
| **Group** | Experimental group label (e.g. `Control`, `Treated`) |

**Minimal example:**

```
SampleID,  Group
Sample_A1, Control
Sample_A2, Control
Sample_B1, Treated
Sample_B2, Treated
```

---

### 3. Enriched Pathways or Functions file

Output from an enrichment tool such as Ingenuity Pathway Analysis, g:Profiler, Enrichr, or similar.
Rows = pathways or functions.

| Column type | Description |
|---|---|
| **Category** | Pathway or function name. Defaults to column 1. |
| **p-value / -log10(p)** | Significance column. Auto-detected by keyword matching (`fdr`, `adj`, `pval`, `-log`, etc.). |
| **Molecules** | Comma-separated list of gene/protein IDs in that pathway. Defaults to the last column. |

**Minimal example:**

```
Pathway,             p-value, Molecules
Cell cycle,          0.001,   TP53,EGFR,CDK2
Apoptosis,           0.005,   TP53,BCL2
DNA repair,          0.020,   BRCA1,TP53
```

> **p-value auto-detection:** The app checks whether values are between 0 and 1
> (raw p-values) or larger (already -log10 transformed) and applies the
> appropriate transformation automatically. A confirmation message is shown in
> the Bar Plot status line.

---

## App Structure

```
app.R
│
├── UI
│   ├── Sidebar (width = 3)
│   │   ├── File uploads
│   │   ├── Column mapping selectors
│   │   ├── Bar Plot Settings
│   │   ├── Chord Diagram Settings
│   │   ├── Heatmap Settings
│   │   └── Boxplot Settings
│   │
│   └── Main panel (width = 9)
│       ├── Tab 1 — Bar Plot
│       ├── Tab 2 — Chord Diagram
│       ├── Tab 3 — Heatmap
│       └── Tab 4 — Boxplot
│
└── Server
    ├── Reactive data loaders
    ├── Dynamic column-selector widgets
    ├── p-value transformation detection
    ├── Group colour inputs
    ├── Bar plot logic
    ├── Chord diagram logic
    ├── Heatmap logic
    └── Boxplot logic
```

---

## Tabs and Features

### Tab 1 — Bar / Bubble Plot

Displays the top N enriched pathways ranked by significance. Switch between
a classic **horizontal bar chart** and an enhanced **bubble chart** using the
plot type toggle in the sidebar.

#### Bar plot mode

| Feature | Detail |
|---|---|
| Ranking | Sorted descending by -log10(p-value) |
| p-value handling | Raw p-values are -log10 transformed automatically; pre-transformed values are used as-is |
| Top N | Configurable from 5 to 200 (default 20) |
| Bar colour | Any R colour name or hex code |
| Status line | Shows plot type, p-value column used, transformation applied, total pathways, and number displayed |
| Download | 300 dpi PNG, 10 × 8 inches, filename: `Barplot_Enriched_Pathways.png` |

#### Bubble plot mode

Bubble plot mode replaces bars with circles whose **size** and **fill colour
intensity** each encode a separate column from the enrichment file, allowing
three variables to be displayed simultaneously: significance (X-axis),
a quantity such as gene count (bubble size), and a secondary metric such as
FDR or fold enrichment (colour intensity).

| Feature | Detail |
|---|---|
| X-axis | Same -log10(p-value) ranking as bar plot mode |
| Bubble size | Maps to any numeric column selected from the enrichment file — auto-detects columns whose name contains keywords such as `count`, `size`, `gene`, `overlap`, `ratio` |
| Bubble colour intensity | Maps to any numeric column selected from the enrichment file — auto-detects FDR/p-value columns by the same keyword logic as the X-axis column |
| Colour auto-transform | If the colour column contains values between 0 and 1, -log10 is applied automatically; the legend label reflects the transformation |
| Size range | Slider to set minimum and maximum bubble diameter (default 3–15) |
| Colour palette | Choose from 10 palettes: Blues, Reds, Purples, Greens, OrRd, YlOrRd, RdYlBu, viridis, magma, plasma |
| Status line | Shows plot type, p-value column, transformation, total pathways, number displayed, bubble size column, and bubble colour column |
| Download | 300 dpi PNG, 10 × 8 inches, filename: `Bubbleplot_Enriched_Pathways.png` |

> **Choosing bubble columns:** Any numeric column from the enrichment file can
> be used for size or colour — they do not have to be the same column as the
> X-axis. The two bubble column selectors are independent, so you can, for
> example, use gene count for size and FDR for colour intensity.

---

### Tab 2 — Chord Diagram

Visualizes the connections between pathways and the molecules they contain.

| Feature | Detail |
|---|---|
| Pathways shown | Top N by significance (default 10) |
| Molecule filter | Minimum number of times a molecule must appear across pathways to be included (default 1) |
| Sector arc width | Proportional to number of connections — wider arc = more molecules (for pathways) or more pathways (for molecules) |
| Ribbon | One ribbon per pathway–molecule connection |
| Sector colours | Pathway family = tomato/red gradient; Molecule family = steelblue gradient. Shade within each family is alphabetical order only and carries no biological meaning |
| Inner radius | Adjustable via slider (0.1–0.7); smaller = more label space |
| Label font sizes | Separately configurable for pathways and molecules |
| Status line | Shows counts of pathways, molecules, and total connections |
| Download | 300 dpi PNG, 10 × 10 inches |

> **Reading the chord diagram:**
> - **Arc width** encodes connectivity — a wide pathway arc means many molecules;
>   a wide molecule arc means that molecule is shared across many pathways (hub protein).
> - **Colour family** (red vs. blue) distinguishes pathways from molecules.
> - **Colour shade** within a family is alphabetical only — it does not encode
>   significance, fold change, or any other variable.

---

### Tab 3 — Heatmap

Shows the Z-score normalized expression pattern for all proteins in a selected
pathway, with samples colour-annotated by group.

| Feature | Detail |
|---|---|
| Category selector | Dropdown populated from all unique pathway/function names in the enrichment file |
| Scaling | Row-wise Z-score (mean = 0, SD = 1 per protein) |
| Clustering | Rows (proteins) clustered; columns (samples) not clustered |
| Colour scale | Dark blue → white → dark red |
| Row labels | Shown when ≤ 100 proteins are displayed |
| Group annotation bar | One colour per group, configurable |
| Status line | Shows category name, gene list size, matched proteins, and sample columns used |
| Download (single) | Selected heatmap — 300 dpi PNG, 8 × 10 inches |
| Download (all) | All categories with ≥ 2 matched proteins — ZIP archive of PNGs |

---

### Tab 4 — Boxplot / Violin Plot

Shows the raw normalized abundance distribution for a single selected
protein/gene, one geometry per experimental group. Switch between a classic
**boxplot** and a **violin plot** using the plot type toggle in the sidebar.

#### Boxplot mode

| Feature | Detail |
|---|---|
| Protein selector | Dropdown populated from all unique molecules across all pathways in the enrichment file |
| Data source | Normalized Data file — same file used by the heatmap |
| Y-axis | Raw normalized values (no Z-score scaling) |
| Points | Individual samples shown as jittered dots overlaid on each box |
| Outlier display | Outlier points suppressed from the box geometry to avoid double-plotting with jitter |
| Group colours | Shared with the Heatmap Settings colour inputs |
| Status line | Shows plot type, selected protein, number of samples found, and group names |
| Download | 300 dpi PNG, 8 × 6 inches, filename: `Boxplot_<protein>.png` |

#### Violin plot mode

Violin plot mode replaces boxes with mirrored kernel density shapes, giving a
clearer view of the full value distribution within each group — particularly
useful when groups contain many samples or when the distribution is multimodal.

| Feature | Detail |
|---|---|
| Geometry | `geom_violin()` with `trim = FALSE` to show full distribution tails |
| Points | Individual samples shown as jittered dots overlaid on each violin |
| Group colours | Shared with the Heatmap Settings colour inputs — same colours as boxplot mode |
| Status line | Shows plot type, selected protein, number of samples found, and group names |
| Download | 300 dpi PNG, 8 × 6 inches, filename: `Violin_<protein>.png` |

> **Choosing between boxplot and violin plot:**
> - Use **Boxplot** when sample sizes are small (< ~10 per group) and you want
>   to clearly see the median, IQR, and whiskers.
> - Use **Violin plot** when sample sizes are larger and you want to visualize
>   the full shape of the distribution, including bimodality or skewness.

---

## Sidebar Controls Reference

### File Uploads

| Input | Description |
|---|---|
| Normalized Data (.csv) | Protein/gene abundance matrix |
| Metadata or Sample Annotation (.csv) | Sample-to-group mapping |
| Enriched Pathway or Functions (.csv) | Enrichment results with molecule lists |

### Column Mapping

| Input | Description |
|---|---|
| Category column | Pathway/function name column in the enrichment file |
| Molecules column | Molecule list column — auto-detects by keyword; defaults to last column |
| Gene/molecule separator | Delimiter used between molecules in the molecules column: comma, slash, or semicolon |
| Gene/Protein ID column | Identifier column in the normalized data file |
| Sample ID column | Sample name column in the annotation file |
| Group column | Group/condition column in the annotation file |

### Bar / Bubble Plot Settings

| Input | Default | Description |
|---|---|---|
| **Plot type** | Bar plot | Toggle between Bar plot and Bubble plot |
| P-value column | Auto-detected | Column used for X-axis ranking (all modes) |
| Show top N pathways | 20 | Number of pathways to display (all modes) |
| Bar fill colour | steelblue | Any R colour name or hex code — *bar mode only* |
| Bubble SIZE column | Auto-detected | Numeric column mapped to bubble area — *bubble mode only* |
| Bubble COLOUR intensity column | Auto-detected | Numeric column mapped to fill colour gradient — *bubble mode only* |
| Bubble size range | 3 – 15 | Min and max rendered bubble diameter — *bubble mode only* |
| Colour palette | Blues | Fill colour gradient palette — *bubble mode only* |
| Plot height (px) | 600 | Display height in the browser |

### Chord Diagram Settings

| Input | Default | Description |
|---|---|---|
| Top N pathways | 10 | Pathways included, selected by significance |
| Minimum appearances | 1 | Minimum pathway count for a molecule to appear |
| Inner circle size | 0.4 | Radius of the central hole (0.1–0.7) |
| Pathway label font size | 0.55 | cex value for pathway sector labels |
| Molecule label font size | 0.70 | cex value for molecule sector labels |
| Pathway sector colour | tomato | Base colour for pathway arcs |
| Molecule sector colour | steelblue | Base colour for molecule arcs |
| Chord plot height (px) | 750 | Display height in the browser |

### Heatmap Settings

| Input | Default | Description |
|---|---|---|
| Select category | First entry | Pathway/function to display |
| Colour per group | purple, darkorange, … | One text input per group; accepts any R colour name or hex code |
| Heatmap height (px) | 700 | Display height in the browser |

### Boxplot / Violin Settings

| Input | Default | Description |
|---|---|---|
| **Plot type** | Boxplot | Toggle between Boxplot and Violin plot |
| Select protein / gene | First molecule | Dropdown from all molecules in the enrichment file |
| Plot height (px) | 500 | Display height in the browser |

> Group colours for the boxplot and violin plot are the same inputs defined
> under **Heatmap Settings**.

---

## Download Outputs

Every visualization has two download buttons — one in the sidebar and one
below the plot in the main panel. Both produce identical files.

| Plot | Format | Size | DPI | Filename |
|---|---|---|---|---|
| Bar Plot | PNG | 10 × 8 in | 300 | `Barplot_Enriched_Pathways.png` |
| Bubble Plot | PNG | 10 × 8 in | 300 | `Bubbleplot_Enriched_Pathways.png` |
| Chord Diagram | PNG | 10 × 10 in | 300 | `Chord_Diagram_Pathways_Molecules.png` |
| Heatmap (selected) | PNG | 8 × 10 in | 300 | `Heatmap_<category>.png` |
| Heatmap (all) | ZIP of PNGs | 8 × 10 in each | 300 | `All_Heatmaps.zip` |
| Boxplot | PNG | 8 × 6 in | 300 | `Boxplot_<protein>.png` |
| Violin Plot | PNG | 8 × 6 in | 300 | `Violin_<protein>.png` |

> The filename adapts automatically to the active plot type — switching to
> Bubble plot or Violin plot mode updates the filename before download without
> any extra steps.
---

## Tips and Troubleshooting

**No proteins matched in heatmap**
> Check that the Gene/Protein ID column in the normalized data file uses the
> same identifiers (e.g. gene symbols) as the Molecules column in the
> enrichment file. Both are case-sensitive.

**Chord diagram is too cluttered**
> Increase *Minimum appearances* in Chord Diagram Settings to show only
> molecules shared across multiple pathways. Also try reducing *Top N pathways*.

**p-value column not auto-detected**
> Manually select the correct column from the dropdown. The status line in the
> Bar / Bubble Plot tab confirms which transformation is being applied.

**Sample columns not found**
> The Sample ID column in the annotation file must contain values that exactly
> match column names in the normalized data file.

**Heatmap shows fewer proteins than expected**
> The app requires at least 2 matched proteins to draw a heatmap. Proteins
> present in the pathway list but absent from the normalized data file are
> silently skipped.

**Boxplot or violin plot shows a flat line or single point**
> Only one sample was matched for that protein. Verify that sample IDs in the
> annotation file match the normalized data column names.

**Bubble plot shows all bubbles the same size**
> The selected size column may contain non-numeric or missing values. Check
> that the column contains valid numbers for all displayed pathways.

**Bubble colour gradient is not visible**
> All values in the colour column may be identical or nearly identical after
> transformation. Try selecting a different column or a higher-contrast palette.

**Violin plot shows a very narrow shape**
> This typically happens with very small groups (n < 4). Consider switching to
> Boxplot mode for small sample sizes, where the median and IQR are easier to
> interpret than a kernel density estimate.

**Download produces a blank or corrupt file**
> Ensure the plot renders correctly in the browser before downloading.
> The download uses the same rendering function as the display.

**Deployment — "App not found" or creates a duplicate**
> Make sure `appName = "EnrichViz"` matches the existing app name exactly,
> including capitalisation.

**Deployment — "Unauthorized" error**
> Re-run `rsconnect::setAccountInfo()` with a fresh token from your
> shinyapps.io dashboard.

**Deployment — old version still shows after update**
> Hard-refresh the browser (`Ctrl + Shift + R` on Windows/Linux,
> `Cmd + Shift + R` on Mac) to clear the cache.

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
