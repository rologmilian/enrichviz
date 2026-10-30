# EnrichViz Shiny App — v1.4.2

Interactive visualization of enriched pathways or biological functions from
normalized proteomics, transcriptomics, or metabolomics data.

Copyright 2026 RGM  
Licensed under the MIT License — see [License](#license) section below.

---

## Table of Contents

1. [Overview](#overview)
2. [What's New in v1.4.2](#whats-new-in-v142)
3. [What's New in v1.4.1](#whats-new-in-v141)
4. [What's New in v1.4.0](#whats-new-in-v140)
5. [What's New in v1.3.1](#whats-new-in-v131)
6. [What's New in v1.3.0](#whats-new-in-v130)
7. [What's New in v1.2.0](#whats-new-in-v120)
8. [Requirements](#requirements)
9. [Installation](#installation)
10. [Input Files](#input-files)
11. [App Structure](#app-structure)
12. [Tabs and Features](#tabs-and-features)
   - [Tab 1 -- Bar Plot](#tab-1--bar-plot)
   - [Tab 2 -- Chord Diagram](#tab-2--chord-diagram)
   - [Tab 3 -- Molecule Frequency](#tab-3--molecule-frequency)
   - [Tab 4 -- Heatmap](#tab-4--heatmap)
   - [Tab 5 -- Boxplot](#tab-5--boxplot)
   - [Tab 6 -- Volcano](#tab-6--volcano)
   - [Tab 7 -- PCA](#tab-7--pca)
   - [Tab 8 -- DE Heatmap](#tab-8--de-heatmap)
13. [Sidebar Controls Reference](#sidebar-controls-reference)
14. [Download Outputs](#download-outputs)
15. [Deployment](#deployment)
16. [Tips and Troubleshooting](#tips-and-troubleshooting)
17. [License](#license)

---

## Overview

EnrichViz (v1.4.2) is an R Shiny application for interactive exploration and
visualization of pathway or functional enrichment results derived from
normalized proteomics, transcriptomics, or metabolomics experiments.

**Citation:** Garcia-Milian R. [EnrichViz: An Interactive R Shiny Application
for Visualization of Pathway Enrichment Results from Omics
Data](https://pubmed.ncbi.nlm.nih.gov/42395511/). bioRxiv [Preprint]. 2026
Jun 24:2026.06.19.733398. doi: 10.64898/2026.06.19.733398. PMID: 42395511;
PMCID: PMC13320841.

The app accepts three required CSV input files (a normalized data matrix, a
sample metadata/annotation table, and an enrichment results table) plus an
optional fourth file with differential analysis results. It produces eight
coordinated views:

- A **Bar Plot** (optionally rendered as a bubble chart) ranking pathways or
  functions by a user-selected metric.
- A **Chord Diagram** showing co-membership of molecules across pathways or
  functions.
- A **Molecule Frequency Plot** displaying how many pathways each gene,
  protein or metabolite appears in, with controls to restrict the analysis to
  all pathways or only the top N most significant ones, and bars optionally
  coloured by log2 fold change from the differential analysis file.
- A **Heatmap** of molecule-level abundances grouped by pathway or function.
- A **Boxplot** (optionally rendered as a violin plot) displaying the
  distribution of abundance values for a selected molecule across groups.
- A **Volcano Plot** of the differential analysis results, coloured by
  fold-change and adjusted p-value cutoffs *(requires the optional DE file)*.
- A **PCA** plot of samples, coloured by group.
- A **DE Heatmap** of the features passing the fold-change and adjusted
  p-value cutoffs *(requires the optional DE file)*.

Throughout the app, "feature" means whatever the rows of your data are:
genes, transcripts, proteins, or metabolites.

All plots can be downloaded as publication-ready PNG files directly
from the sidebar or from within each tab.

The page header reads *"Heatmap, Bar/Bubble, Chord, Molecule Frequency,
Box/Violin, Volcano, and PCA Plots of Pathway and Differential Analysis
Results"*, and the browser tab title is *"EnrichViz v1.4.2"*.

---

## What's New in v1.4.2

### Molecule Frequency (renamed from Gene Frequency), coloured by fold change

- The **Gene Frequency** tab is renamed **Molecule Frequency**, since rows can
  be genes, proteins or metabolites. The plot title, x-axis ("Molecule"),
  sidebar labels, status line and download file name
  (`Molecule_Frequency_YYYYMMDD_HHMMSS.png`) follow. Internal input and output
  IDs (`freq_*`, `gene_freq_plot`) are unchanged.
- **Colour bars by log2 fold change** (new checkbox, on by default): when the
  optional differential analysis file is loaded, each bar is coloured by the
  molecule's log2 fold change on a diverging scale centred at 0, using the
  Volcano plot's Down / Up colours. Molecules not found in the DE file are
  grey, and the status line reports how many of the shown molecules were
  matched. Without a DE file, bars keep a single colour.

### Chord diagram labels and colours

- Long pathway names are cut to **Max characters in pathway labels** (new
  input, default 35) and end in `…`, and the plot leaves more space around the
  circle, so labels are no longer clipped at the edges.
- Sector colours now range from a darker to a lighter shade of the chosen
  colour. Previously they faded to light yellow, so the last pathways and their
  links were nearly white.

### Box / violin plots

- The x-axis is labelled with the group column name.
- Violins end at the lowest and highest observed values (`trim = TRUE`), so
  count data no longer show density tails below zero.
- The status line says "Selected molecule" (was "Selected protein").

---

## What's New in v1.4.1

### Chord diagram with many molecules

- **Fix:** the chord diagram no longer fails with *"Maybe your `gap.degree` is
  too large…"* when many pathways and molecules are shown. The gap after each
  sector was a fixed 2°, so more than ~180 sectors needed more than the full
  circle (found with a GSEA result: 9 pathways and 187 molecules at the default
  settings). Gaps now shrink as sectors are added and use at most a quarter of
  the circle.
- The Chord Diagram status line adds a tip when more than 100 molecules are
  shown: raise **Minimum appearances** for readable labels.

---

## What's New in v1.4.0

### Differential analysis: Volcano, PCA and DE Heatmap

A fourth, **optional** upload accepts a differential analysis results table
(DESeq2, edgeR, limma, MSstats, MetaboAnalyst, ...) and enables three new tabs.
The app works with genes, proteins or metabolites ("features").

- **Differential Analysis panel** (green, in the sidebar): auto-detected
  feature ID, log2 fold-change and adjusted p-value columns; cutoffs
  (default |log2FC| >= 1, padj < 0.05); live Up / Down count; and a
  **log2(x + 1)** toggle for raw-scale data (untick for already-logged
  proteomics or metabolomics data).
- **Volcano** tab: Up / Down / not-significant colouring, cutoff lines,
  top-N and user-typed labels (non-overlapping via `ggrepel`).
- **PCA** tab: all features, top N most variable, or DE-significant
  features; optional scaling; choice of PCs; group colours; optional sample
  labels. Works without the DE file.
- **DE Heatmap** tab: Z-scored heatmap of significant features (up, down or
  both; top N by padj), sharing the pathway heatmap's clustering, scale and
  colour settings.

See [Tab 6](#tab-6--volcano), [Tab 7](#tab-7--pca) and
[Tab 8](#tab-8--de-heatmap) for details.

### Batch heatmap download limited to filtered pathways

**Download ALL Heatmaps** is replaced by **Download Filtered Heatmaps (N)
(.zip)**. It is only available while a keyword filter is active, exports
only the matching pathways, and is capped at 100 to avoid very long
downloads and server time-outs. Pathways that cannot be drawn are skipped
instead of aborting the archive. See [Batch heatmap
download](#batch-heatmap-download).

### Faster large dropdowns

The heatmap pathway selector and boxplot molecule selector use server-side
search, so files with thousands of pathways or features no longer slow the
browser. Both keep the current selection when the list changes.

### Changes and fixes

- **Keyword filter matching:** pathway-name keywords now match as plain text
  (characters like `(` or `+` no longer cause errors); molecule keywords
  must match a whole member symbol (`MAPK1` no longer matches `MAPK10`).
- **Status lines** for the Bar Plot and Gene Frequency tabs now report
  keyword-filtered pathway counts. The Gene Frequency status notes that the
  metadata filter does not apply to it.
- **Enrichment file detection** now recognises exports from IPA,
  clusterProfiler, g:Profiler, Enrichr, DAVID and fgsea: the pathway-name,
  member and significance columns (adjusted p-values preferred) and the
  member separator are chosen automatically, and lines above the header
  (e.g. IPA's copyright line) are skipped. Files without a member column
  (e.g. MetaboAnalyst) load for the Bar Plot, and the other views explain
  what is missing. A wrong column choice shows a message instead of
  disconnecting the app. See [Supported enrichment
  tools](#supported-enrichment-tools).
- **Heatmaps** hide sample names when there are more than 60 samples.
- **Long pathway names** on the bar and bubble plots are wrapped (or
  truncated) so a single very long name no longer squeezes the plot.
- **Bar Plot colouring:** bars can be coloured by any numeric column
  (z-score, NES, set size, FDR, ...), with a diverging scale centred at 0 for
  signed values and a gradient otherwise.
- **Molecule keyword filter** accepts `*` wildcards (`SYN*` matches `SYN1`,
  `SYN3`, ...).
- **Gene Frequency** can use either the keyword-filtered pathways or all
  pathways (new **Apply the pathway keyword filter** option); the status
  line reports which.
- **Chord Diagram:** "Inner circle size" is renamed **Link curve depth**, with
  help text explaining that it only changes how links curve (it is not a
  filter). Range is now 0.2-0.9.
- **Heatmap colour scale limit** has help text explaining that lower values
  darken the heatmap and higher values lighten it.
- The sidebar help under **Select Pathway / Functions -- Column** now
  describes the automatic column detection.
- The bar-plot download no longer uses `type = "cairo-png"`, so it works on
  any server.
- **Requirements:** the app loads `tidyverse` (not the individual packages
  previously listed) and now also needs `ggrepel`.

---

## What's New in v1.3.1

### Dual Keyword Filter for Pathways

This release adds a **Keyword Filter — Pathways** panel, styled in yellow and
placed just below the **Bar / Bubble Plot Settings** header in the sidebar.
It allows users to restrict the pathway or function table used by all
downstream plots without modifying the uploaded enrichment file.

Two independent text inputs are provided, intentionally kept separate to avoid
ambiguous matches:

| Text box | Searches in | Example input |
|---|---|---|
| **Filter by pathway name** | Pathway / function name column (`col_category`) | `immune, signaling, metabolism` |
| **Filter by molecule / gene symbol** | Molecule members column (`col_molecules`) | `MAPK1, TP53, EGFR` |

#### Why two separate boxes?

Combining both searches into a single field would cause a gene symbol such as
`MAPK1` to match pathways where that symbol appears anywhere — including
pathways whose **name** does not relate to MAPK signalling but whose member
list happens to contain `MAPK1`. Keeping the two filters separate gives full
control: search by pathway concept, by gene membership, or by both
simultaneously.

#### Logic

- Keywords within each box are separated by commas and matched with **OR
  logic** (case-insensitive).
- **Pathway name** keywords are matched as plain text anywhere in the name
  (e.g. `signal` matches "Signaling by WNT"). Special characters such as
  `(`, `+` or `.` are treated literally, not as regular expressions.
- **Molecule** keywords must match a **whole member symbol** — `MAPK1` matches
  pathways containing `MAPK1` but not ones that only contain `MAPK10` or
  `MAPK11`. Members are split using the selected molecule separator.
- Add `*` as a **wildcard** to match a symbol root: `SYN*` matches `SYN1`,
  `SYN2`, `SYN3` -- and also any other symbol starting with "SYN", such as
  `SYNPO` or `SYNGAP1`. `*` can appear anywhere (e.g. `HLA-D*`, `*KIN`).
- If **both** boxes contain text, they combine with **AND logic**: a pathway
  must satisfy both conditions to be retained.
- Leaving a box blank imposes no constraint from that box.
- Leaving **both** boxes blank restores the full pathway table — behaviour is
  identical to previous versions.

#### Examples

| Pathway name filter | Molecule filter | Result |
|---|---|---|
| `MAPK` | *(blank)* | Only pathways whose **name** contains "MAPK" |
| *(blank)* | `MAPK1, MAPK3` | All pathways that **contain** MAPK1 or MAPK3 as members |
| *(blank)* | `SYN*` | All pathways containing any member whose symbol starts with "SYN" |
| `signaling` | `TP53` | Pathways **named** with "signaling" **AND** containing TP53 as a member |
| *(both blank)* | *(both blank)* | All pathways — no filtering |

#### Live summary badge

A summary badge immediately below the two text inputs shows how many pathways
are currently retained out of the total — for example:

```
✅  12 / 87 pathways retained — pathway name: signaling  AND  molecule: TP53
ℹ️  Using all 87 pathways (no keyword filter active).
```

#### What is affected

The filter is applied through a single `ipa_funct_filtered()` reactive that
all downstream consumers use:

| Plot / control | Effect |
|---|---|
| Bar / Bubble Plot | Only filtered pathways are ranked and displayed |
| Chord Diagram | Only filtered pathways contribute edges and molecules |
| Gene Frequency Plot | Frequency counts are computed from filtered pathways only |
| Heatmap category selector | Dropdown lists only the filtered pathway names |
| Batch heatmap download | Enabled only while a filter is active; exports only the filtered pathways (max 100) |

#### What is not affected

- Column-selector widgets (`col_pvalue`, `col_molecules`, `col_category`, etc.)
  always read from the unfiltered `ipa_funct()` so that column detection is
  never disrupted by an active filter.

The keyword filter also controls the batch heatmap download: **Download
Filtered Heatmaps (.zip)** is only available while a keyword filter is
active, and exports only the matching pathways (see [Download
Outputs](#download-outputs)).

#### Implementation notes

- The reactive `ipa_funct_filtered()` is defined immediately after the three
  file-loader reactives and before all plot reactives.
- Each keyword box splits its input on commas and lower-cases the keywords.
- Pathway-name keywords are tested with `grepl(..., fixed = TRUE)` on the
  lower-cased name (literal substring match, no regex).
- Molecule keywords are compared with `%in%` against the lower-cased,
  trimmed member symbols after splitting on `mol_separator` (exact match).
  Keywords containing `*` are converted to an anchored regular expression
  (`*` -> `.*`, all other characters matched literally).
- The two filter steps are applied sequentially (pathway name first, then
  molecule), so the AND combination requires no additional logic.
- The summary badge reactive `ui_keyword_filter_summary` mirrors the pattern
  used by the existing metadata filter badge (`ui_filter_summary`).

---

## What's New in v1.3.0

### Gene / Molecule Frequency by Pathway Plot

The main addition in this release is a dedicated **Gene Frequency** tab that
displays a bar chart of how many pathways each gene or molecule belongs to.
This allows users to quickly identify hub molecules that are shared across many
biological processes.

#### How it works

1. Upload the enrichment results CSV file.
2. Navigate to the **Gene Frequency** tab, or configure the plot using the
   **Gene Frequency Plot Settings** section in the sidebar.
3. Choose the **pathway scope** — either use all available pathways or restrict
   the analysis to the **top N most significant pathways** (ranked by the
   selected p-value column).
4. Set a **minimum pathway count** to filter out genes that appear in only a
   few pathways and reduce visual clutter.
5. Set the **top N genes** to cap the number of bars displayed, showing only
   the most highly connected molecules.
6. The plot updates reactively whenever any of the above controls change.

#### Controls

| Control | Description | Default |
|---|---|---|
| Pathway scope | Use all pathways or restrict to Top N | All pathways |
| Top N pathways | Number of top significant pathways to include (visible only when scope = Top N) | 20 |
| Minimum # pathways | Only show genes appearing in at least this many pathways | 5 |
| Show top N genes | Maximum number of genes displayed after filtering | 50 |
| Plot height (px) | Adjustable rendering height | 600 |

#### Status line

The status line above the plot reports:

- Which pathway scope is active (all or top N, with the p-value column named)
- Total number of genes/molecules found across selected pathways
- Number of genes meeting the minimum pathway threshold
- Number of genes actually displayed

#### Download

The frequency plot can be downloaded as a PNG (12 x 7 in, 300 DPI) using the
**Download Frequency Plot (.png)** button in the sidebar or the
**Download Plot (.png)** button at the bottom of the Gene Frequency tab.

File name format:

```
Gene_Frequency_YYYYMMDD_HHMMSS.png
```

#### Implementation notes

- The reactive `gene_freq_data()` computes, for each molecule, the number of
  **distinct pathways** it is assigned to (duplicate pathway-molecule pairs are
  removed with `distinct()` before counting).
- When scope is set to **Top N**, the same p-value column and transformation
  logic used by the Bar Plot and Chord Diagram is applied (raw p-values are
  auto-transformed with -log10).
- The molecule separator selected in the **Pathway / Functions** section of the
  sidebar is respected.
- The plot is rendered with `ggplot2` (`geom_col`) and uses `theme_classic()`
  with 90 degree rotated x-axis labels for readability with large gene sets.

---

## What's New in v1.2.0

### Metadata Slice Filter

The most significant addition in v1.2.0 is a **Filter Samples by
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
   samples are retained out of the total (e.g. `24 / 48 samples retained
   -- sex: male`).
7. Choose `(none)` in the filter column drop-down to deactivate all filtering
   and restore the full sample set.

#### What is affected

The filter propagates automatically to **all plots**:

| Plot | Effect |
|---|---|
| Bar Plot | Status line reports the active filter |
| Chord Diagram | Status line reports the active filter |
| Gene Frequency | Not affected — frequencies depend only on the pathway table. The status line notes the active filter and that it is not applied |
| Heatmap | Only filtered samples appear as columns; filter is shown in the plot title |
| Boxplot / Violin | Only filtered samples contribute data points; filter is shown in the plot subtitle |

#### Typical use case

A metadata file contains columns `sex` (`male` / `female`) and `diagnosis`
(`disease` / `healthy`). To compare disease vs. healthy **within males only**:

1. Set **Filter column** -> `sex`.
2. Uncheck `female`, keep `male`.
3. Set **Group column** -> `diagnosis`.
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
| tidyverse | >= 2.0.0 (provides ggplot2, dplyr, tidyr, readr, stringr) |
| circlize | >= 0.4.15 |
| pheatmap | >= 1.0.12 |
| ggrepel | >= 0.9.0 (non-overlapping labels on volcano and PCA plots) |

The app loads these five packages with `library()`. Colour palettes
(viridis, ColorBrewer) are provided through ggplot2's own dependencies, so
no extra palette packages need to be installed.

Install all R dependencies at once:

```r
install.packages(c("shiny", "tidyverse", "circlize", "pheatmap", "ggrepel"))
```

> **Note:** v1.2.0 replaced `ComplexHeatmap` (Bioconductor) with `pheatmap`
> (CRAN), removing the Bioconductor dependency and simplifying installation.

---

## Installation

1. Clone or download this repository.
2. Place `app.R` (v1.4.2) in a dedicated project folder.
3. Install all dependencies listed above.
4. Launch the app from R or RStudio:

```r
shiny::runApp("path/to/project/folder")
```

Or open `app.R` in RStudio and click **Run App**.

---

## Input Files

EnrichViz expects **three CSV files** (`.csv`), plus an **optional fourth**
file for differential analysis. All files must have a header row. Maximum
upload size is 500 MB per file.

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
  supported and can be used as **metadata slice filters** (introduced in v1.2.0).

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

- One column contains **pathway or function names**.
- One column contains a **p-value, adjusted p-value, or -log10(p-value)**;
  the app detects whether values are raw (0-1) or already -log10.
- One column contains **member molecules / genes** as a delimited string
  (comma, slash, semicolon, or pipe). This column is optional: without it
  only the Bar Plot works.

The app detects all three columns and the separator automatically; every
choice can be changed in the sidebar. Any lines above the header row (such
as IPA's copyright line) are skipped.

Example layout:

```
Pathway            pvalue   Molecules
Glycolysis         0.001    ProtA,ProtB,ProtC
Oxidative_Phospho  0.003    ProtB,ProtD
```

> **Note:** Molecule identifiers must match exactly (case-sensitive) the
> identifiers in the normalized data file.

#### Supported enrichment tools

Exports from these tools are recognised without changing any setting
(the file must be saved as `.csv`):

| Tool | Pathway name | Members (separator) | Significance column used |
|---|---|---|---|
| IPA (Canonical Pathways export) | `Ingenuity Canonical Pathways` | `Molecules` (`,`) | `-log(p-value)` |
| clusterProfiler, over-representation | `Description` | `geneID` (`/`) | `p.adjust` |
| clusterProfiler, GSEA | `Description` | `core_enrichment` (`/`) | `p.adjust` |
| g:Profiler | `term_name` | `intersections` (`,`) | `adjusted_p_value` |
| Enrichr | `Term` | `Genes` (`;`) | `Adjusted P-value` |
| DAVID (Functional Annotation Chart) | `Term` | `Genes` (`,`) | `FDR` |
| fgsea | `pathway` | `leadingEdge` (`,`) | `padj` |
| MetaboAnalyst (enrichment) | first column | *none in export* -- Bar Plot only | `FDR` |

How detection works:

- **Pathway name:** a text column named like `Description`, `term_name`,
  `Pathway`, `Term`, `gene_set` or `Name`; otherwise the text column with
  the most distinct values. Columns that hold member lists are never
  chosen.
- **Members:** a text column named like `core_enrichment`, `geneID`,
  `Molecules`, `members`, `intersections`, `leadingEdge` or `Genes` that
  actually contains delimited lists; otherwise the text column with the most
  separators. If none exists, the selector shows **(none)** and the views that
  need members explain why they are empty.
- **Significance:** adjusted p-values (`padj`, `p.adjust`, `Adjusted
  P-value`, `FDR`, `Benjamini`, `qvalue`, ...) are preferred over raw
  p-values, then -log10 columns. Only numeric columns are considered.
- **Separator:** the character (`,` `/` `;` `|`) found most often in the
  members column.

Other tools generally work too if their export has a pathway-name column, a
p-value column and a delimited member column. Tab-delimited exports
(e.g. Enrichr `.txt`, IPA `.xls`) must first be saved as CSV.

### 4. Differential Analysis Results File (optional)

A results table from a differential expression / abundance analysis (e.g.
DESeq2, edgeR, limma, MSstats, MetaboAnalyst), one row per feature, where:

- One column contains **feature IDs** (gene symbols, protein accessions,
  metabolite names, ...). These must match the IDs in the selected
  Gene/Protein ID column of the normalized data file.
- One column contains the **log2 fold change**.
- One column contains the **adjusted p-value** (padj / FDR / q-value).
- Other columns are allowed and ignored.

Example layout (DESeq2):

```
symbol   baseMean  log2FoldChange  lfcSE  stat    pvalue    padj
NPAS4    63.78     -3.97           0.34   -11.74  7.6e-32   1.7e-27
FOSB     229.78    -2.46           0.24   -10.41  2.4e-25   2.6e-21
```

> The three columns are auto-detected from common names (`symbol`,
> `protein`, `metabolite`, `log2FoldChange`, `logFC`, `padj`, `FDR`,
> `adj.P.Val`, `qvalue`, ...) and can be overridden in the sidebar. If the
> normalized data's ID column name also exists in the DE file, it is selected
> by default. Rows with a missing fold change or adjusted p-value (e.g. DESeq2
> independent filtering) are dropped.

---

## App Structure

The app uses a standard Shiny `fluidPage` layout with a fixed left sidebar
(`sidebarPanel`, width = 3) and a fixed main content area (`mainPanel`,
width = 9). Both panels scroll independently.

A fixed blurb bar below the title summarises the three required input files.

The main panel contains a `tabsetPanel` with **eight tabs**:

| Tab index | Tab label |
|---|---|
| 1 | Bar Plot |
| 2 | Chord Diagram |
| 3 | Molecule Frequency |
| 4 | Heatmap |
| 5 | Boxplot |
| 6 | Volcano |
| 7 | PCA |
| 8 | DE Heatmap |

The sidebar is divided into clearly labelled sections:

- **Upload Data Files** -- three required `fileInput` widgets plus one
  optional widget for differential analysis results.
- **Filter Samples by Metadata** *(introduced in v1.2.0)* -- column picker,
  value checkboxes, Select all / Clear all buttons, live sample count badge.
- **Select Pathway / Functions -- Column** -- category and molecule column
  selectors, separator radio buttons.
- **Normalized counts or protein abundance -- Column** -- gene/protein ID
  column selector.
- **Sample Annotation -- Column** -- sample ID and group column selectors.
- **Bar / Bubble Plot Settings** -- plot type, p-value column, top N,
  colour / palette controls, and optional colouring of bars by a numeric
  column.
- **Keyword Filter — Pathways** *(new in v1.3.1)* -- two independent text
  inputs for filtering by pathway name and by molecule/gene symbol; live
  count badge; affects all downstream plots.
- **Chord Diagram Settings** -- top N, minimum frequency, link curve depth,
  label sizes, maximum pathway-label length, sector colours.
- **Molecule Frequency Plot Settings** *(new in v1.3.0; renamed in v1.4.2)* --
  pathway scope toggle, optional top N pathways, minimum pathway count, top N
  molecules, plot height, colour by log2 fold change, download button.
- **Heatmap Settings** -- category selector, per-group colour inputs, height,
  single-heatmap download, and the batch download of keyword-filtered
  heatmaps,
  download buttons.
- **Boxplot / Violin Settings** -- molecule selector, plot type, height,
  download button.
- **Differential Analysis** (green panel) -- DE column selectors,
  fold-change and adjusted p-value cutoffs, up/down count badge, and the
  log2(x + 1) toggle used by PCA and the DE heatmap.
- **Volcano Plot Settings** -- labels, colours, height, download button.
- **PCA Settings** -- feature set, scaling, PCs to plot, sample labels,
  height, download button.
- **DE Heatmap Settings** -- direction, maximum features, height, download
  button.

### Large dropdown lists

The heatmap pathway selector and the boxplot molecule selector use
**server-side search** (`updateSelectizeInput(..., server = TRUE)`): the full
list stays on the server, and only the entries matching what you type are sent
to the browser. This keeps the app responsive with enrichment files of
thousands of gene sets and data sets of tens of thousands of features. Click a
selector and start typing to find an entry.

---

## Tabs and Features

### Tab 1 -- Bar Plot

Displays pathways or functions ranked by a chosen p-value metric.

| Feature | Details |
|---|---|
| Plot type toggle | Switch between a horizontal bar chart and a bubble chart (`radioButtons`) |
| P-value column | Select the significance column from the enrichment file (`selectInput`); raw p-values (0-1) are auto-transformed with -log10 |
| Auto-transform indicator | `helpText` confirms whether -log10 will be applied |
| Top N selector | Cap the number of displayed pathways (`numericInput`, default 20) |
| Long pathway names | **Wrap** (default) splits names longer than *Max characters* (default 50) over up to 3 lines, breaking at spaces, `_` or `-` so GSEA-style names (`GOBP_..._...`) wrap too; **Truncate** cuts them to one line ending in `…`; **Show in full** leaves them unchanged. Applies to bar and bubble plots and their downloads. Only the axis labels change -- pathways whose shortened names coincide are still drawn as separate bars |
| Colour bars by | Optional numeric column used to colour the bars, like the bubble plot (`selectInput`, default **None**). Useful for z-score (IPA), NES (GSEA), set size / count, or an adjusted p-value when the x-axis shows something else |
| Colour scale | **Auto** (default): diverging blue-white-orange centred at 0 (white = 0) when the column contains negative values, e.g. z-score or NES; otherwise a light-to-dark gradient. **Gradient** or **Diverging** can be forced (`radioButtons`) |
| Gradient palette / Negative & positive colours | Palette for gradients (`selectInput`, default `Blues`); end colours for diverging scales (`textInput`, defaults `steelblue` / `darkorange`) |
| P-value colour columns | Columns named like p-values / FDR / q-values with values in 0-1 are shown as -log10, so darker = more significant |
| Missing values | Pathways without a value in the colour column are drawn in grey |
| Bar colour | Single fill colour used when **Colour bars by** = None (`textInput`, default `steelblue`) |
| Bubble size column | Column mapped to bubble area (`selectInput`) |
| Bubble colour column | Column mapped to bubble fill intensity (`selectInput`) |
| Bubble size range | Min/max point size (`sliderInput`, default 3-15) |
| Bubble palette | Sequential or diverging colour palette (`selectInput`) |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 600) |
| Keyword filter | Pathway name and/or molecule filter applied upstream via `ipa_funct_filtered()` *(new in v1.3.1)* |
| Status line | Reports active metadata filter, plot type, p-value column, transformation applied, total and displayed pathway counts |

### Tab 2 -- Chord Diagram

Renders a chord diagram (via `circlize`) showing molecule co-membership across
the top N pathways or functions.

| Feature | Details |
|---|---|
| Top N pathways | Number of highest-ranked pathways to include (`numericInput`, default 10) |
| Minimum frequency | Show only molecules appearing in at least N pathways (`numericInput`, default 2) |
| Link curve depth | How far the links bend toward the centre (`sliderInput`, 0.2-0.9, default 0.4). Lower values keep links near the outer ring, which can look crowded; higher values arc links through the centre so they are easier to follow. Appearance only -- it does not filter pathways or molecules. (Previously labelled "Inner circle size") |
| Pathway label size | Font size for pathway sector labels (`numericInput`, default 0.55) |
| Molecule label size | Font size for molecule labels (`numericInput`, default 0.70) |
| Max characters in pathway labels | Longer pathway names are cut and end in `…` so they fit around the circle (`numericInput`, default 35) *(new in v1.4.2)* |
| Pathway sector colour | Hex or named colour for pathway arcs (`textInput`, default `tomato`); sectors are shaded from a darker to a lighter tone of it |
| Molecule sector colour | Hex or named colour for molecule arcs (`textInput`, default `steelblue`) |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 750) |
| Keyword filter | Pathway name and/or molecule filter applied upstream via `ipa_funct_filtered()` *(new in v1.3.1)* |
| Status line | Reports active metadata filter, pathway count, molecule count, total connections, minimum frequency, and link curve depth; suggests raising the minimum frequency when more than 100 molecules are shown |

### Tab 3 -- Molecule Frequency

Displays a bar chart of how many pathways each molecule (gene, protein or
metabolite) belongs to, allowing rapid identification of hub molecules shared
by many biological processes. Renamed from "Gene Frequency" in v1.4.2.

| Feature | Details |
|---|---|
| Pathway scope | Toggle between using **all pathways** or only the **top N most significant** (`radioButtons`) |
| Apply keyword filter | When ticked (default), frequencies use only the pathways retained by the Keyword Filter; untick to use all pathways even while a keyword filter is active (`checkboxInput`) |
| Top N pathways | Number of top-ranked pathways to use when scope = Top N (`numericInput`, default 20; only visible when scope = Top N) |
| Minimum pathway count | Only display molecules appearing in at least this many pathways (`numericInput`, default 5) |
| Top N molecules | Maximum number of molecules shown after filtering (`numericInput`, default 50) |
| Colour by log2 fold change | When the DE file is loaded, bars are coloured by log2 fold change on a diverging scale centred at 0 (Volcano Down / Up colours); molecules not in the DE file are grey (`checkboxInput`, default on) *(new in v1.4.2)* |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 600) |
| Keyword filter | Pathway name and/or molecule filter applied upstream via `ipa_funct_filtered()` *(new in v1.3.1)* |
| Status line | Reports active pathway scope, total molecule count, molecules meeting the minimum threshold, number displayed, and how many shown molecules were matched in the DE file |
| Download | PNG (12 x 7 in, 300 DPI); file name includes a timestamp |

### Tab 4 -- Heatmap

Renders a clustered heatmap (via `pheatmap`) of Z-score normalised molecule
abundances for a selected pathway or function, using only the filtered samples.
User can select clustering or not for the columns/samples.

| Feature | Details |
|---|---|
| Category selector | Choose one pathway or function to visualise (`selectizeInput` with type-ahead search); list is restricted to keyword-filtered pathways *(new in v1.3.1)*. The label shows how many pathways are available, and the current choice is kept when the filter changes if it is still in the list |
| Per-group colours | One `textInput` per group (hex or name) for the column annotation bar |
| Row clustering | Hierarchical clustering of molecules (`checkboxInput`, default on) |
| Column clustering | Optional (`checkboxInput`, default off); when off, columns are ordered by group for comparability |
| Colour scale limit | Symmetric ± Z-score limit for the colour scale (`sliderInput`, 1-4, default 2). Z-scores beyond the limit get the darkest colour: **lower** values give a darker, higher-contrast heatmap where small differences stand out; **higher** values give a lighter heatmap with more shades between extremes. Shared with the DE Heatmap |
| Row name display | Shown when <= 100 molecules are matched |
| Column name display | Sample names shown when <= 60 samples (hidden above that to stay readable) |
| Z-score scaling | Applied row-wise before rendering |
| Plot title | Includes the selected category name and active metadata filter |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 700) |
| Batch download | **Download Filtered Heatmaps (N) (.zip)** -- one PNG per keyword-filtered pathway; see below |
| Status line | Reports active metadata filter, samples used, category name, gene list size, matched proteins, and sample column names |

### Tab 5 -- Boxplot

Displays the distribution of normalized abundance values for a single selected
molecule across sample groups, using only the filtered samples.

| Feature | Details |
|---|---|
| Molecule selector | `selectizeInput` with type-ahead search; lists all molecules present in both the enrichment file and the normalized data, with the matched count in the label |
| Plot type toggle | Switch between boxplot and violin plot (`radioButtons`); violins end at the observed minimum and maximum (`trim = TRUE`) |
| Axes | x-axis labelled with the group column name; y-axis "Normalized abundance / counts" |
| Data points overlay | Individual sample points overlaid via `geom_jitter` (always on) |
| Colour by group | Per-group colours shared with the Heatmap Settings section |
| Plot subtitle | Shows the active metadata filter when one is set |
| Plot height | Adjustable rendering height in pixels (`numericInput`, default 500) |
| Status line | Reports active metadata filter, plot type, selected molecule, sample count, groups, and value range |

### Tab 6 -- Volcano

Plots every tested feature from the differential analysis file: log2 fold
change on the x-axis and -log10(adjusted p-value) on the y-axis. Requires the
optional DE file.

| Feature | Details |
|---|---|
| Cutoffs | Features with \|log2FC\| >= cutoff **and** padj < cutoff are coloured Up / Down; all others are grey. Dashed lines mark the cutoffs |
| Legend | Shows the number of Up and Down features |
| Labels | Top N significant features by padj (`numericInput`, default 10) plus any IDs typed in the label box (case-insensitive); non-overlapping via `ggrepel` |
| Colours | Up and Down colours (`textInput`, defaults `firebrick` / `steelblue`) |
| padj = 0 | Capped at the smallest representable number so the point is plotted instead of being lost as infinite |
| Plot height | Adjustable (`numericInput`, default 600) |
| Status line | Reports the DE columns used, cutoffs, number of features tested (and rows removed for missing padj), and Up / Down counts |

The metadata filter does not apply here: the volcano plot shows the DE results
as uploaded.

### Tab 7 -- PCA

Principal component analysis of the samples in the normalized data file,
using only samples retained by the metadata filter. Does **not** require the
DE file unless the DE-significant feature set is chosen.

| Feature | Details |
|---|---|
| Features used | All features, the top N most variable features (default 500), or only DE-significant features passing the cutoffs (`radioButtons`) |
| log2(x + 1) | Applied when the shared toggle in the Differential Analysis panel is on (default). Turn it off for data that are already log-transformed, as is common for proteomics and metabolomics. The app stops with a message if the data contain negative values while the toggle is on |
| Scaling | Features are always centred; optionally scaled to unit variance (`checkboxInput`, default off) |
| PCs | Choose which components go on the x and y axes (default PC1 vs PC2); axis labels show % variance explained |
| Colours | Points coloured by the group column, using the per-group colours from Heatmap Settings |
| Sample labels | Optional, non-overlapping via `ggrepel` (`checkboxInput`, default off) |
| Plot height | Adjustable (`numericInput`, default 600) |
| Status line | Reports active metadata filter, samples, features used, transformation, and variance explained by PC1-PC5 |

Features with missing values or zero variance are removed before PCA. At
least 3 samples are required.

### Tab 8 -- DE Heatmap

A heatmap of the features passing the fold-change and adjusted p-value
cutoffs, drawn with the same engine as the pathway Heatmap tab. Requires the
optional DE file.

| Feature | Details |
|---|---|
| Direction | Up & down, Up only, or Down only (`radioButtons`) |
| Max features | The top N significant features by padj that are present in the normalized data (`numericInput`, default 50) |
| Values | log2(x + 1) (if the shared toggle is on), then row Z-score |
| Shared settings | Row/column clustering, colour-scale limit and group colours come from Heatmap Settings |
| Samples | Only samples retained by the metadata filter |
| Plot height | Adjustable (`numericInput`, default 800) |
| Status line | Reports active metadata filter, cutoffs, direction, number of significant features, number shown, samples used, and transformation |

---

## Sidebar Controls Reference

### File Input Controls

| Input ID | Widget type | Label | Accepts |
|---|---|---|---|
| `file_norm_prot` | `fileInput` | Normalized Data (.csv) | `.csv` |
| `file_annotation` | `fileInput` | Metadata or Sample Annotation (.csv) | `.csv` |
| `file_ipa_funct` | `fileInput` | Enriched Pathway or Functions (.csv) | `.csv` |
| `file_de` | `fileInput` | Differential Analysis Results (.csv) -- optional | `.csv` |

### Metadata Filter Controls *(introduced in v1.2.0)*

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `filter_col` | `selectInput` | Filter column | `(none)` |
| `filter_vals` | `checkboxGroupInput` | Keep values in `<col>` | All values selected |
| `filter_select_all` | `actionButton` | Select all | -- |
| `filter_select_none` | `actionButton` | Clear all | -- |

### Pathway / Function Column Controls

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `col_category` | `selectInput` | Category column | Detected pathway-name column (see [Supported enrichment tools](#supported-enrichment-tools)) |
| `col_molecules` | `selectInput` | Molecules column | Text column matching a keyword (`core_enrichment`, `geneID`, `molecules`, `members`, `genes`, ...); otherwise the text column with the most separators. Numeric columns are never auto-selected; **(none)** if the file has no member list |
| `mol_separator` | `radioButtons` | Gene/molecule separator | The separator found most often in the molecules column (re-detected whenever that column changes); `,` if none found |

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
| `col_pvalue` | `selectInput` | P-value column | Detected significance column; adjusted p-values preferred |
| `bar_top_n` | `numericInput` | Show top N pathways | `20` |
| `label_mode` | `selectInput` | Long pathway names (Wrap / Truncate / Show in full) | `wrap` |
| `label_width` | `numericInput` | Max characters (per line when wrapping) | `50` |
| `bar_fill_col` | `selectInput` | Colour bars by | None (single colour); lists numeric columns |
| `bar_fill_scale` | `radioButtons` | Colour scale | `auto` |
| `bar_palette` | `selectInput` | Gradient palette | `Blues` |
| `bar_neg_color` / `bar_pos_color` | `textInput` | Negative / Positive colour (diverging) | `steelblue` / `darkorange` |
| `bar_color` | `textInput` | Bar fill colour (when not colouring by a column) | `steelblue` |
| `bubble_size_col` | `selectInput` | Bubble SIZE column | Numeric column matching a keyword (`count`, `size`, `ratio`, ...) |
| `bubble_color_col` | `selectInput` | Bubble COLOUR column | Same as the detected significance column |
| `bubble_size_range` | `sliderInput` | Bubble size range | `3` - `15` |
| `bubble_palette` | `selectInput` | Colour palette | `Blues` |
| `bar_plot_height` | `numericInput` | Plot height (px) | `600` |

### Keyword Filter — Pathways *(new in v1.3.1)*

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `pathway_keyword_filter` | `textInput` | Filter by pathway name | *(blank — no filter)* |
| `molecule_keyword_filter` | `textInput` | Filter by molecule / gene symbol | *(blank — no filter)* |

### Chord Diagram Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `chord_top_n` | `numericInput` | Top N pathways | `10` |
| `chord_min_freq` | `numericInput` | Minimum molecule frequency | `2` |
| `chord_inner_radius` | `sliderInput` | Link curve depth | `0.4` (range 0.2-0.9) |
| `chord_label_size_path` | `numericInput` | Pathway label font size | `0.55` |
| `chord_label_size_mol` | `numericInput` | Molecule label font size | `0.70` |
| `chord_label_width` | `numericInput` | Max characters in pathway labels | `35` |
| `chord_pathway_color` | `textInput` | Pathway sector colour | `tomato` |
| `chord_molecule_color` | `textInput` | Molecule sector colour | `steelblue` |
| `chord_plot_height` | `numericInput` | Chord plot height (px) | `750` |

### Molecule Frequency Plot Settings *(new in v1.3.0; renamed in v1.4.2)*

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `freq_pathway_scope` | `radioButtons` | Pathways used to compute molecule frequency | `all` |
| `freq_use_keyword_filter` | `checkboxInput` | Apply the pathway keyword filter (if one is active) | `TRUE` |
| `freq_top_n_pathways` | `numericInput` | Top N pathways (visible when scope = top) | `20` |
| `freq_min_pathways` | `numericInput` | Minimum # pathways a molecule must appear in | `5` |
| `freq_top_n_genes` | `numericInput` | Show top N molecules (after filtering) | `50` |
| `freq_plot_height` | `numericInput` | Frequency plot height (px) | `600` |
| `freq_color_lfc` | `checkboxInput` | Colour bars by log2 fold change (needs the differential analysis file) | `TRUE` |

### Heatmap Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `selected_category` | `selectizeInput` (server-side) | Select: `<category column>` (N) | First category in keyword-filtered enrichment table; current choice kept if still available |
| `heatmap_cluster_rows` | `checkboxInput` | Cluster rows | `TRUE` |
| `heatmap_cluster_cols` | `checkboxInput` | Cluster columns | `FALSE` |
| `heatmap_scale_limit` | `sliderInput` | Colour scale limit (± Z-score); lower = darker, higher = lighter | `2` |
| `color_group_<i>` | `textInput` | Colour -- `<group>` | Cycles through preset palette |
| `plot_height` | `numericInput` | Heatmap height (px) | `700` |

### Boxplot / Violin Settings

| Input ID | Widget type | Label | Default / Notes |
|---|---|---|---|
| `plot_type` | `radioButtons` | Plot type | `box` |
| `boxplot_protein` | `selectizeInput` (server-side) | Search protein / gene / metabolite (N matched) | First matched molecule; current choice kept if still available |
| `boxplot_height` | `numericInput` | Plot height (px) | `500` |

### Differential Analysis

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `de_col_gene` | `selectInput` | Feature ID column | Same name as the normalized data ID column if present, else auto-detected |
| `de_col_lfc` | `selectInput` | log2 fold-change column | Auto-detected by keyword |
| `de_col_padj` | `selectInput` | Adjusted p-value column | Auto-detected by keyword |
| `de_lfc_cut` | `numericInput` | \|log2 fold change\| >= | `1` |
| `de_padj_cut` | `numericInput` | Adjusted p-value < | `0.05` |
| `de_log_transform` | `checkboxInput` | Apply log2(x + 1) to normalized data (PCA & DE heatmap) | `TRUE` |

### Volcano Plot Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `volcano_top_labels` | `numericInput` | Label top N significant features | `10` |
| `volcano_label_genes` | `textInput` | Also label these features | *(blank)* |
| `volcano_up_color` | `textInput` | Up-regulated colour | `firebrick` |
| `volcano_down_color` | `textInput` | Down-regulated colour | `steelblue` |
| `volcano_height` | `numericInput` | Plot height (px) | `600` |

### PCA Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `pca_gene_set` | `radioButtons` | Features used for PCA | `top` |
| `pca_top_n` | `numericInput` | Top N most variable features (visible when set = top) | `500` |
| `pca_scale` | `checkboxInput` | Scale features to unit variance | `FALSE` |
| `pca_x` / `pca_y` | `numericInput` | X-axis PC / Y-axis PC | `1` / `2` |
| `pca_show_labels` | `checkboxInput` | Show sample labels | `FALSE` |
| `pca_height` | `numericInput` | Plot height (px) | `600` |

### DE Heatmap Settings

| Input ID | Widget type | Label | Default |
|---|---|---|---|
| `de_heatmap_direction` | `radioButtons` | Features to show | `both` |
| `de_heatmap_top_n` | `numericInput` | Max features (top N by padj) | `50` |
| `de_heatmap_height` | `numericInput` | Plot height (px) | `800` |

---

## Download Outputs

Each tab exposes two download buttons (sidebar and in-tab footer) for its
respective plot:

| Output ID | Label | Format |
|---|---|---|
| `download_barplot` / `download_barplot_main` | Download Plot (.png) | PNG, 4200 x 2400 px, 300 DPI |
| `download_chord` / `download_chord_main` | Download Chord Diagram (.png) | PNG, 10 x 10 in, 300 DPI |
| `download_gene_freq` / `download_gene_freq_main` | Download Frequency Plot (.png) | PNG, 12 x 7 in, 300 DPI |
| `download_single` / `download_single_main` | Download Selected Heatmap (.png) | PNG, 8 x 10 in, 300 DPI |
| `download_all` / `download_all_main` | Download Filtered Heatmaps (N) (.zip) | ZIP of PNG files (8 x 10 in, 300 DPI), one per keyword-filtered pathway; see below |
| `download_boxplot_sidebar` / `download_boxplot_main` | Download Plot (.png) | PNG, 8 x 6 in, 300 DPI |
| `download_volcano` / `download_volcano_main` | Download Volcano Plot (.png) | PNG, 9 x 7 in, 300 DPI |
| `download_pca` / `download_pca_main` | Download PCA Plot (.png) | PNG, 9 x 7 in, 300 DPI |
| `download_de_heatmap` / `download_de_heatmap_main` | Download DE Heatmap (.png) | PNG, 8 x 10 in, 300 DPI |

File names are auto-generated in the format:

```
Barplot_Enriched_Pathways_YYYYMMDD_HHMMSS.png
Bubbleplot_Enriched_Pathways_YYYYMMDD_HHMMSS.png
Chord_Diagram_Pathways_Molecules.png
Molecule_Frequency_YYYYMMDD_HHMMSS.png
Heatmap_<Category>.png
Filtered_Heatmaps_YYYYMMDD_HHMMSS.zip   (contains Heatmap_<Category>.png files)
Boxplot_<Molecule>.png
Violin_<Molecule>.png
Volcano_YYYYMMDD_HHMMSS.png
PCA_YYYYMMDD_HHMMSS.png
DE_Heatmap_YYYYMMDD_HHMMSS.png
```

### Batch heatmap download

Exporting a heatmap for every pathway in a large enrichment file (often
1,000+ gene sets) takes a long time and can time out on hosted servers, so
the batch download is restricted:

- It is **only available while a keyword filter is active** (pathway name
  and/or molecule). Without a filter, a note explains how to enable it.
- It exports **only the keyword-filtered pathways**, and the button shows how
  many, e.g. *Download Filtered Heatmaps (24) (.zip)*.
- It is capped at **100 pathways**. Above that, the button is replaced by a
  message asking you to narrow the filter.
- Pathways that cannot be drawn (fewer than 2 matched features, or zero
  variance) are skipped.

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

### Shiny Server / Posit Connect

Copy the project folder to the server apps directory and follow the standard
server deployment documentation.

---

## Tips and Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Chord diagram is empty | No molecules exceed the minimum frequency threshold | Lower `chord_min_freq` (default is 2) |
| Molecule Frequency plot is empty | No molecules meet the minimum pathway count | Lower `freq_min_pathways` (default is 5) or switch scope to "All pathways" |
| Molecule Frequency plot shows fewer molecules than expected | Top N pathways scope is too restrictive | Increase `freq_top_n_pathways` or switch to "All pathways" |
| All plots show no data after typing a keyword | Keyword does not match any pathway name or molecule | Check spelling or clear the keyword filter boxes |
| Keyword filter retains fewer pathways than expected | Both filter boxes are active (AND logic) | Clear one box to broaden the match |
| Heatmap category dropdown is empty after filtering | Active keyword filter excluded all pathways | Broaden or clear the keyword filter |
| Heatmap is blank | Selected pathway contains no molecules matching the data matrix | Check identifier matching between input files |
| Heatmap shows fewer columns than expected | Metadata filter is active | Check the filter panel badge; broaden the selection or choose `(none)` |
| Boxplot shows a flat line or single point | Only one sample per group after filtering | Broaden the metadata filter or check group assignments |
| Filter checkbox group is empty | Selected filter column contains only `NA` / blank values | Choose a different filter column |
| App crashes on upload | Wrong file format or encoding | Ensure files are comma-delimited UTF-8 CSV |
| Molecule not found in boxplot selector | Identifier mismatch between files | Check for trailing spaces or case differences |
| Download produces a blank PNG | Plot has not rendered yet | Switch to the desired tab and wait for the plot before downloading |
| `there is no package called 'tidyverse'` on launch | Dependencies installed individually | Run `install.packages("tidyverse")` |
| Volcano / DE Heatmap tabs show nothing | DE file not uploaded | Upload the optional differential analysis file |
| DE Heatmap: "Fewer than 2 significant features found" | Cutoffs too strict, or DE IDs differ from the normalized data IDs (e.g. Entrez vs symbol, UniProt vs gene name) | Relax the cutoffs, or pick matching ID columns in both the DE and Normalized data selectors |
| PCA / DE Heatmap: "Data contain negative values" | Data are already log-transformed | Untick **Apply log2(x + 1)** in the Differential Analysis panel |
| PCA looks dominated by a few features | Very different feature scales (common in metabolomics) | Tick **Scale features to unit variance** |
| No "Download Filtered Heatmaps" button | No keyword filter active, or more than 100 pathways match | Enter a pathway-name or molecule keyword and narrow it until <= 100 pathways match |
| Boxplot / heatmap find no molecules after upload | Wrong Molecules column or separator selected | Check the detected Molecules column and separator in **Select Pathway / Functions -- Column** |
| Molecule keyword filter misses a gene | Molecule filter needs the whole symbol | Type the full symbol (e.g. `MAPK1`), use a wildcard (`MAPK*`), or check the molecule separator setting |
| Bar / bubble plot squeezed to the right | A very long pathway name | Set **Long pathway names** to Wrap or Truncate and lower **Max characters**; increase plot height if wrapped labels look crowded |
| Heatmap looks washed out / too dark | Colour scale limit too high / too low | Lower the **colour scale limit** for more contrast, raise it for a lighter heatmap |
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
