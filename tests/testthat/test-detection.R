# Goal: Check that the columns of each supported enrichment-tool export are recognized.
# Why:  Table 2 of the manuscript and the README list these layouts; a change to the keyword lists
#       must not break them. The tables mimic the column names and separators of each tool.

detect_all <- function(df) {
  mol <- detect_molecules_col(df)
  list(category  = detect_category_col(df),
       molecules = mol,
       pvalue    = detect_pvalue_col(df),
       separator = if (is.null(mol)) NA else names(which.max(separator_density(df[[mol]]))))
}

test_that("clusterProfiler GSEA export is recognized", {
  df <- data.frame(ID = c("GO:1", "GO:2"), Description = c("synapse organization", "cytokine production"),
                   setSize = c(50, 60), NES = c(2.1, -1.8), pvalue = c(1e-5, 1e-4),
                   p.adjust = c(1e-4, 1e-3), leading_edge = c("tags=40%", "tags=30%"),
                   core_enrichment = c("A/B/C", "C/D/E"), check.names = FALSE)
  r <- detect_all(df)
  expect_equal(r$category, "Description"); expect_equal(r$molecules, "core_enrichment")
  expect_equal(r$pvalue, "p.adjust");      expect_equal(r$separator, "/")
})

test_that("clusterProfiler ORA export is recognized", {
  df <- data.frame(ID = "hsa1", Description = "pathway one", GeneRatio = "3/20", BgRatio = "10/500",
                   pvalue = 0.001, p.adjust = 0.01, qvalue = 0.01, geneID = "A/B/C", Count = 3)
  r <- detect_all(df)
  expect_equal(r$molecules, "geneID"); expect_equal(r$pvalue, "p.adjust"); expect_equal(r$separator, "/")
})

test_that("g:Profiler web export is recognized", {
  df <- data.frame(source = "GO:BP", term_name = c("glycolysis", "immune response"), term_id = c("GO:1", "GO:2"),
                   adjusted_p_value = c(0.001, 0.01), negative_log10_of_adjusted_p_value = c(3, 2),
                   term_size = c(50, 200), query_size = 100, intersection_size = c(3, 4),
                   intersections = c("A,B,C", "B,C,D,E"))
  r <- detect_all(df)
  expect_equal(r$category, "term_name"); expect_equal(r$molecules, "intersections")
  expect_equal(r$pvalue, "adjusted_p_value"); expect_equal(r$separator, ",")
})

test_that("Enrichr export is recognized", {
  df <- data.frame(Term = c("t1", "t2"), Overlap = c("3/50", "2/40"), `P-value` = c(1e-4, 1e-3),
                   `Adjusted P-value` = c(1e-3, 1e-2), `Odds Ratio` = c(5, 3),
                   `Combined Score` = c(40, 20), Genes = c("A;B;C", "B;D"), check.names = FALSE)
  r <- detect_all(df)
  expect_equal(r$category, "Term"); expect_equal(r$molecules, "Genes")
  expect_equal(r$pvalue, "Adjusted P-value"); expect_equal(r$separator, ";")
})

test_that("DAVID Functional Annotation Chart export is recognized", {
  df <- data.frame(Category = "GOTERM_BP_DIRECT", Term = c("GO:1~synapse", "GO:2~immune"), Count = c(3, 2),
                   PValue = c(1e-4, 1e-3), Genes = c("A, B, C", "B, D"), Benjamini = c(1e-3, 1e-2),
                   FDR = c(1e-3, 1e-2))
  r <- detect_all(df)
  expect_equal(r$category, "Term"); expect_equal(r$molecules, "Genes"); expect_equal(r$separator, ",")
  expect_true(r$pvalue %in% c("Benjamini", "FDR"))   # both are adjusted p-values
})

test_that("fgsea export is recognized", {
  df <- data.frame(pathway = c("p1", "p2"), pval = c(1e-4, 1e-3), padj = c(1e-3, 1e-2), NES = c(2, -2),
                   size = c(30, 40), leadingEdge = c("A,B,C", "D,E"))
  r <- detect_all(df)
  expect_equal(r$category, "pathway"); expect_equal(r$molecules, "leadingEdge"); expect_equal(r$pvalue, "padj")
})

test_that("IPA canonical pathways: the B-H column is taken as the significance column", {
  df <- data.frame(`Ingenuity Canonical Pathways` = c("Pathway 1", "Pathway 2"),
                   `-log(p-value)` = c(4, 3), `-log(B-H p-value)` = c(2.5, 1.8), Ratio = c(0.1, 0.05),
                   `z-score` = c(2, NA), Molecules = c("A,B,C", "B,D"), check.names = FALSE)
  r <- detect_all(df)
  expect_equal(r$category, "Ingenuity Canonical Pathways"); expect_equal(r$molecules, "Molecules")
  expect_equal(r$pvalue, "-log(B-H p-value)")
  expect_true(is_numeric_col(df$`z-score`))            # blanks do not hide a numeric column
})

test_that("IPA upstream regulators: B-H corrected p-value is preferred over the raw p-value", {
  df <- data.frame(`Upstream Regulator` = c("TNF", "IL6"), `Expr Log Ratio` = c(1.2, NA),
                   `Activation z-score` = c(2.1, 1.5), `p-value of overlap` = c(1e-6, 1e-4),
                   `B-H corrected p-value` = c(1e-4, 1e-2),
                   `Target Molecules in Dataset` = c("A,B,C", "C,D"), check.names = FALSE)
  r <- detect_all(df)
  expect_equal(r$pvalue, "B-H corrected p-value"); expect_equal(r$molecules, "Target Molecules in Dataset")
})

test_that("WebGestalt export: overlapId (Entrez IDs) is detected; userId must be chosen by the user", {
  df <- data.frame(geneSet = "hsa04080", description = "Neuroactive ligand-receptor interaction",
                   size = 370, overlap = 3, enrichmentRatio = 3.2, pValue = 1e-12, FDR = 1e-9,
                   overlapId = "10874;1136;1142", userId = "NMU;CHRNA3;CHRNB3")
  r <- detect_all(df)
  expect_equal(r$category, "description"); expect_equal(r$pvalue, "FDR"); expect_equal(r$molecules, "overlapId")
})

test_that("MetaboAnalyst enrichment export has no molecule list", {
  df <- data.frame(Pathway = c("Glycolysis", "TCA cycle"), total = c(20, 30), hits = c(3, 2),
                   `Raw p` = c(0.001, 0.01), FDR = c(0.01, 0.05), check.names = FALSE)
  r <- detect_all(df)
  expect_null(r$molecules); expect_equal(r$pvalue, "FDR")
})

test_that("a plain CSV is read unchanged and a copyright line above the header is skipped", {
  df <- data.frame(Term = c("t1", "t2"), FDR = c(0.01, 0.02), Genes = c("A;B", "C;D"))
  f  <- tmp_csv(df)
  expect_equal(names(read_enrichment_csv(f)), names(df))
  g  <- tempfile(fileext = ".csv")
  writeLines(c("(c) 2000-2026 QIAGEN. All rights reserved.", readLines(f)), g)
  expect_equal(names(read_enrichment_csv(g)), names(df))
})

test_that("version suffix is added before the file extension", {
  expect_equal(with_version("Heatmap_X.png"), paste0("Heatmap_X_v", app_version, ".png"))
  expect_equal(with_version("Filtered.zip"),  paste0("Filtered_v", app_version, ".zip"))
})

test_that("chord labels stay distinct when truncation makes them identical", {
  x   <- c("POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION", "POSITIVE_REGULATION_OF_CELL_ADHESION",
           "CYTOKINE_PRODUCTION")
  plain <- shorten_labels(x, 22, "truncate")
  expect_equal(plain[1], plain[2])                         # the problem: identical labels
  lab <- unique_pathway_labels(x, 22, "truncate")
  expect_equal(length(unique(lab)), 3)
  expect_match(lab[1], "PRODUCTION$"); expect_match(lab[2], "ADHESION$")
  expect_equal(lab[3], "CYTOKINE_PRODUCTION")              # short names are unchanged
})
