# Goal: Start the app with a small synthetic dataset and check that every view is drawn and that
#       the status lines report the expected settings.
# Why:  Guards the behavior described in the manuscript (log2 option, warning for log data, group
#       order, version in status lines and file names) against regressions.
library(shinytest2)

# 30 genes x 8 samples of linear counts; group "Ctrl" listed first in the annotation
make_inputs <- function(log_scale = FALSE) {
  set.seed(1)
  genes <- paste0("G", 1:30)
  m <- matrix(rpois(30 * 8, 200), nrow = 30, dimnames = list(genes, paste0("S", 1:8)))
  m[1:10, 5:8] <- m[1:10, 5:8] * 3                      # G1-G10 higher in "Case"
  if (log_scale) m <- log2(m + 1)
  norm <- data.frame(gene_symbol = genes, m, check.names = FALSE)
  ann  <- data.frame(sample = paste0("S", 1:8), group = rep(c("Ctrl", "Case"), each = 4))
  path <- c("POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY", "POSITIVE_REGULATION_OF_SYNAPTIC_PLASTICITY",
            "CYTOKINE_PRODUCTION")
  enr  <- data.frame(Description = path, NES = c(2, 1.5, -1.8), p.adjust = c(0.001, 0.001, 0.01),
                     core_enrichment = c("G1/G2/G3/G4/G5", "G1/G2/G3/G6", "G4/G5/G7/G8"))
  de   <- data.frame(gene_symbol = genes, log2FoldChange = c(rep(1.6, 10), rep(0, 20)),
                     padj = c(rep(0.001, 10), rep(0.9, 20)))
  list(norm = tmp_csv(norm), ann = tmp_csv(ann), enr = tmp_csv(enr), de = tmp_csv(de))
}

start_app <- function(f) {
  app <- AppDriver$new(normalizePath(file.path(test_path(), "..", "..")), name = "enrichviz",
                       height = 1000, width = 1400, load_timeout = 60 * 1000, timeout = 60 * 1000)
  app$upload_file(file_norm_prot = f$norm); app$upload_file(file_annotation = f$ann)
  app$upload_file(file_ipa_funct = f$enr);  app$upload_file(file_de = f$de)
  app$wait_for_idle(1000)
  app
}
rendered <- function(app, id) !is.null(app$wait_for_value(output = id, ignore = list(NULL, ""))$src)

test_that("all views are drawn and status lines report version, transformation and ties", {
  skip_on_cran()
  app <- start_app(make_inputs())
  on.exit(app$stop(), add = TRUE)
  tabs <- c("Bar Plot" = "barplot", "Chord Diagram" = "chord_plot", "Molecule Frequency" = "gene_freq_plot",
            "Heatmap" = "heatmap_plot", "Boxplot" = "boxplot_plot", "Volcano" = "volcano_plot",
            "PCA" = "pca_plot", "DE Heatmap" = "de_heatmap_plot")
  app$set_inputs(freq_min_pathways = 1, wait_ = FALSE)
  for (tb in names(tabs)) {
    app$set_inputs(main_tabs = tb, wait_ = FALSE)
    expect_true(rendered(app, tabs[[tb]]), label = tb)
  }
  st <- app$get_value(output = "status")
  expect_match(st, paste0("EnrichViz        : v", app_version), fixed = TRUE)
  expect_match(st, "Transformation   : log2(x + 1), then row Z-score", fixed = TRUE)
  expect_false(grepl("WARNING", st))
  expect_match(app$get_value(output = "boxplot_status"), "Transformation   : log2(x + 1)", fixed = TRUE)
  expect_match(app$get_value(output = "barplot_status"), "Tied pathways", fixed = TRUE)
  expect_match(app$get_value(output = "barplot_status"), paste0("v", app_version), fixed = TRUE)
  expect_match(basename(app$get_download("download_volcano_main")), paste0("_v", app_version, ".png"), fixed = TRUE)
})

test_that("log-scale data with the log2 option on give a warning", {
  skip_on_cran()
  app <- start_app(make_inputs(log_scale = TRUE))
  on.exit(app$stop(), add = TRUE)
  app$set_inputs(main_tabs = "Heatmap", wait_ = FALSE)
  expect_true(rendered(app, "heatmap_plot"))
  expect_match(app$get_value(output = "status"), "WARNING: data look already log-transformed", fixed = TRUE)
  app$set_inputs(de_log_transform = FALSE)
  expect_match(app$get_value(output = "status"), "Transformation   : none", fixed = TRUE)
})
