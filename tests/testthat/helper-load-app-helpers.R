# Goal: Load the top-level helper functions of app.R (column detection, CSV reading) without
#       starting the app.
# Why:  The helpers do not use Shiny inputs, so they can be tested directly; ui, server and
#       shinyApp() are skipped.
app_file <- normalizePath(file.path(testthat::test_path(), "..", "..", "app.R"))
for (ex in parse(app_file)) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) &&
      !as.character(ex[[2]]) %in% c("ui", "server"))
    eval(ex, envir = globalenv())
}

# Writes a data frame to a temporary CSV and returns the path
tmp_csv <- function(df) {
  f <- tempfile(fileext = ".csv")
  write.csv(df, f, row.names = FALSE)
  f
}
