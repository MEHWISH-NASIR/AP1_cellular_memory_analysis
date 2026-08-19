# render_report.R
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop(
    "The rmarkdown package is required. ",
    "Run renv::restore() first."
  )
}

if (!rmarkdown::pandoc_available()) {
  stop(
    "Pandoc is not available."
  )
}

if (!file.exists("analysis_report.Rmd")) {
  stop(
    "analysis_report.Rmd was not found."
  )
}

dir.create(
  "docs",
  showWarnings = FALSE,
  recursive = TRUE
)

rmarkdown::render(
  input = "analysis_report.Rmd",
  output_file = "index.html",
  output_dir = "docs",
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)

output_file <- "docs/index.html"

if (!file.exists(output_file)) {
  stop(
    "Rendering completed but ",
    output_file,
    " was not created."
  )
}

cat(
  "\nHTML report created successfully:\n",
  output_file,
  "\n",
  sep = ""
)