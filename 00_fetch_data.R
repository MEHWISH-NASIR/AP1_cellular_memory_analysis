dir.create("data", showWarnings = FALSE)


ap1_base_url <- paste0(
  "https://raw.githubusercontent.com/",
  "arjunrajlaboratory/cellularmemory/main/",
  "extractedData/burninrnaseq/202208/"
)

urls <- c(
  "log2CPM_matrix.txt" = paste0(ap1_base_url, "log2CPM_matrix.txt"),
  "studydesign.tsv"    = paste0(ap1_base_url, "studydesign.tsv"),
  "dex_log2FC.txt"     = paste0(ap1_base_url, "dex_log2FC.txt"),
  "dextable.txt"       = paste0(ap1_base_url, "dextable.txt"),

  # Human kinase reference list
  "KinHubKinaseList.csv" =
    "https://raw.githubusercontent.com/openkinome/kinodata/master/data/KinHubKinaseList.csv"
)

# ------------------------------------------------------------
# Download files only when absent
# ------------------------------------------------------------

for (f in names(urls)) {

  destination <- file.path("data", f)

  if (!file.exists(destination)) {

    message("Downloading: ", f)

    temp_file <- paste0(destination, ".tmp")

    download.file(
      urls[[f]],
      destfile = temp_file,
      mode = "wb",
      quiet = FALSE
    )

    if (!file.exists(temp_file) ||
        file.info(temp_file)$size == 0) {

      unlink(temp_file)

      stop(
        "Download failed or produced an empty file: ",
        f
      )
    }

    if (!file.rename(temp_file, destination)) {
      unlink(temp_file)

      stop(
        "Could not move downloaded file to: ",
        destination
      )
    }

    message("Downloaded successfully: ", f)

  } else {

    message("Already exists: ", f)
  }
}

message("Input data setup complete.")
