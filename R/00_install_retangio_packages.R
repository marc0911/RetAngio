options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = 8)

lib_user <- Sys.getenv("R_LIBS_USER")
if (!nzchar(lib_user)) {
  stop("R_LIBS_USER is not set.")
}
dir.create(lib_user, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_user, .libPaths()))

message("Using .libPaths():")
print(.libPaths())

install_if_missing <- function(pkgs) {
  ip <- rownames(installed.packages())
  to_install <- setdiff(pkgs, ip)
  if (length(to_install) > 0) {
    message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install)
  } else {
    message("All requested CRAN packages already installed.")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install_bioc_if_missing <- function(pkgs) {
  ip <- rownames(installed.packages())
  to_install <- setdiff(pkgs, ip)
  if (length(to_install) > 0) {
    message("Installing Bioconductor packages: ", paste(to_install, collapse = ", "))
    BiocManager::install(to_install, ask = FALSE, update = FALSE)
  } else {
    message("All requested Bioconductor packages already installed.")
  }
}

cran_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "yaml",
  "Matrix",
  "hdf5r",
  "ggplot2",
  "patchwork",
  "cowplot",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "purrr",
  "stringr",
  "forcats",
  "glue",
  "scales",
  "viridis",
  "future",
  "future.apply",
  "clustree",
  "harmony",
  "janitor"
)

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scDblFinder",
  "scater",
  "scran",
  "scuttle",
  "batchelor"
)

install_if_missing(cran_pkgs)
install_bioc_if_missing(bioc_pkgs)

message("Final package check:")
check_pkgs <- c(
  "Seurat", "SeuratObject", "yaml", "Matrix", "hdf5r",
  "SingleCellExperiment", "scDblFinder", "scater", "scran", "scuttle",
  "batchelor", "clustree", "harmony", "future", "future.apply"
)

status <- vapply(check_pkgs, requireNamespace, logical(1), quietly = TRUE)
print(status)

if (!all(status)) {
  missing_pkgs <- names(status)[!status]
  stop("Some packages are still missing: ", paste(missing_pkgs, collapse = ", "))
}

message("RetAngio package bootstrap completed successfully.")
