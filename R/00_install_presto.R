#!/usr/bin/env Rscript

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste(..., collapse = "")))
}

user_lib <- Sys.getenv("R_LIBS_USER")
if (identical(user_lib, "")) {
  stop("R_LIBS_USER is empty.", call. = FALSE)
}

dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

log_message("R version: ", R.version.string)
log_message("R_LIBS_USER: ", user_lib)
log_message("Effective .libPaths(): ", paste(.libPaths(), collapse = " | "))

cran_repo <- "https://cloud.r-project.org"

cran_needed <- c("remotes", "Rcpp", "RcppArmadillo", "Matrix")
cran_missing <- cran_needed[!vapply(cran_needed, requireNamespace, logical(1), quietly = TRUE)]

if (length(cran_missing) > 0) {
  log_message("Installing missing CRAN packages: ", paste(cran_missing, collapse = ", "))
  install.packages(cran_missing, repos = cran_repo, dependencies = TRUE)
} else {
  log_message("Required CRAN bootstrap packages already available")
}

log_message("Installing presto from GitHub")
remotes::install_github(
  repo = "immunogenomics/presto",
  upgrade = "never",
  dependencies = TRUE
)

if (!requireNamespace("presto", quietly = TRUE)) {
  stop("presto installation failed. Package still unavailable after install_github.", call. = FALSE)
}

log_message("presto installed successfully")
log_message("presto version: ", as.character(utils::packageVersion("presto")))

sessionInfo()
