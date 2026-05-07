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

log_message("Installing readxl from CRAN")
install.packages("readxl", repos = cran_repo, dependencies = TRUE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("readxl installation failed. Package still unavailable.", call. = FALSE)
}

log_message("readxl installed successfully")
log_message("readxl version: ", as.character(utils::packageVersion("readxl")))
