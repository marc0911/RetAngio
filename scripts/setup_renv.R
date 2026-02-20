if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

if (file.exists("renv.lock")) {
  renv::restore(prompt = FALSE)
} else {
  renv::init(bare = TRUE)
  renv::install(c(
    "Seurat", "SeuratObject", "tidyverse", "yaml", "fs", "digest", "readxl", "janitor",
    "edgeR", "DESeq2", "fgsea", "clusterProfiler", "scDblFinder", "SoupX", "clustree"
  ))
  renv::snapshot(prompt = FALSE)
}
