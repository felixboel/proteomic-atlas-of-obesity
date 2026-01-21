get_project_root <- function() {
  normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
}

load_processed_data <- function(data_dir = NULL) {
  if (is.null(data_dir)) {
    data_dir <- file.path(get_project_root(), "data", "processed")
  }

  paths <- list(
    protein_map = file.path(data_dir, "uniprot2reactome.rds"),
    pathway_list = file.path(data_dir, "reactome_pathways.rds"),
    pathway_hierarchy = file.path(data_dir, "reactome_pathways_relation.rds"),
    diff_expr = file.path(data_dir, "combined_all.rds")
  )

  missing <- names(paths)[!file.exists(unlist(paths))]
  if (length(missing) > 0) {
    stop(
      paste0(
        "Missing processed data files: ",
        paste(missing, collapse = ", "),
        ". Run scripts/prepare_data.R to generate them."
      ),
      call. = FALSE
    )
  }

  list(
    protein_map = readRDS(paths$protein_map),
    pathway_list = readRDS(paths$pathway_list),
    pathway_hierarchy = readRDS(paths$pathway_hierarchy),
    diff_expr = readRDS(paths$diff_expr)
  )
}
