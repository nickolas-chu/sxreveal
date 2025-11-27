#' quickAdd: Quick Sex Assignment for Seurat
#'
#' Used to add Xist expression, summed Y chromosome gene expression,
#' probabilities being female, and sex assignments to the metadata of a
#' Seurat object
#'
#' @param Probabilities dataframe of probabilities
#' @param Seuratobj A target Seurat object
#' @export
#' @author Nickolas C. Chu
quickAdd <- function(Probabilities, Seuratobj ){
  # Ensure cell_id column exists
  if (!"cell_id" %in% colnames(Probabilities)) {
    stop("Probabilities dataframe must contain a 'cell_id' column")
  }
  cols_to_add <- setdiff(colnames(Probabilities), c("cell_id", "cluster"))
  
  # Align by cell_id without reordering Seurat metadata
  meta <- Seuratobj@meta.data
  # Match rows by cell_id
  idx <- match(rownames(meta), Probabilities$cell_id)
  # For each column, add values into metadata
  for (col in cols_to_add) {
    meta[[col]] <- Probabilities[[col]][idx]
  }
  # Replace metadata
  Seuratobj@meta.data <- meta
  
  return(Seuratobj)

}
