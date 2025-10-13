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
quickAdd <- function(Probabilities, Seuratobj )
{
  Seuratobj[["Ygenes2"]] <- Probabilities$Ygenes
  Seuratobj[["Xist2"]] <- Probabilities$Xist
  Seuratobj[["ProbUni"]] <- Probabilities$ProbFemaleUni
  Seuratobj[["SexUni"]] <- Probabilities$SexUni
  Seuratobj[["ProbMulti"]] <- Probabilities$ProbFemaleMulti
  Seuratobj[["SexMulti"]] <- Probabilities$SexMulti
  Seuratobj[["ProbMultinCount"]] <- Probabilities$ProbFemaleMultinCount
  Seuratobj[["SexMultinCount"]] <- Probabilities$SexMultinCount
  Seuratobj[["ProbMultiXY"]] <- Probabilities$ProbFemaleXY
  Seuratobj[["SexMultiXY"]] <- Probabilities$SexMultiXY

  return(Seuratobj)

}
