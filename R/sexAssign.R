#' sexAssign: Sex assigning function
#'
#' Used to assign sex to each cell based on probability of being female
#' Default cutoffs:
#' probability >= 0.8  -> female
#' probability <= 0.2 -> male
#' else -> soup
#'
#' @param Probabilities dataframe
#' @param femalecut selected probability of being female used to assign females
#' @param malecut selected probability of being female used to assign males
#' @export
#' @author Nickolas C. Chu
sexAssign <- function(Probabilities, femalecut = 0.8, malecut = 0.2) {
  for (i in 1:nrow(Probabilities)) {
    val <- Probabilities$ProbFemaleXY[i]
    if (!is.na(val) && val >= femalecut) {
      Probabilities[["SexMultiXY"]][i] <- "female"
    } else if (!is.na(val) && val <= malecut) {
      Probabilities[["SexMultiXY"]][i] <- "male"
    } else {
      Probabilities[["SexMultiXY"]][i] <- "soup"
    }
  }
  for (i in 1:nrow(Probabilities)) {
    val <- Probabilities$ProbFemaleUniY[i]
    if (!is.na(val) && val >= femalecut) {
      Probabilities[["SexUniY"]][i] <- "female"
    } else if (!is.na(val) && val <= malecut) {
      Probabilities[["SexUniY"]][i] <- "male"
    } else {
      Probabilities[["SexUniY"]][i] <- "soup"
    }
  }
  for (i in 1:nrow(Probabilities)) {
    val <- Probabilities$ProbFemaleMulti[i]
    if (!is.na(val) && val >= femalecut) {
      Probabilities[['SexMulti']][i] <- "female"
    } else if (!is.na(val) && val <= malecut) {
      Probabilities[['SexMulti']][i] <- "male"
    } else {
      Probabilities[['SexMulti']][i] <- "soup"
    }
  }
  for (i in 1:nrow(Probabilities)) {
    val <- Probabilities$ProbFemaleUni[i]
    if (!is.na(val) && val >= femalecut) {
      Probabilities[['SexUniX']][i] <- "female"
    } else if (!is.na(val) && val <= malecut) {
      Probabilities[['SexUniX']][i] <- "male"
    } else {
      Probabilities[['SexUniX']][i] <- "soup"
    }
  }
  for (i in 1:nrow(Probabilities)) {
    val <- Probabilities$ProbFemaleMultinCount[i]
    if (!is.na(val) && val >= femalecut) {
      Probabilities[['SexMultinCount']][i] <- "female"
    } else if (!is.na(val) && val <= malecut) {
      Probabilities[['SexMultinCount']][i] <- "male"
    } else {
      Probabilities[['SexMultinCount']][i] <- "soup"
    }
  }
  return(Probabilities)
}

