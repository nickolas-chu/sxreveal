#' addProbs: Sex assigning function
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
addProbs <- function(Probabilities, femalecut = 0.8, malecut = 0.2)
{
  for (i in 1:nrow(Probabilities)) {
    if (Probabilities$ProbFemaleMulti[i] >= 0.8) {
      Probabilities[['SexMulti']][i] <- "female"
    }else if(Probabilities$ProbFemaleMulti[i] <= 0.2) {
      Probabilities[['SexMulti']][i] <- "male"
    }else {
      Probabilities[['SexMulti']][i] <- "soup"
    }
  }
  #uni
  for (i in 1:nrow(Probabilities)) {
    if (Probabilities$ProbFemaleUni[i] >= 0.8) {
      Probabilities[['SexUni']][i] <- "female"
    }else if(Probabilities$ProbFemaleUni[i] <= 0.2) {
      Probabilities[['SexUni']][i] <- "male"
    }else {
      Probabilities[['SexUni']][i] <- "soup"
    }
  }
  #multincount
  for (i in 1:nrow(Probabilities)) {
    if (Probabilities$ProbFemaleMultinCount[i] >= 0.8) {
      Probabilities[['SexMultinCount']][i] <- "female"
    }else if(Probabilities$ProbFemaleMultinCount[i] <= 0.2) {
      Probabilities[['SexMultinCount']][i] <- "male"
    }else {
      Probabilities[['SexMultinCount']][i] <- "soup"
    }
  }
  return(Probabilities)

}
