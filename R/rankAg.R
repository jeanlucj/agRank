#' Ranking aggregation of triple comparisons
#'
#' Ranking aggregation using Bradley-Terry model, Plackett-Luce model, Thurstone model,
#' and linear model.
#' This is the wrapper function for targetBT, targetPL, targetTH, and rankLM
#'
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items to rank;
#' each row vector is a partial ranking (triple comparisons),
#' with i-th element being the rank assigned to item i;
#' if an item is not ranked its partial ranking is NA
#' @param K the additive relationship matrix;
#' all methods must specify K except LM
#' @param method one of "BT", "PL", "TH", "LM"
#' The following parameters are only relevant for the ranking models
#' @param rate the learning rate
#' @param maxiter the maximum number of iterations
#' @param tol the tolerance
#' @param startVar initial guess at variance of scores
#' @param startScores initial guesses at scores
#' @param decay how fast the learning rate decays when log post doesn't improve
#'
#' @return Return a list with three components:
#'     \item{ranks}{a vector where the i-th element is the rank assigned to the i-th item.}
#'     \item{scores}{a vector of the estimated scores of the varieties.}
#'     \item{scoreVar}{a scalar of the estimated score variance. This variance is method specific.}
#'
#' @examples
#' #synthetic ranking data
#' data <- rbind(c(0,3,1,2), c(3,0,1,2), c(3,2,0,1), c(2,3,1,0), c(0,3,2,1))
#' #use identity relationship matrix
#' K <- diag(1, 4)
#'
#' #rank aggregation
#' rankAg(data, K, method='TH')
#'
#' @export

rankAg <- function(data, K=diag(ncol(data)), method="TH", rate=0.01, maxiter=1000, tol=1e-8, startVar=1, startScores=NULL, decay=0.9){
  if (method %in% c("BT", "PL", "TH")){
    targetFN <- switch(method, BT=targetBT, PL=targetPL, TH=targetTH)
    res <- gradientDescent(data, targetFN, K, rate, maxiter, tol, startVar, startScores, decay)
  } else{ # Use linear model
    res <- rankLM(data, K)
  }
  scores <- res$scores
  scoreVar <- res$scoreVar
  ranks <- rank(scores)

  return(list(ranks=ranks, scores=scores, scoreVar=scoreVar))
}
