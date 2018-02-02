#' Ranking aggregation of triple comparisons
#'
#'
#' Ranking aggregation using Bradley-Terry model, Plackett-Luce model, Thurstone model, Multinomial Preference model
#' and linear model. This is the wrapper function for rankLM, sgdBT, sgdPL, sgdThurs and sgdMPM.
#'
#'
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items to rank;
#' each row vector is a partial ranking (triple comparisons),
#' with i-th element being the rank assined to item i;
#' the entry where that item is not ranked in the partial ranking is replaced by 0
#' @param K the additive relationship matrix;
#' all methods must specify K except LM
#' @param method one of "BT", "PL", "TH", "LM"
#' @param rate the learning rate
#' @param maxiter the maximum number of iterations
#' @param tol the tolerance
#' @param startVar initial guess at variance of scores
#' @param startScores initial guesses at scores
#' @param decay how fast the learning rate decays when log post doesn't improve
#' @param nTechRep for methods that may not converge perfectly, try more than once and pick the best
#'
#' @return Return a list with two components:
#'     \item{ranks}{a vector where the i-th element is the rank assigned to the i-th item.}
#'     \item{ranking}{a vector where the i-th element is the item ranked in the i-th place}
#'
#' @examples
#' #synthetic ranking data
#' data <- rbind(c(3,1,2), c(3,1,2), c(3,2,1))
#' #use identity relationship matrix
#' K <- diag(1, 3)
#'
#' #rank aggregation
#' rankAg(data, K, method = 'TH')
#'
#' @export

rankAg <- function(data, K = diag(ncol(data)), method = "TH", rate=0.01, maxiter=5000, tol=1e-8, startVar=1, startScores=rnorm(ncol(data)), decay=1.1, nTechRep=2){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  nVarieties <- ncol(data)

  sigma <- diag(nVarieties)

  if (method %in% c("BT", "PL", "TH")){
    if(!is.matrix(K)){
      stop('relationship matrix must be specified for BT model')
    }
    
    value <- Inf
    for (techRep in 1:nTechRep){
      fails <- -1
      doOver <- TRUE
      while (doOver){
        fails <- fails + 1
        res <- try(switch(method,
          BT = {
            sgdBT(data, K, rate, maxiter, tol, startVar, startScores, decay)
          },
          PL = {
            sgdPL(data, K, rate, maxiter, tol, startVar, startScores, decay)
          },
          TH = {
            sgdThurs(data, K, rate, maxiter, tol, startVar, startScores, decay)
          }
        ), silent=T)#END switch
        doOver <- class(res) == "try-error"
      }
      if(res$value < value){
        scores <- res$scores
        value <- res$value
      } 
    }#END techRep
    
    names(scores) <- 1:nVarieties #assign labels
    ranking <- as.numeric(names(sort(scores, decreasing = T)))
    ranks <- match(1:nVarieties, ranking)
  } else { # Use linear model
    if(!is.matrix(K)){
      res <- rankLM(data)
      ranks <- res$ranks
      ranking <- res$ranking
    } else{
      res <- rankLM(data, K)
      ranks <- res$ranks
      ranking <- res$ranking
    }
  }

  return(list(ranks = ranks, ranking = ranking))
}

