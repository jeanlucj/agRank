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
#' @param method one of "BT", "PL", "TH", "MPM", "LM"
#' @param rate the learning rate
#' @param maxiter the maximum number of iterations
#' @param tol the tolerance
#' @param start initial guesses at scores
#' @param decay how fast the learning rate decays when log post doesn't improve
#'
#' @return Return a list with two components:
#'     \item{ranks}{a vector where the i-th element is the rank assigned to the i-th item.}
#'     \item{ranking}{a vector where the i-th element is the item ranked in the i-th place}
#'
#'
#'
#' @examples
#' #synthetic ranking data
#' data = rbind(c(3,1,2), c(3,1,2), c(3,2,1))
#' #use identity relationship matrix
#' K = diag(1, 3)
#'
#' #rank aggregation
#' rankAg(data, K, method = 'TH')
#'
#'
#'
#'
#' @export

rankAg = function(data, K = diag(ncol(data)), method = "TH"){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0


  nVarieties = ncol(data)


  mu = rep(1, nVarieties) #mean vector of the normal prior on scores
  rate = 1
  maxiter = 5000
  tol = 1e-8
  #starting point for parameters, the first nvar elements are for scores,
  #the next nobs elements are for adherences
  start = rnorm(ncol(data),10,1)
  decay = 1.1
  param=start
  sigma <- diag(nVarieties)

  if(method == 'BT'){
    if(!is.matrix(K)){
      stop('relationship matrix must be specified for BT model')
    }
    score = sgdBT(data, mu, K, rate, maxiter, tol, start, decay)$scores
    names(score) = 1:nVarieties #assign labels
    ranking = as.numeric(names(sort(score, decreasing = T)))
    ranks = match(1:nVarieties, ranking)

  } else if(method == 'PL'){
    if(!is.matrix(K)){
      stop('relationship matrix must be specified for PL model')
    }

    score = sgdPL(data, mu, K, rate, maxiter, tol, start, decay)$score
    names(score) = 1:nVarieties #assign labels
    ranking = as.numeric(names(sort(score, decreasing = T)))
    ranks = match(1:nVarieties, ranking)

  } else if(method == 'TH'){

    if(!is.matrix(K)){
      stop('relationship matrix must be specified for TH model')
    }

    score = sgdThurs(data, mu, K, rate, maxiter, tol, start, decay)$scores
    names(score) = 1:nVarieties #assign labels
    ranking = as.numeric(names(sort(score, decreasing = T)))
    ranks = match(1:nVarieties, ranking)

  } else if(method == 'MPM'){

    if(!is.matrix(K)){
      stop('relationship matrix must be specified for MPM model')
    }
    start = c(start, rep(1, nVarieties))
    score = sgdMPM(data, mu, K, rate, maxiter, tol, start, decay)$score
    names(score) = 1:nVarieties #assign labels
    ranking = as.numeric(names(sort(score, decreasing = T)))
    ranks = match(1:nVarieties, ranking)


  } else if(method == 'LM'){

    if(!is.matrix(K)){

      res = rankLM(data)
      ranks = res$ranks
      ranking = res$ranking

    } else{

      res = rankLM(data, K)
      ranks = res$ranks
      ranking = res$ranking

    }

  } else{

    print('method must be one of BT, PL, TH, MPM, LM')

  }


  return(list(ranks = ranks, ranking = ranking))

}

